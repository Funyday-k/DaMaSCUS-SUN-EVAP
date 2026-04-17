import os
import glob
import time
import numpy as np
from multiprocessing import Pool, Value

# ============================================================
# 常量
# ============================================================
rSun = 6.957E5  # km
bins = np.linspace(0, 2.0 * rSun, 2001)
NUM_BINS = len(bins) - 1
N_COLS = 6  # 格式: time[s], r[km], vx[km/s], vy[km/s], vz[km/s], E[eV]

DTYPE = np.float32


# ============================================================
# 核心：单批次处理函数 —— 一次读取同时计算 histogram 和 v²
# ============================================================
# 全局共享计数器（在 worker 初始化时设置）
_shared_counter = None
_shared_start_time = None
_total_files = None

def _init_worker(counter, start_t, total_n):
    global _shared_counter, _shared_start_time, _total_files
    _shared_counter = counter
    _shared_start_time = start_t
    _total_files = total_n


def process_batch(file_list):
    """处理一批文件，在 worker 内部累加，返回 (hist, v2_tsum, t_sum, captured_count, total_count)。
    只统计被捕获的轨迹（能量 E <= 0 出现过）。
    减少 IPC 开销：每个 worker 只返回一次结果。"""
    global _shared_counter, _shared_start_time, _total_files
    hist    = np.zeros(NUM_BINS, dtype=np.float64)
    v2_tsum = np.zeros(NUM_BINS, dtype=np.float64)
    t_sum   = np.zeros(NUM_BINS, dtype=np.float64)
    captured_count = 0
    total_count = 0

    for file_path in file_list:
        try:
            raw = np.fromfile(file_path, dtype=DTYPE)
            if raw.size == 0 or raw.size % N_COLS != 0:
                continue
            mat = raw.reshape(-1, N_COLS)
            del raw
            if len(mat) < 2:
                continue

            total_count += 1

            # 检测是否被捕获：能量列 (index 5) 是否出现 <= 0
            energy = mat[:, 5]
            if not np.any(energy <= 0):
                # 未被捕获，跳过此轨迹
                del mat
                continue

            captured_count += 1

            # 新格式：col0=time, col1=r, col2-4=vx/vy/vz, col5=E
            t  = mat[:, 0].copy()
            r  = mat[:, 1].copy()
            v2 = mat[:, 2]**2 + mat[:, 3]**2 + mat[:, 4]**2
            del mat

            # 梯形时间权重
            dt = np.empty_like(t)
            dt[0]    = (t[1] - t[0]) * 0.5
            dt[-1]   = (t[-1] - t[-2]) * 0.5
            dt[1:-1] = (t[2:] - t[:-2]) * 0.5
            del t

            # bin 索引（searchsorted 比 digitize 快）
            idx = np.searchsorted(bins, r, side='right') - 1
            mask = (idx >= 0) & (idx < NUM_BINS)
            idx_m  = idx[mask]
            dt_m   = dt[mask]
            v2dt_m = v2[mask] * dt_m
            del r, v2, dt, idx, mask

            # bincount 比 np.add.at 快 5-10 倍
            # hist 和 t_sum 完全相同（都是 Σdt per bin），只算一次
            dt_bc   = np.bincount(idx_m, weights=dt_m,   minlength=NUM_BINS)[:NUM_BINS]
            v2dt_bc = np.bincount(idx_m, weights=v2dt_m, minlength=NUM_BINS)[:NUM_BINS]
            del idx_m, dt_m, v2dt_m

            hist    += dt_bc
            t_sum   += dt_bc      # 与 hist 相同
            v2_tsum += v2dt_bc

        except Exception as e:
            print(f"Error reading {file_path}: {e}", flush=True)

        # 更新共享进度计数器
        if _shared_counter is not None:
            with _shared_counter.get_lock():
                _shared_counter.value += 1
                cnt = _shared_counter.value
            if cnt % 10000 == 0:
                elapsed = time.time() - _shared_start_time.value
                rate = cnt / elapsed if elapsed > 0 else 0
                total = _total_files.value
                eta = (total - cnt) / rate if rate > 0 else float('inf')
                print(f"  [进度] {cnt}/{total} 文件 "
                      f"({100*cnt/total:.1f}%), "
                      f"速率: {rate:.0f} files/s, "
                      f"ETA: {eta/60:.1f} min", flush=True)

    return hist, v2_tsum, t_sum, captured_count, total_count


# ============================================================
# 配置
# ============================================================
input_dir = "/lustre/scratch/kennyng/DaMaSCUS_OUT/results_0.000000_-40.000000"
output_base = "/project/kennyng/backup_DM/bincount"
# !! 关键修改：减少并行进程数以降低 Lustre MDS 压力 !!
# Lustre 元数据服务器不适合大量并发小文件 I/O
# 16 进程 → 8 进程可大幅减少 MDS 拥堵
num_processes = 8

if __name__ == '__main__':
    # 获取文件列表（os.listdir 比 glob 在大目录上快很多）
    print(f"正在扫描目录 {input_dir} ...", flush=True)
    t0 = time.time()
    all_entries = os.listdir(input_dir)
    txt_files = sorted([os.path.join(input_dir, f) for f in all_entries if f.endswith('.dat')])
    del all_entries
    num_files = len(txt_files)
    print(f"Found {num_files} trajectory files to process in {input_dir} "
          f"(扫描耗时 {time.time()-t0:.1f}s)", flush=True)

    input_dir_name = os.path.basename(os.path.normpath(input_dir))

    # ============================================================
    # 先做小规模基准测试，估算总耗时
    # ============================================================
    BENCHMARK_N = min(200, num_files)
    print(f"\n===== 基准测试：处理 {BENCHMARK_N} 个文件 =====", flush=True)
    t_bench_start = time.time()
    bench_result = process_batch(txt_files[:BENCHMARK_N])
    t_bench = time.time() - t_bench_start
    per_file = t_bench / BENCHMARK_N if BENCHMARK_N > 0 else 0
    est_serial = per_file * num_files
    est_parallel = est_serial / num_processes
    print(f"  基准: {BENCHMARK_N} 文件耗时 {t_bench:.2f}s, "
          f"每文件 {per_file*1000:.1f}ms", flush=True)
    print(f"  估算: 串行 {est_serial/60:.0f} min, "
          f"{num_processes} 进程并行 ~{est_parallel/60:.0f} min", flush=True)

    # ============================================================
    # 将文件列表均分为批次，每个 worker 处理一批
    # 使用更多更小的批次，让进度更新更频繁
    # ============================================================
    batch_size = max(1, (num_files + num_processes * 4 - 1) // (num_processes * 4))
    batches = [txt_files[i:i + batch_size] for i in range(0, num_files, batch_size)]
    num_batches = len(batches)
    print(f"\n分为 {num_batches} 个批次，每批约 {batch_size} 个文件", flush=True)

    # ============================================================
    # 单次并行：同时计算 histogram 和 v²
    # ============================================================
    print("\n===== 并行计算 histogram + v² (单次遍历) =====", flush=True)
    total_hist    = np.zeros(NUM_BINS, dtype=np.float64)
    total_v2_tsum = np.zeros(NUM_BINS, dtype=np.float64)
    total_t_sum   = np.zeros(NUM_BINS, dtype=np.float64)
    total_captured = 0
    total_valid    = 0

    # 共享进度计数器
    shared_counter = Value('i', 0)
    shared_start = Value('d', time.time())
    shared_total = Value('i', num_files)

    t_start = time.time()
    with Pool(processes=num_processes,
              initializer=_init_worker,
              initargs=(shared_counter, shared_start, shared_total)) as pool:
        result_iter = pool.imap_unordered(process_batch, batches)
        for i in range(num_batches):
            try:
                # 不设硬超时 —— 之前的超时导致作业被错误终止
                # Lustre I/O 本身就慢，设超时只会浪费 SLURM 机时
                hist, v2_tsum, t_sum, captured_count, batch_total = result_iter.next()
            except StopIteration:
                print(f"\n[警告] 迭代器提前结束 (batch {i+1})", flush=True)
                break
            except Exception as e:
                print(f"\n[错误] 批次 {i+1} 异常: {type(e).__name__}: {e}", flush=True)
                continue  # 跳过失败批次，继续处理其余
            total_hist    += hist
            total_v2_tsum += v2_tsum
            total_t_sum   += t_sum
            total_captured += captured_count
            total_valid    += batch_total
            elapsed = time.time() - t_start
            processed = shared_counter.value
            rate = processed / elapsed if elapsed > 0 else 0
            capture_ratio = total_captured / total_valid if total_valid > 0 else 0
            print(f"  批次 {i+1}/{num_batches} 完成，"
                  f"累计 {processed}/{num_files} 文件，"
                  f"有效: {total_valid}，捕获: {total_captured} ({capture_ratio:.4%})，"
                  f"速率: {rate:.0f} files/s，"
                  f"已用: {elapsed/60:.1f}min",
                  flush=True)

    capture_ratio_final = total_captured / total_valid if total_valid > 0 else 0
    print(f"\n全部完成，有效文件: {total_valid}/{num_files}，"
          f"捕获: {total_captured} ({capture_ratio_final:.6%})", flush=True)

    # ============================================================
    # 输出 histogram
    # ============================================================
    total_output_file = os.path.join(output_base, f"total_bincount_{input_dir_name}.txt")
    if os.path.exists(total_output_file):
        print(f"[覆盖] 检测到已有文件 {total_output_file}，将覆盖", flush=True)
    with open(total_output_file, 'w') as f:
        f.write(f"# Ntot(captured) = {total_captured}, Ntot(valid) = {total_valid}, "
                f"capture_ratio = {capture_ratio_final:.8e}\n")
        for i in range(NUM_BINS):
            f.write(f"{i}\t{int(total_hist[i])}\n")
    print(f"Histogram 结果已保存到 {total_output_file}", flush=True)

    # ============================================================
    # 输出 capture ratio 单独文件
    # ============================================================
    ratio_output_file = os.path.join(output_base, f"capture_ratio_{input_dir_name}.txt")
    with open(ratio_output_file, 'w') as f:
        f.write(f"# Capture ratio for {input_dir_name}\n")
        f.write(f"total_valid_files\t{total_valid}\n")
        f.write(f"captured_files\t{total_captured}\n")
        f.write(f"capture_ratio\t{capture_ratio_final:.8e}\n")
    print(f"Capture ratio 已保存到 {ratio_output_file}", flush=True)

    # ============================================================
    # 输出 v²
    # ============================================================
    mean_v2 = np.zeros(NUM_BINS, dtype=np.float64)
    nonzero = total_t_sum > 0
    mean_v2[nonzero] = total_v2_tsum[nonzero] / total_t_sum[nonzero]

    v2_output_file = os.path.join(output_base, f"total_v2_{input_dir_name}.txt")
    if os.path.exists(v2_output_file):
        print(f"[覆盖] 检测到已有文件 {v2_output_file}，将覆盖", flush=True)
    with open(v2_output_file, 'w') as f:
        f.write("# bin_index\tv2_sum\t\t\tcount\t\tmean_v2\n")
        for i in range(NUM_BINS):
            count = int(total_t_sum[i]) if total_t_sum[i] > 0 else 0
            f.write(f"{i}\t{total_v2_tsum[i]:.6e}\t{count}\t{mean_v2[i]:.6e}\n")
    print(f"v^2 统计结果已保存到 {v2_output_file}", flush=True)

print("全部完成!", flush=True)
