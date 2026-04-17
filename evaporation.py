import os
import numpy as np
from typing import List, Optional
import multiprocessing

N_COLS = 8
_DTYPE = np.float64  # 默认值，由 analyze_trajectory_folder 检测后设置

def detect_dtype(file_path):
    """检测文件的数据类型，返回 np.float64 或 np.float32 或 None"""
    file_size = os.path.getsize(file_path)
    if file_size == 0:
        return None
    row_f64 = N_COLS * 8
    row_f32 = N_COLS * 4
    if file_size % row_f64 == 0:
        raw = np.fromfile(file_path, dtype=np.float64)
        mat = raw.reshape(-1, N_COLS)
        if len(mat) >= 2 and mat[0, 0] >= 0 and mat[1, 0] > mat[0, 0]:
            return np.float64
    if file_size % row_f32 == 0:
        raw = np.fromfile(file_path, dtype=np.float32)
        mat = raw.reshape(-1, N_COLS)
        if len(mat) >= 2 and mat[0, 0] >= 0 and mat[1, 0] > mat[0, 0]:
            return np.float32
    return None

def calculate_bound_time_for_file(filepath: str) -> Optional[dict]:
    """
    分析单个轨迹文件，计算能量为负的持续时间间隔。

    假设文件格式如下（二进制 .dat，8列 float64）：
    - 第0列: time[s]
    - 第1-3列: x, y, z[km]
    - 第4-6列: vx, vy, vz[km/s]
    - 第7列: E[eV]（负值表示粒子被重力捕获）

    参数:
    filepath (str): 轨迹文件的路径。

    返回:
    Optional[dict]: 包含文件名、时间间隔和状态的字典，或 None 表示失败。
    """
    N_COLS = 8
    try:
        raw = np.fromfile(filepath, dtype=_DTYPE)
        if raw.size == 0 or raw.size % N_COLS != 0:
            return {'filename': os.path.basename(filepath), 'interval': None, 'status': 'no_valid_data'}
        mat = raw.reshape(-1, N_COLS)
        if _DTYPE == np.float32:
            mat = mat.astype(np.float64)
        time_arr   = mat[:, 0]  # time[s]
        energy_arr = mat[:, 7]  # E[eV]
    except Exception as e:
        return {
            'filename': os.path.basename(filepath),
            'interval': None,
            'status': f'error: {str(e)}'
        }

    # 找出能量为负的最早和最晚时刻
    negative_mask = energy_arr < 0
    if not np.any(negative_mask):
        return {
            'filename': os.path.basename(filepath),
            'interval': None,
            'status': 'no_negative_energy'
        }

    t_start_negative = time_arr[negative_mask][0]
    t_end_negative   = time_arr[negative_mask][-1]
    time_interval = t_end_negative - t_start_negative

    return {
        'filename': os.path.basename(filepath),
        'interval': time_interval,
        'status': 'success'
    }

def analyze_trajectory_folder(folder_path: str, max_samples: int = 1000000, num_processes: int = 8) -> float:
    """
    遍历指定文件夹中的所有trajectory文件，计算并返回平均的能量为负时间间隔。

    参数:
    folder_path (str): 包含轨迹文件的文件夹路径。
    max_samples (int): 最大处理的样本数量，默认为100000。
    num_processes (int): 并行处理的进程数量，默认为8。

    返回:
    float: 所有文件的时间间隔的平均值。如果无有效数据则返回0.0。
    """
    if not os.path.isdir(folder_path):
        print(f"错误: 文件夹 '{folder_path}' 不存在。请检查路径。", flush=True)
        return 0.0

    # 收集所有 trajectory*.dat 文件路径
    file_list = []
    for filename in os.listdir(folder_path):
        if "trajectory" in filename and filename.endswith(".dat"):
            filepath = os.path.join(folder_path, filename)
            if os.path.isfile(filepath):
                file_list.append(filepath)
                if len(file_list) >= max_samples:
                    break
    
    if not file_list:
        print("未找到任何trajectory文件。", flush=True)
        return 0.0

    # 用第一个文件检测 dtype
    global _DTYPE
    detected = detect_dtype(file_list[0])
    if detected is None:
        print(f"错误: 无法从第一个文件检测数据类型: {file_list[0]}", flush=True)
        return 0.0
    _DTYPE = detected
    print(f"检测到数据类型: {'float64' if _DTYPE == np.float64 else 'float32'} (基于 {os.path.basename(file_list[0])})", flush=True)
    
    print(f"\n==================================================", flush=True)
    print(f"开始并行扫描文件夹: '{folder_path}'", flush=True)
    print(f"文件数量: {len(file_list)}", flush=True)
    print(f"使用进程数量: {num_processes}", flush=True)
    print(f"==================================================\n", flush=True)
    
    # 使用流式处理，每完成一个立即输出，避免占用过多内存
    all_intervals = []
    success_count = 0
    error_count = 0
    processed_count = 0
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        # 使用 imap_unordered 进行流式处理，哪个完成先处理哪个
        for result in pool.imap_unordered(calculate_bound_time_for_file, file_list, chunksize=1):
            processed_count += 1
            
            if result is None:
                continue
            
            filename = result['filename']
            interval = result['interval']
            status = result['status']
            
            if status == 'success':
                print(f"--- 处理文件: {filename} ---", flush=True)
                print(f"成功: 时间间隔为: {interval:.4f}", flush=True)
                all_intervals.append(interval)
                success_count += 1
            elif status != 'no_negative_energy':
                error_count += 1
            
            # 每处理50个文件输出一次进度
            # if processed_count % 50 == 0:
            #     print(f"\n[进度] 已处理 {processed_count}/{len(file_list)} 个文件，有效间隔: {len(all_intervals)} 个\n", flush=True)
    
    print(f"\n==================================================", flush=True)
    print(f"所有文件处理完毕。共处理 {len(file_list)} 个文件。", flush=True)
    print(f"==================================================\n", flush=True)

    # 计算平均时间间隔
    if not all_intervals:
        print("\n最终结果: 未在任何文件中计算出有效的时间间隔，无法计算平均值。", flush=True)
        return 0.0
    
    # average_time = sum(all_intervals) / len(all_intervals)
    
    print(f"\n共计算出 {len(all_intervals)} 个有效的时间间隔。", flush=True)
    # print(f"它们的平均值为: {average_time:.4f}", flush=True)
    
    return 0.0

# --- 主程序入口和自测试示例 ---
if __name__ == "__main__":
    import sys
    from datetime import datetime
    
    target_folder_path = "/project/kennyng/backup_DM/Trajectoies/results_-0.301030_-32.000000"
    
    # 从路径中提取参数信息
    folder_name = os.path.basename(target_folder_path)
    params = folder_name.replace("results_", "")
    
    # 创建带参数的日志文件
    log_file = f"/project/kennyng/backup_DM/evaporation/evaporation_{params}.log"
    with open(log_file, 'a') as f:
        f.write(f"\n{'='*60}\n")
        f.write(f"程序启动时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"处理文件夹: {target_folder_path}\n")
        f.write(f"{'='*60}\n")
    
    # 重定向 stdout 和 stderr 到日志文件和控制台
    class Logger:
        def __init__(self, filename):
            self.terminal = sys.stdout
            self.log = open(filename, "a")
        
        def write(self, message):
            self.terminal.write(message)
            self.terminal.flush()
            self.log.write(message)
            self.log.flush()
        
        def flush(self):
            self.terminal.flush()
            self.log.flush()
    
    sys.stdout = Logger(log_file)
    
    final_average_time = analyze_trajectory_folder(target_folder_path, num_processes=8)

    # --- 3. 打印最终的平均时间 ---
    # print(f"\n程序执行完毕。最终计算出的平均时间间隔为: {final_average_time:.4f}", flush=True)
    print(f"日志文件已保存到: {log_file}", flush=True)
    