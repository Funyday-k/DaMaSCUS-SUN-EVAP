import math

vals = []
with open('/project/kennyng/backup_DM/Trajectoies/Data/results_0.000000_-38.000000/evaporation_summary.txt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        if len(parts) >= 2:
            vals.append(float(parts[1]))

zeros = sum(1 for v in vals if v == 0.0)
print('Total:', len(vals), 'Zeros:', zeros)

logvals = [math.log10(v) for v in vals if v > 0]
bins = {}
for lv in logvals:
    b = int(lv) if lv >= 0 else int(lv) - 1
    bins[b] = bins.get(b, 0) + 1
for k in sorted(bins):
    print('  log10 [%d,%d): %d' % (k, k+1, bins[k]))
print('Non-zero count:', len(logvals))

if logvals:
    logvals_sorted = sorted(logvals)
    n = len(logvals_sorted)
    print('Median log10:', logvals_sorted[n//2])
    print('Min log10:', logvals_sorted[0])
    print('Max log10:', logvals_sorted[-1])
    mean_log = sum(logvals) / len(logvals)
    print('Mean log10:', mean_log)
