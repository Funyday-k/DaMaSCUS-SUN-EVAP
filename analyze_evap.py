import math
import argparse
import statistics

DEFAULT_PATH = '/project/kennyng/backup_DM/Trajectoies/Data/results_0.000000_-38.000000/evaporation_summary.txt'


def parse_row(parts):
    if len(parts) >= 4:
        return float(parts[2]), bool(int(float(parts[3])))
    if len(parts) >= 2:
        return float(parts[1]), False
    raise ValueError('expected at least 2 columns')


parser = argparse.ArgumentParser(description='Summarize positive evaporation times from an evaporation_summary file.')
parser.add_argument('path', nargs='?', default=DEFAULT_PATH)
parser.add_argument('--include-truncated', action='store_true')
args = parser.parse_args()

vals = []
total_rows = 0
skipped_non_evaporation = 0
skipped_truncated = 0

with open(args.path) as f:
    for line_number, line in enumerate(f, start=1):
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        if not parts:
            continue
        try:
            t_evap, truncated = parse_row(parts)
        except ValueError as exc:
            raise ValueError('Could not parse line %d: %s' % (line_number, line.rstrip())) from exc
        total_rows += 1
        if truncated and not args.include_truncated:
            skipped_truncated += 1
            continue
        if not math.isfinite(t_evap) or t_evap <= 0.0:
            skipped_non_evaporation += 1
            continue
        vals.append(t_evap)

print('Input:', args.path)
print('Rows:', total_rows)
print('Positive finite evaporation times:', len(vals))
print('Skipped non-evaporation/invalid:', skipped_non_evaporation)
print('Skipped truncated:', skipped_truncated)

logvals = [math.log10(v) for v in vals]
bins = {}
for lv in logvals:
    b = int(lv) if lv >= 0 else int(lv) - 1
    bins[b] = bins.get(b, 0) + 1
for k in sorted(bins):
    print('  log10 [%d,%d): %d' % (k, k+1, bins[k]))

if logvals:
    logvals_sorted = sorted(logvals)
    print('Median t_evap [s]:', statistics.median(vals))
    print('Median log10:', statistics.median(logvals_sorted))
    print('Min log10:', logvals_sorted[0])
    print('Max log10:', logvals_sorted[-1])
    mean_log = sum(logvals) / len(logvals)
    print('Mean log10:', mean_log)
else:
    print('No positive finite evaporation times to summarize.')
