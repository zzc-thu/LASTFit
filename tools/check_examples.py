"""Audit the scalar example settings and sequential Fortran restart records."""
import argparse
import json
import re
import struct
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description='Check example restart shapes and interior physical values.')
parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
parser.add_argument('--report', type=Path, help='Optional JSON report')
args = parser.parse_args()
REPO = args.repo.resolve()


def settings(path):
    # These example namelists contain numeric scalar assignments only.
    clean = '\n'.join(line.split('!')[0] for line in path.read_text().splitlines())
    pairs = re.findall(r'(\w+)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?)', clean)
    return {key.lower(): float(value.lower().replace('d', 'e')) for key, value in pairs}


def records(path):
    result = []
    with path.open('rb') as stream:
        while header := stream.read(4):
            if len(header) != 4:
                raise ValueError('Incomplete record marker')
            size, = struct.unpack('<i', header)
            if size < 0 or size > path.stat().st_size:
                raise ValueError('Unsupported record marker or endianness')
            data = stream.read(size)
            trailer = stream.read(4)
            if len(data) != size or trailer != header:
                raise ValueError('Record trailer/payload mismatch')
            result.append(data)
    return result


def stats(a):
    finite = np.isfinite(a)
    return {'count': int(a.size), 'nonfinite': int((~finite).sum()),
            'min': float(a[finite].min()) if finite.any() else None,
            'max': float(a[finite].max()) if finite.any() else None}


results = []
for config in sorted((REPO / 'example').rglob('*.cfg')):
    c = settings(config)
    nx, ny, nz, px, pz = [int(c[k]) for k in ('nx', 'ny', 'nz', 'npx0', 'npz0')]
    entry = {'config': config.relative_to(REPO).as_posix(), 'settings': c, 'files': []}
    folder = config.parent / 'restart'
    if config.name == 'Config_steady.cfg' or config.name == 'Config_E4_fine.cfg':
        entry['files_not_applicable'] = True
    else:
        for rank in range(px * pz):
            lx = nx // px + (rank % px < nx % px)
            lz = nz // pz + (rank // px < nz % pz)
            n = (lx + 6) * ny * (lz + 6)
            s = (lx + 6) * (lz + 6)
            for suffix, expected in [('flowsfg', [5*n*8, n*8, n*8, n*8]), ('shksfg', [s*8, s*8])]:
                hits = list(folder.rglob(f'SF_Results{rank:05d}.{suffix}')) if folder.exists() else []
                item = {'rank': rank, 'kind': suffix, 'expected_record_bytes': expected}
                if len(hits) != 1:
                    item['error'] = f'Expected one file, found {len(hits)}'
                else:
                    path = hits[0]
                    item['path'] = path.relative_to(REPO).as_posix()
                    try:
                        rr = records(path)
                        item['actual_record_bytes'] = [len(r) for r in rr]
                        item['format_matches'] = item['actual_record_bytes'] == expected
                        item['records_stats'] = [stats(np.frombuffer(r, dtype='<f8')) for r in rr]
                        if suffix == 'flowsfg' and len(rr[0]) == expected[0]:
                            cv = np.frombuffer(rr[0], dtype='<f8').reshape((5, lx+6, ny, lz+6), order='F')[:, 3:3+lx, :, 3:3+lz]
                            rho = cv[0]
                            with np.errstate(divide='ignore', invalid='ignore'):
                                p = (c['gamma'] - 1) * (cv[4] - (cv[1:4]**2).sum(axis=0) / (2*rho))
                            item['interior_density'] = stats(rho)
                            item['interior_pressure'] = stats(p)
                            item['nonpositive_density'] = int((rho <= 0).sum())
                            item['nonpositive_pressure'] = int((p <= 0).sum())
                            item['interior_cv'] = stats(cv)
                            item['wall_temperature'] = stats(c['gamma'] * c['mach_ref']**2 * p[:, 0, :] / rho[:, 0, :])
                            bad = np.argwhere(~np.isfinite(np.frombuffer(rr[0], dtype='<f8').reshape((5, lx+6, ny, lz+6), order='F')))
                            if bad.size:
                                item['nonfinite_index_bounds_zero_based'] = [bad.min(axis=0).tolist(), bad.max(axis=0).tolist()]
                    except (ValueError, IndexError) as exc:
                        item['error'] = str(exc)
                entry['files'].append(item)
    results.append(entry)

errors = 0
for entry in results:
    c = entry['settings']
    problems = []
    warnings = []
    if 'hifire5_acoustic' in entry['config'] and c['modeltype'] != 4:
        problems.append('HIFiRE requires ModelType=4 in the current source')
    for item in entry['files']:
        if item.get('error') or not item.get('format_matches', False):
            problems.append(f"rank {item['rank']} {item['kind']}: record mismatch or missing file")
        if item.get('interior_cv', {}).get('nonfinite', 0):
            problems.append(f"rank {item['rank']}: non-finite interior flow")
        if item.get('nonpositive_density', 0) or item.get('nonpositive_pressure', 0):
            problems.append(f"rank {item['rank']}: non-positive density or pressure")
        if any(a['nonfinite'] for a in item.get('records_stats', [])):
            warnings.append(f"rank {item['rank']} {item['kind']}: non-finite stored values; inspect halos")
    if c['analysistype'] == 2 and c['if_continue_calculate'] == 1:
        warnings.append('Verify complete per-rank .time files for continuous restart; legacy data are warm starts')
    if c['analysistype'] == 3:
        warnings.append('LNS uses LNS_dt, not FixedDeltaTime; input validity is not an operator regression')
    entry['errors'] = problems
    entry['warnings'] = sorted(set(warnings))
    errors += len(problems)
    print(('FAIL' if problems else 'PASS'), entry['config'])
    for msg in problems:
        print('  ERROR:', msg)
    for msg in entry['warnings']:
        print('  WARNING:', msg)
if args.report:
    args.report.write_text(json.dumps(results, indent=2), encoding='utf-8')
raise SystemExit(1 if errors else 0)
