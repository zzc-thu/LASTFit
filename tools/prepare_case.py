"""Stage an isolated short run without modifying reference case data."""
import argparse
import re
import shutil
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('case', choices=('circular_cylinder', 'parabolic_nonlinear',
                                    'parabolic_lns', 'blunt_cone_aoa1',
                                    'hifire5_acoustic', 'hifire5_steady'))
parser.add_argument('--run-dir', type=Path, required=True)
parser.add_argument('--solver', type=Path, required=True)
parser.add_argument('--steps', type=int, default=10, help='Last step for a new short run (default: 10)')
parser.add_argument('--allow-unverified', action='store_true', help='Allow the currently failing LNS diagnostic')
args = parser.parse_args()
repo = Path(__file__).resolve().parents[1]
if args.case == 'parabolic_lns' and not args.allow_unverified:
    parser.error('LNS currently fails the finite-value regression; use --allow-unverified for diagnosis only')
folder = 'hifire5_acoustic' if args.case == 'hifire5_steady' else args.case
case_dir = repo / 'example' / folder
config_name = 'Config_steady.cfg' if args.case == 'hifire5_steady' else 'Config.cfg'
run = args.run_dir.resolve()
if args.steps < 4:
    parser.error('--steps must be at least 4 for LNS startup')
if run.exists():
    parser.error('--run-dir must not exist; existing runs are never overwritten')
if run.is_relative_to(repo / 'example') or run.is_relative_to(repo / 'src'):
    parser.error('run directory must not be inside reference examples or source')
if not args.solver.is_file():
    parser.error('--solver must point to a compiled executable')
config = (case_dir / config_name).read_text(encoding='ascii')
for key in ('TotalTimeSteps', 'LNS_STEPS_MAX', 'StepsWriteData'):
    config = re.sub(r'(?im)^(\s*' + key + r'\s*=\s*)\d+',
                    lambda m: m[1] + str(args.steps), config)
config = re.sub(r'(?im)^(\s*StepsOutputScreen\s*=\s*)\d+', lambda m: m[1] + '1', config)

copies = []
def checkpoint(source, target):
    pairs = sorted(source.glob('SF_Results*.flowsfg'))
    if not pairs:
        parser.error(f'No flow checkpoint in {source}')
    for flow in pairs:
        shock = flow.with_suffix('.shksfg')
        if not shock.is_file():
            parser.error(f'Missing paired shock checkpoint: {shock}')
        copies.extend((p, target / p.name) for p in (flow, shock))
        time = flow.with_suffix('.time')
        if time.is_file():
            copies.append((time, target / time.name))

if args.case in ('circular_cylinder', 'parabolic_nonlinear'):
    checkpoint(case_dir / 'restart', run / 'RESU')
if args.case in ('parabolic_lns', 'parabolic_nonlinear'):
    checkpoint(repo / 'example' / 'parabolic_lns' / 'restart', run / 'RESU_STEADY')
if args.case == 'hifire5_acoustic':
    checkpoint(case_dir / 'restart' / 'unsteady', run / 'RESU_STEADY')

for name in ('INIT', 'CheckFiles', 'JACO', 'RESU', 'RESU_STEADY', 'LNSResults', 'Pert'):
    (run / name).mkdir(parents=True, exist_ok=True)
(run / 'Config.cfg').write_text(config, encoding='ascii')
shutil.copy2(args.solver, run / args.solver.name)
for source, target in copies:
    shutil.copy2(source, target)
px = int(re.search(r'(?im)^\s*npx0\s*=\s*(\d+)', config)[1])
pz = int(re.search(r'(?im)^\s*npz0\s*=\s*(\d+)', config)[1])
print(f'Staged {args.case}: {run}')
print(f'MPI ranks: {px*pz}; final step: {args.steps}')
print(f'Run from that directory: mpirun -np {px*pz} ./{args.solver.name}')
print('A short run checks startup only, not convergence or manuscript reproducibility.')
