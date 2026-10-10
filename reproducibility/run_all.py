"""Run the eight exact checks accompanying Clipped Gram energies and zero counts.

Use Python with assertions enabled and Node on PATH or pass --node PATH.
The runner preserves checker stdout and reports each process exit status.
"""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--node', help='Path to Node; defaults to node on PATH.')
    args = parser.parse_args()
    if not __debug__ or sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE', '') not in ('', '0'):
        raise RuntimeError('Run with assertions enabled; remove -O, -OO and PYTHONOPTIMIZE.')
    node = args.node or shutil.which('node')
    if not node:
        raise RuntimeError('Node is required for the two-lattice check; pass --node PATH.')
    here = Path(__file__).resolve().parent
    jobs = [
        (sys.executable, 'verify_clipped_gram_constants.py', []),
        (sys.executable, 'positive_elevation_check.py', []),
        (sys.executable, 'positive_elevation_pressure_check.py', []),
        (sys.executable, 'verify_positive_elevation_independent.py', []),
        (node, 'structural_union_exact_check.js', []),
        (sys.executable, 'verify_structural_tex_scalars.py', []),
        (sys.executable, 'verify_twelve_fibre.py', ['--certify']),
        (sys.executable, 'verify_twelve_fibre_independent.py', []),
    ]
    missing = [name for _, name, _ in jobs if not (here / name).is_file()]
    if missing:
        raise RuntimeError('Missing package checks: ' + ', '.join(missing))
    outcomes = []
    for executable, name, switches in jobs:
        print('\nRunning ' + name + (' ' + ' '.join(switches) if switches else ''), flush=True)
        started = time.perf_counter()
        result = subprocess.run([executable, str(here / name), *switches], cwd=here, check=False)
        elapsed = time.perf_counter() - started
        outcomes.append((name, result.returncode, elapsed))
        print(f'Exit {result.returncode}; elapsed {elapsed:.3f} seconds.', flush=True)
    print('\nExact-check process summary:', flush=True)
    for name, code, elapsed in outcomes:
        print(f'{name}: exit={code}, seconds={elapsed:.3f}', flush=True)
    passed = sum(code == 0 for _, code, _ in outcomes)
    print(f'{passed}/{len(jobs)} processes exited successfully.', flush=True)
    print('These checks supplement the manuscript proofs and the separate Lean development.', flush=True)
    return 0 if passed == len(jobs) else 1


if __name__ == '__main__':
    sys.exit(main())
