#!/usr/bin/env python3
"""Bounded SMS CLI/sampling tests; optional same-seed upstream regression.

Uses only the Python standard library and fresh temporary directories.
Usage: python3 tests/sms_integration.py build/pimc.e [--baseline UPSTREAM_EXE]
"""
import argparse
import math
from pathlib import Path
import subprocess
import tempfile

BASE = ['--canonical', '--geometry=prism', '--size=10000', '--number_particles=4',
        '--temperature=7e-9', '--mass=155.8', '--number_time_slices=16',
        '--update_length=4', '--action=primitive', '--external=free',
        '--number_eq_steps=0', '--bin_size=2', '--number_bins_stored=2',
        '--update=center of mass', '--update=bisection', '--estimator=energy',
        '--no_save_state', '--process=42']
SMS = ['--interaction=sms', '--sms_dipole_debye=4.6', '--sms_omega_mhz=4',
       '--sms_delta_mhz=22', '--sms_xi_fraction=0.3']


def replace(options, key, value):
    return [s for s in options if not s.startswith('--'+key+'=')] + ['--'+key+'='+str(value)]


def run(exe, folder, options, error=None):
    folder.mkdir()
    result = subprocess.run([str(exe), *options], cwd=folder, capture_output=True,
                            text=True, timeout=45)
    output = result.stdout + result.stderr
    if error is not None:
        if result.returncode == 0 or error not in output:
            raise AssertionError(f'Expected rejection {error!r}:\n{output}')
        if (folder/'OUTPUT').exists():
            raise AssertionError('Invalid options created simulation output')
    elif result.returncode or 'Measurement complete.' not in output:
        raise AssertionError(output)
    return folder


def one_file(folder, kind):
    files = list((folder/'OUTPUT').glob(f'*-{kind}-*.dat'))
    if len(files) != 1:
        raise AssertionError(f'Expected one {kind} in {folder}: {files}')
    return files[0]


def energy(folder):
    rows = [list(map(float, line.split())) for line in one_file(folder, 'estimator').read_text().splitlines()
            if line.strip() and not line.startswith('#')]
    if len(rows) != 2 or any(not math.isfinite(v) for row in rows for v in row):
        raise AssertionError(f'Invalid energy rows: {rows}')
    return rows


def state(folder):
    text = one_file(folder, 'state').read_text()
    # Counters precede the bead array. Compare coordinates, links, mask and RNG.
    return text[text.index('('):]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    exe = args.executable.resolve()
    with tempfile.TemporaryDirectory(prefix='sms-test-') as temporary:
        root = Path(temporary)
        good = run(exe, root/'prism', BASE+SMS)
        energy(good)
        log = one_file(good, 'log').read_text()
        for marker in ('SMS C3 [K A^3]', 'SMS C6 [K A^6]', 'no images or tail correction'):
            if marker not in log:
                raise AssertionError(f'Missing metadata: {marker}')
        cylinder = replace(BASE, 'geometry', 'cylinder')+['--radius=2000']
        energy(run(exe, root/'cylinder', cylinder+SMS))

        for i, option in enumerate(SMS[1:]):
            run(exe, root/f'missing{i}', BASE+[s for s in SMS if s != option], error='SMS requires')
        for i, (key, value) in enumerate([
                ('sms_omega_mhz', 0), ('sms_dipole_debye', -1), ('sms_delta_mhz', -1),
                ('sms_xi_fraction', 1), ('sms_xi_fraction', -0.1),
                ('sms_dipole_debye', 'nan'), ('sms_omega_mhz', 'inf'),
                ('sms_dipole_debye', '1e300')]):
            run(exe, root/f'invalid{i}', BASE+replace(SMS, key, value), error='ERROR:')
        for action in ('pair_product', 'gsf', 'li_broughton'):
            run(exe, root/action, replace(BASE, 'action', action)+SMS, error='SMS requires')
        run(exe, root/'virial', replace(BASE, 'estimator', 'virial')+SMS, error='radial virial')
        for i, cutoff in enumerate(('0', '-1', 'nan', '1e300')):
            run(exe, root/f'cutoff{i}', BASE+SMS+['--potential_cutoff='+cutoff], error='SMS cutoff')

        # Almost static two-particle paths test the potential adapter and its cutoff
        # through the actual primitive energy measurement, not a copied kernel.
        xyz = root/'pair.xyz'
        xyz.write_text('2\nangstrom\nX -1000 0 0\nX 1000 0 0\n')
        fixed = replace(BASE, 'update', 'center of mass')+['--com_delta=1e-12', '--xyz='+str(xyz)]
        c3, c6 = 408.69644992609824, 87552088116.64389
        s = math.sin(0.3*math.pi/2)
        expected = c3*(-1+3*s)/2000**3+c6*(1-s*s)/2000**6
        for label, cutoff, value in [('inside', 3000, expected), ('outside', 1000, 0)]:
            rows = energy(run(exe, root/label, fixed+SMS+['--potential_cutoff='+str(cutoff)]))
            for row in rows:
                for column in (1, 3):  # total V and operator V_int
                    if not math.isclose(row[column], value, rel_tol=2e-5, abs_tol=1e-18):
                        raise AssertionError(f'{label}: V={row[column]}, expected {value}')

        # Across the periodic x face, separation 9000 A must reduce to 1000 A.
        xyz.write_text('2\nangstrom\nX -4500 0 0\nX 4500 0 0\n')
        expected = c3*(-1+3*s)/1000**3+c6*(1-s*s)/1000**6
        for row in energy(run(exe, root/'periodic', fixed+SMS+['--potential_cutoff=1500'])):
            if not math.isclose(row[3], expected, rel_tol=2e-5):
                raise AssertionError('Minimum-image interaction mismatch')

        # Old-format text checkpoints must still load into a fresh SMS segment.
        restart = replace(BASE, 'start_with_state', one_file(good, 'state'))+SMS
        energy(run(exe, root/'restart', restart))
        if args.baseline:
            for interaction in ('free', 'harmonic'):
                options = BASE+['--interaction='+interaction]
                if interaction == 'harmonic':
                    options += ['--omega=1e-9']
                before = run(args.baseline.resolve(), root/(interaction+'-before'), options)
                after = run(exe, root/(interaction+'-after'), options)
                if state(before) != state(after) or energy(before) != energy(after):
                    raise AssertionError(f'Non-SMS regression: {interaction}')
                # A state written by unmodified upstream is accepted as an SMS seed.
                energy(run(exe, root/(interaction+'-old-state'), BASE+SMS+
                           ['--start_with_state='+str(one_file(before, 'state'))]))
            print('Same-seed upstream free/harmonic states and energy match exactly')
    print('SMS prism/cylinder, metadata, rejection, cutoff, and checkpoint tests passed')


if __name__ == '__main__':
    main()
