"""Synthetic input-format regressions; not solver-accuracy validation."""
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

CHECK = Path(__file__).with_name('check_examples.py')

class InputChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.repo = Path(self.temp.name)
        self.case = self.repo / 'example' / 'circular_cylinder'
        self.restart = self.case / 'restart'
        self.restart.mkdir(parents=True)
        (self.case / 'Config.cfg').write_text(
            '&InitialSettings\nModelType=1, AnalysisType=1, Nx=8, Ny=3, Nz=8, '
            'npx0=1, npz0=1, Gamma=1.4, Mach_ref=2, '
            'IF_Continue_Calculate=1, IF_Continue_LNS=0\n/\n', encoding='ascii')
        self.cv = np.zeros((5,14,3,14), order='F')
        self.cv[0] = 1
        self.cv[4] = 3
        self.write_flow()
        self.write_records(self.restart / 'SF_Results00000.shksfg',
                           [np.ones((14,14)), np.zeros((14,14))])

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def write_records(path, arrays):
        with path.open('wb') as stream:
            for a in arrays:
                data = a.astype('<f8').tobytes(order='F')
                marker = struct.pack('<i', len(data))
                stream.write(marker + data + marker)

    def write_flow(self, coordinates=True):
        arrays = [self.cv]
        if coordinates:
            arrays += [np.zeros((14,3,14)) for _ in range(3)]
        self.write_records(self.restart / 'SF_Results00000.flowsfg', arrays)

    def check(self):
        return subprocess.run([sys.executable, str(CHECK), '--repo', str(self.repo)],
                              capture_output=True, text=True)

    def test_complete_pair_passes(self):
        self.assertEqual(self.check().returncode, 0)

    def test_missing_coordinates_fail(self):
        self.write_flow(coordinates=False)
        self.assertNotEqual(self.check().returncode, 0)

    def test_missing_shock_fails(self):
        (self.restart / 'SF_Results00000.shksfg').unlink()
        self.assertNotEqual(self.check().returncode, 0)

    def test_nonfinite_interior_fails(self):
        self.cv[0,3,1,3] = np.nan
        self.write_flow()
        self.assertNotEqual(self.check().returncode, 0)

    def test_nonpositive_pressure_fails(self):
        self.cv[4,3,1,3] = -1
        self.write_flow()
        self.assertNotEqual(self.check().returncode, 0)

    def test_corrupt_record_marker_fails(self):
        path = self.restart / 'SF_Results00000.flowsfg'
        with path.open('r+b') as stream:
            stream.seek(4 + self.cv.nbytes)
            stream.write(struct.pack('<i', 1))
        self.assertNotEqual(self.check().returncode, 0)

if __name__ == '__main__':
    unittest.main()
