"""Synthetic checks for LNS record validation and comparison."""
from pathlib import Path
import importlib.util
import struct
import tempfile
import unittest

import numpy as np

spec = importlib.util.spec_from_file_location('check_lns_results',Path(__file__).with_name('check_lns_results.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
compare, load = module.compare, module.load


def fixture(folder,px=1,value=0.0):
    folder.mkdir()
    (folder/'RestartMetadata.dat').write_text(f'11 6e-6 8 7 8 {px} 1')
    (folder/'IterationNumber.dat').write_bytes(struct.pack('<iii',4,11,4))
    for layer in range(4):
        for rank in range(px):
            lx = 8//px
            cv = np.zeros((5,lx+6,7,14))
            cv[:,3:3+lx,:,3:11] = value
            shock = np.zeros((lx+6,14))
            for suffix,rr in [('flowsfg',[cv]),('shksfg',[shock,shock])]:
                payload = bytearray()
                for array in rr:
                    raw = array.astype('<f8').tobytes(order='F')
                    marker = struct.pack('<i',len(raw))
                    payload.extend(marker+raw+marker)
                (folder/f'Pert_Results_t{layer}_{rank:05d}.{suffix}').write_bytes(payload)


class Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_zero(self):
        fixture(self.root/'case')
        report, arrays = load(self.root/'case')
        self.assertEqual(len(report['records']),12)
        self.assertTrue(all(not value.any() for value in arrays.values()))

    def test_ghost_nan(self):
        fixture(self.root/'case')
        path = self.root/'case/Pert_Results_t0_00000.flowsfg'
        raw = bytearray(path.read_bytes())
        struct.pack_into('<d',raw,4,float('nan'))
        path.write_bytes(raw)
        with self.assertRaisesRegex(ValueError,'non-finite'):
            load(self.root/'case')

    def test_corrupt_record(self):
        fixture(self.root/'case')
        path = self.root/'case/Pert_Results_t0_00000.flowsfg'
        raw = bytearray(path.read_bytes())
        struct.pack_into('<i',raw,len(raw)-4,8)
        path.write_bytes(raw)
        with self.assertRaisesRegex(ValueError,'trailer'):
            load(self.root/'case')

    def test_scaling(self):
        fixture(self.root/'full',value=2.0)
        fixture(self.root/'half',value=1.0)
        a, b = (load(self.root/name)[1] for name in ('full','half'))
        self.assertTrue(all(item['passed'] for item in compare(a,b,scale=2.0,atol=0,rtol=0)))
        self.assertFalse(all(item['passed'] for item in compare(a,b,atol=0,rtol=0)))

    def test_partition(self):
        fixture(self.root/'serial',value=1.0)
        fixture(self.root/'parallel',px=2,value=1.0)
        a, b = (load(self.root/name)[1] for name in ('serial','parallel'))
        self.assertTrue(all(item['passed'] for item in compare(a,b,atol=0,rtol=0)))


if __name__ == '__main__':
    unittest.main()
