"""Audit LNS history records and compare linearity or restart results."""
import argparse
import json
from pathlib import Path
import struct

import numpy as np


def records(path):
    result = []
    with path.open('rb') as stream:
        while header := stream.read(4):
            if len(header) != 4:
                raise ValueError(f'{path.name}: incomplete record marker')
            length, = struct.unpack('<i', header)
            if length <= 0 or length % 8 or length > path.stat().st_size:
                raise ValueError(f'{path.name}: invalid record size')
            data = stream.read(length)
            if len(data) != length or stream.read(4) != header:
                raise ValueError(f'{path.name}: record trailer/payload mismatch')
            result.append(np.frombuffer(data, dtype='<f8'))
    return result


def load(folder):
    fields = (folder / 'RestartMetadata.dat').read_text().split()
    if len(fields) < 7:
        raise ValueError('Incomplete restart metadata')
    step, dt = int(fields[0]), float(fields[1].lower().replace('d','e'))
    nx, ny, nz, px, pz = map(int, fields[2:7])
    if step < 4 or not np.isfinite(dt) or dt <= 0 or min(nx,ny,nz,px,pz) < 1 or px > nx or pz > nz:
        raise ValueError('Invalid restart metadata')
    iteration = (folder / 'IterationNumber.dat').read_bytes()
    if len(iteration) != 12 or struct.unpack('<iii',iteration) != (4,step,4):
        raise ValueError('Iteration number does not match restart metadata')
    arrays = {}
    stats = []
    for layer in range(4):
        cv = np.empty((5,nx,ny,nz))
        height = np.empty((nx,nz))
        velocity = np.empty((nx,nz))
        z_offset = 0
        for rz in range(pz):
            lz = nz//pz + (rz < nz%pz)
            x_offset = 0
            for rx in range(px):
                lx = nx//px + (rx < nx%px)
                rank = rz*px+rx
                x_slice = slice(x_offset,x_offset+lx)
                z_slice = slice(z_offset,z_offset+lz)
                for suffix, expected in [('flowsfg',1),('shksfg',2)]:
                    path = folder / f'Pert_Results_t{layer}_{rank:05d}.{suffix}'
                    rr = records(path)
                    if len(rr) != expected:
                        raise ValueError(f'{path.name}: wrong record count')
                    for index, values in enumerate(rr):
                        shape = (5,lx+6,ny,lz+6) if suffix == 'flowsfg' else (lx+6,lz+6)
                        if values.size != np.prod(shape):
                            raise ValueError(f'{path.name}: shape does not match metadata')
                        field = values.reshape(shape,order='F')
                        nonfinite = int((~np.isfinite(field)).sum())
                        stats.append(dict(file=path.name,record=index,nonfinite=nonfinite))
                        if nonfinite:
                            raise ValueError(f'{path.name}: {nonfinite} non-finite stored values')
                        if suffix == 'flowsfg':
                            cv[:,x_slice,:,z_slice] = field[:,3:3+lx,:,3:3+lz]
                        else:
                            (height if index == 0 else velocity)[x_slice,z_slice] = field[3:3+lx,3:3+lz]
                x_offset += lx
            z_offset += lz
        arrays[f'cv_t{layer}'] = cv
        arrays[f'height_t{layer}'] = height
        arrays[f'velocity_t{layer}'] = velocity
    return dict(step=step,dt=dt,grid=[nx,ny,nz],ranks=[px,pz],records=stats), arrays


def compare(candidate, reference, scale=1.0, atol=1e-12, rtol=1e-12):
    result = []
    for key, value in candidate.items():
        expected = scale*reference[key]
        difference = np.abs(value-expected)
        result.append(dict(field=key,max_difference=float(difference.max()),
                           passed=bool(np.all(difference <= atol+rtol*np.abs(expected)))))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results',type=Path,help='Directory containing the four history levels')
    parser.add_argument('--reference',type=Path,help='Another run, possibly with a different MPI partition')
    parser.add_argument('--scale',type=float,default=1.0,help='Multiplier applied to the reference histories')
    parser.add_argument('--atol',type=float,default=1e-12)
    parser.add_argument('--rtol',type=float,default=1e-12)
    parser.add_argument('--zero',action='store_true',help='Require exactly zero interior perturbations')
    parser.add_argument('--report',type=Path)
    args = parser.parse_args()
    if not np.isfinite([args.scale,args.atol,args.rtol]).all() or min(args.atol,args.rtol) < 0:
        parser.error('Scale and tolerances must be finite; tolerances must be non-negative')
    try:
        report, arrays = load(args.results)
        report['passed'] = True
        if args.zero:
            report['zero'] = all(not np.any(value != 0) for value in arrays.values())
            report['passed'] &= report['zero']
        if args.reference:
            ref_report, ref_arrays = load(args.reference)
            if any(report[key] != ref_report[key] for key in ('step','dt','grid')):
                raise ValueError('Reference time/grid differs')
            report['comparisons'] = compare(arrays,ref_arrays,args.scale,args.atol,args.rtol)
            report['passed'] &= all(item['passed'] for item in report['comparisons'])
    except (OSError,ValueError,struct.error) as error:
        report = dict(passed=False,error=str(error))
    if args.report:
        args.report.write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report,indent=2))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
