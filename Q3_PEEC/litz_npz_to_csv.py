#!/usr/bin/env python3
"""Temporary numpy converter: litz_q3 .npz -> CSV for MATLAB.

Only needs numpy (no scipy/h5py). Exports, per case:
  litz_<case>_centerlines.csv : station,strand,x,y,z   (SI, m; z resampled to ~500 stations)
  litz_<case>_currents.csv    : strand,Re,Im           (COMSOL per-strand complex current)
Run: python3 litz_npz_to_csv.py
"""
import os
import numpy as np

SRC = os.path.join(os.path.dirname(__file__), '..', 'litz_q3', 'data')
OUT = os.path.dirname(__file__)

cases = [
    ('rigid_twist',    'rigid_twist_centerlines'),
    ('shell_exchange', 'shell_exchange_centerlines'),
    ('full_exchange',  'full_exchange_centerlines'),
    ('final_484',      'final_484_centerlines'),
]

meta_rows = []
for tag, f in cases:
    d = np.load(os.path.join(SRC, f + '.npz'), allow_pickle=True)
    xyz = d['xyz_m']                      # (nz, N, 3)
    nz, N, _ = xyz.shape
    step = max(1, int(round(nz / 500.0)))
    xyz = xyz[::step]
    nzs = xyz.shape[0]
    st = np.repeat(np.arange(nzs), N)
    sd = np.tile(d['strand_ids'], nzs)
    flat = xyz.reshape(-1, 3)
    out = np.column_stack([st, sd, flat])
    np.savetxt(os.path.join(OUT, f'litz_{tag}_centerlines.csv'), out,
               delimiter=',', header='station,strand,x,y,z', comments='', fmt='%.10g')
    meta_rows.append((tag, N, int(nzs), float(d['period_m']), str(d['topology'])))
    print(f'{tag:16s} N={N:4d} stations={nzs:4d} (step {step}) period={float(d["period_m"]):.3f} m topology={d["topology"]}')

# COMSOL per-strand currents for the case that matches final_484 (P=0.5, N=484)
for tag, path in [('final_484', 'raw/full_exchange_484_verified/result.npz'),
                  ('initial_400', 'raw/initial_rigid_400/result.npz'),
                  ('fullexch_400', 'raw/full_exchange_400_verified/result.npz')]:
    d = np.load(os.path.join(SRC, path), allow_pickle=True)
    cur = d['currents'].ravel()
    np.savetxt(os.path.join(OUT, f'litz_{tag}_currents.csv'),
               np.column_stack([np.arange(len(cur)), cur.real, cur.imag]),
               delimiter=',', header='strand,Re,Im', comments='', fmt='%.12g')
    print(f'litz_{tag}_currents.csv  n={len(cur)}  sum={cur.sum():.6g}')

with open(os.path.join(OUT, 'litz_cases_meta.csv'), 'w') as fh:
    fh.write('case,N,stations,period_m,topology\n')
    for r in meta_rows:
        fh.write(f'{r[0]},{r[1]},{r[2]},{r[3]:.6g},{r[4]}\n')
print('done ->', OUT)
