#!/usr/bin/env python3
# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: AGPL-3.0-or-later

"""TDDFT regression against pyscf (refs.json, made by make_refs.py).
Usage: python3 run_tests.py [path/to/Direwolf]   (run from anywhere; ~1 min)
Gas phase: 3e-3 eV. Solvent (C-PCM, non-equilibrium): the cavities differ, 6e-2 eV."""
import json, os, re, subprocess, sys, tempfile
here = os.path.dirname(os.path.abspath(__file__))
exe = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(here, '..', '..', 'Direwolf')
refs = json.load(open(os.path.join(here, 'refs.json')))
H2CO = "C 0.0 0.0 0.0\nO 0.0 0.0 1.22\nH 0.0 0.94 -0.58\nH 0.1 -0.94 -0.58"
OH = "O 0.0 0.0 0.0\nH 0.0 0.0 0.97"
def inp(atoms, mult, extra, solvent, pc=False):
    s = ("&molecule\n imult=%d\n icharge=0\n functional='PBE0'\n baselabel='def2-svp'\n verbose=1\n"
         " J='exact'\n K='exact'\n tddft_nroots=3\n calc_force=.false.\n%s&end\n" % (mult, extra))
    if solvent:
        s += "&cosmo\n cosmo_on=.true.\n cosmo_solvent='water'\n&end\n"
    s += "&atoms\n" + atoms + "\n&end\n"
    if pc:   # two external point charges (Angstrom), as in make_refs.py
        s += "&pointcharges\n npc=2\n pc_q=0.5,-0.5\n pc_x=3.0,-3.0\n pc_y=0.0,0.0\n pc_z=2.0,0.0\n&end\n"
    return s
# OH: state 1 is the Goldstone mode of the degenerate pi pair (exactly 0 for a rotation-invariant grid;
# 0.003 eV in pyscf, 0.03 eV here, 0.045 eV when the molecule is tilted), so only states 2-3 are compared (gas)
# and state 2 (water; state 3 is a diffuse state, 0.14 eV per 0.1 of cavity radius scale)
use = {'oh_gas_uks': [1, 2], 'oh_water_uks': [1], 'oh_gas_uks_pc': [1, 2]}
# name, atoms, multiplicity, extra keys, solvent, block header, tolerance (eV)
cases = [
  ('h2co_gas_singlet',  H2CO, 1, '', False, 'Singlet excited states', 3e-3),
  ('h2co_gas_triplet',  H2CO, 1, ' tddft_triplets=.true.\n', False, 'Triplet excited states', 3e-3),
  ('h2co_gas_tda',      H2CO, 1, ' tddft_tda=.true.\n', False, 'Singlet excited states', 3e-3),
  ('h2co_water_singlet', H2CO, 1, '', True, 'Singlet excited states', 6e-2),
  ('h2co_water_triplet', H2CO, 1, ' tddft_triplets=.true.\n', True, 'Triplet excited states', 6e-2),
  ('oh_gas_uks',        OH, 2, '', False, 'UKS linear response', 3e-3),
  ('oh_water_uks',      OH, 2, '', True, 'UKS linear response', 6e-2),
  ('h2co_gas_singlet_pc', H2CO, 1, '', False, 'Singlet excited states', 3e-3),
  ('oh_gas_uks_pc',     OH, 2, '', False, 'UKS linear response', 3e-3),
  ('h2co_water_singlet_pc', H2CO, 1, '', True, 'Singlet excited states', 6e-2),
]
def energies(text, header):
    i = text.index(header)
    rows = re.findall(r'^\s+(\d+)\s+(-?\d+\.\d{4})\s+-?\d+\.\d{4}\s+\S+\s+\S+', text[i:], re.M)
    return [float(r[1]) for r in rows[:3]]
bad = 0
with tempfile.TemporaryDirectory() as d:
    for name, atoms, mult, extra, solv, header, tol in cases:
        f = os.path.join(d, name + '.inp')
        open(f, 'w').write(inp(atoms, mult, extra, solv, name.endswith('_pc')))
        subprocess.run([exe, f], cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            e = energies(open(os.path.join(d, name + '.out')).read(), header)
        except Exception as ex:
            print('FAIL %-20s no result (%s)' % (name, ex)); bad += 1; continue
        idx = use.get(name, [0, 1, 2])
        err = max(abs(e[i] - refs[name][i]) for i in idx) if len(e) == 3 else 9.9
        ok = len(e) == 3 and err < tol
        bad += not ok
        print('%s %-20s max|dE| %.1e eV (tol %.0e)' % ('ok  ' if ok else 'FAIL', name, err, tol))
sys.exit(1 if bad else 0)
