#!/usr/bin/env python3
# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Polarizability of water PBE0/def2-SVP with a solvent model and/or external point charges.
Static diagonal against the second derivative of the energy by finite field in pyscf (field 2e-3 a.u., conv_tol 1e-12,
grid level 5; fd.py in this directory prints the values), dynamic tensor against the sum over all TDDFT states of the same
run settings (alpha_ii(w) = sum_n 2 w_n mu_n,i^2 / (w_n^2 - w^2); water/def2-SVP has 5 x 19 = 95 states).
Tolerances: 2e-3 a.u. (C-PCM, point charges: cavity and grid), 1.5e-2 (SMD: the two SMD cavities differ), 2e-3 (sum over states:
4-decimal output).   Usage: python3 run_tests.py [path/to/Direwolf]   (~1 min)"""
import os, re, subprocess, sys, tempfile
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
exe = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(here, '..', '..', 'Direwolf')
ATOMS = "O 0.0000 0.0000 0.1173\nH 0.0000 0.7572 -0.4692\nH 0.0000 -0.7572 -0.4692"
PC = "&pointcharges\n npc=2\n pc_q=0.5,-0.5\n pc_x=3.0,-3.0\n pc_y=0.0,0.0\n pc_z=2.0,0.0\n&end\n"
SOLV = {'none': '', 'cpcm': "&cosmo\n solvent='water'\n convention='cpcm'\n&end\n", 'smd': "&cosmo\n smd=.true.\n solvent='water'\n&end\n"}
# name: (solvent, point charges, reference alpha_xx,yy,zz from pyscf finite field, tolerance)
REF = {
  'gas':      ('none', False, (3.1256, 7.0262, 5.2587), 2e-3),
  'cpcm':     ('cpcm', False, (3.3358, 7.6195, 5.5985), 2e-3),
  'pc':       ('none', True,  (3.1087, 6.9738, 5.2159), 2e-3),
  'cpcm_pc':  ('cpcm', True,  (3.3164, 7.5535, 5.5486), 2e-3),
  'smd':      ('smd',  False, (3.4480, 7.7353, 5.7447), 1.5e-2),
}
def inp(solv, pc, extra):
    s = ("&molecule\n imult=1\n icharge=0\n functional='PBE0'\n baselabel='def2svp'\n verbose=1\n J='exact'\n K='exact'\n"
         " calc_force=.false.\n%s&end\n&atoms\n%s\n&end\n" % (extra, ATOMS))
    return s + SOLV[solv] + (PC if pc else '')
def run(d, name, text):
    f = os.path.join(d, name + '.inp'); open(f, 'w').write(text)
    subprocess.run([exe, f], cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return open(os.path.join(d, name + '.out')).read()
def tensors(t):
    out = []
    for b in re.split(r'Static dipole polarizability|Dynamic dipole polarizability', t)[1:]:
        rows = re.findall(r'^\s+([xyz])\s+(-?\d+\.\d{4})\s+(-?\d+\.\d{4})\s+(-?\d+\.\d{4})', b, re.M)[:3]
        out.append(np.array([[float(x) for x in r[1:]] for r in rows]))
    return out
bad = 0
with tempfile.TemporaryDirectory() as d:
    for name, (solv, pc, ref, tol) in REF.items():
        try:
            diag = np.diag(tensors(run(d, name, inp(solv, pc, " runtype='polar'\n")))[0])
            err = max(abs(diag - np.array(ref)))
        except Exception as ex:
            print('FAIL %-8s no result (%s)' % (name, ex)); bad += 1; continue
        ok = err < tol; bad += not ok
        print('%s %-8s alpha_xx,yy,zz %s  max|d| %.1e a.u. (tol %.0e)' % ('ok  ' if ok else 'FAIL', name, np.round(diag, 4), err, tol))
    # dynamic tensor (optical solvent response) against the sum over all states, with the solvent and the point charges
    for name, solv, pc in (('dyn_cpcm', 'cpcm', False), ('dyn_cpcm_pc', 'cpcm', True)):
        try:
            om = 0.05
            dyn = tensors(run(d, name, inp(solv, pc, " runtype='polar'\n polar_omega='%g'\n" % om)))[1]
            t2 = run(d, name + '_td', inp(solv, pc, " tddft_nroots=95\n"))
            rr = re.findall(r'^\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+\S+\s+\S+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', t2, re.M)
            W = np.array([float(r[2]) for r in rr]); mu = np.array([[float(r[3]), float(r[4]), float(r[5])] for r in rr])
            sos = np.array([[np.sum(2*W*mu[:, i]*mu[:, j]/(W**2 - om**2)) for j in range(3)] for i in range(3)])
            err = abs(dyn - sos).max()
        except Exception as ex:
            print('FAIL %-8s no result (%s)' % (name, ex)); bad += 1; continue
        ok = err < 2e-3 and len(W) >= 90; bad += not ok
        print('%s %-11s omega %.2f: direct vs sum over %d states, max|d| %.1e a.u. (tol 2e-3)' % ('ok  ' if ok else 'FAIL', name, om, len(W), err))
sys.exit(1 if bad else 0)
