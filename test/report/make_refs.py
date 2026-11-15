# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: AGPL-3.0-or-later

# Regenerates the pyscf part of refs.json (not needed to run the tests): python3 make_refs.py
# UKS PBE0/def2-SVP of square cyclobutadiene, followed to the symmetry-broken solution (stability = 2), Hessian
# frequencies. A closed-shell Hessian on this solution gives -3100 cm-1 for the lowest mode.
import json, os
import numpy as np
from pyscf import gto, dft
from pyscf.hessian import thermo
CBD = ("C 0.720000 0.720000 0.0;C -0.720000 0.720000 0.0;C -0.720000 -0.720000 0.0;C 0.720000 -0.720000 0.0;"
       "H 1.483651 1.483651 0.0;H -1.483651 1.483651 0.0;H -1.483651 -1.483651 0.0;H 1.483651 -1.483651 0.0")
mol = gto.M(atom=CBD, basis='def2-svp', verbose=0)
mf = dft.UKS(mol, xc='pbe0'); mf.grids.level = 5; mf.conv_tol = 1e-11
mf.kernel()
for _ in range(8):                      # rotate along the unstable mode until the UKS solution is stable
    mo, stable = mf.stability(return_status=True)[:2]
    if stable is True: break
    mf.kernel(dm0=mf.make_rdm1(mo, mf.mo_occ))
f = thermo.harmonic_analysis(mol, mf.Hessian().kernel())['freq_wavenumber'].real
here = os.path.dirname(os.path.abspath(__file__))
p = os.path.join(here, 'refs.json')
refs = json.load(open(p)) if os.path.exists(p) else {}
refs['cbd_uks_freq'] = [float(x) for x in f]
json.dump(refs, open(p, 'w'), indent=1)
print('E %.8f  <S^2> %.3f  lowest %.1f' % (mf.e_tot, mf.spin_square()[0], f.min()))
