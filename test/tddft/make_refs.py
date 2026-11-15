# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: AGPL-3.0-or-later

# Regenerates refs.json with pyscf (not needed to run the tests): python3 make_refs.py
import json
from pyscf import gto, dft, qmmm
import numpy as np
EV = 27.211386245988
H2CO = "C 0 0 0;O 0 0 1.22;H 0 0.94 -0.58;H 0.1 -0.94 -0.58"
OH = "O 0 0 0;H 0 0 0.97"
PC_XYZ = np.array([[3.0, 0.0, 2.0], [-3.0, 0.0, 0.0]])   # Angstrom
PC_Q = [0.5, -0.5]
def run(atom, spin, solvent, mode, pc=False):
    mol = gto.M(atom=atom, basis='def2-svp', spin=spin, verbose=0)
    mf = (dft.UKS if spin else dft.RKS)(mol, xc='pbe0')
    if solvent:
        mf = mf.PCM(); mf.with_solvent.method = 'C-PCM'; mf.with_solvent.eps = 78.4
    if pc:
        mf = qmmm.mm_charge(mf, PC_XYZ, PC_Q, unit='Angstrom')
    mf.grids.level = 5; mf.kernel()
    if mode == 'tda':
        td = mf.TDA()
    else:
        td = mf.TDDFT()
    if mode == 'triplet':
        td.singlet = False
    td.nstates = 3
    if solvent:
        td.with_solvent.equilibrium_solvation = False
    td.kernel()
    return [float(x*EV) for x in td.e]
# NOTE: pyscf's UKS TDDFT for OH may or may not return the Goldstone mode (state 1, ~0.02 eV) of the degenerate pi pair;
# when it is missing the list is shifted by one state and run_tests.py compares the wrong states. After regenerating,
# check that 'oh_gas_uks', 'oh_water_uks' and 'oh_gas_uks_pc' start with that state (as in the committed refs.json).
refs = {
  'h2co_gas_singlet':  run(H2CO, 0, False, 'full'),
  'h2co_gas_triplet':  run(H2CO, 0, False, 'triplet'),
  'h2co_gas_tda':      run(H2CO, 0, False, 'tda'),
  'h2co_water_singlet': run(H2CO, 0, True, 'full'),
  'h2co_water_triplet': run(H2CO, 0, True, 'triplet'),
  'oh_gas_uks':        run(OH, 1, False, 'full'),
  'oh_water_uks':      run(OH, 1, True, 'full'),
  'h2co_gas_singlet_pc': run(H2CO, 0, False, 'full', True),
  'oh_gas_uks_pc':     run(OH, 1, False, 'full', True),
  'h2co_water_singlet_pc': run(H2CO, 0, True, 'full', True),
}
json.dump(refs, open('refs.json', 'w'), indent=1)
print(refs)
