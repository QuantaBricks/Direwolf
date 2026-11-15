# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: AGPL-3.0-or-later

import numpy as np
from pyscf import gto, dft
m=gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692",basis='def2-svp',unit='Angstrom',verbose=0)
def energy(F, pcm):
    mf=dft.RKS(m); mf.xc='pbe0'; mf.grids.level=5; mf.conv_tol=1e-12; mf.max_cycle=100
    if pcm:
        mf=mf.PCM(); mf.with_solvent.method='C-PCM'; mf.with_solvent.eps=78.4
    h=mf.get_hcore()
    r=m.intor_symmetric('int1e_r',comp=3)
    mf.get_hcore=lambda *a,F=F,h=h,r=r: h + np.einsum('x,xij->ij',F,r)   # E = ... + F.r  (electron charge -1 => -(-r)F)
    return mf.kernel()
h=2e-3
for pcm in (False,True):
    e0=energy(np.zeros(3),pcm); a=[]
    for ax in (0,1,2):
        F=np.zeros(3);F[ax]=h
        a.append(-(energy(F,pcm)+energy(-F,pcm)-2*e0)/h**2)
    print('pyscf FD', 'C-PCM' if pcm else 'gas ', np.round(a,4), 'iso', round(sum(a)/3,4))
