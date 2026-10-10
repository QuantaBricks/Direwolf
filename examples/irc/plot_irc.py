#!/usr/bin/env python3
# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Plot an IRC written by Direwolf: energy profile and the distances that change most.

usage: plot_irc.py <base>        reads <base>.irc.dat and <base>.irc.xyz, writes <base>.irc.png
       plot_irc.py <base> -n 4   show the 4 most-changing atom pairs (default 3)

<base>.irc.dat   s (amu^1/2 bohr; negative = backward), E (Ha), E-E(TS) (kcal/mol), gmax
<base>.irc.xyz   one frame per line of the .dat, backward end ... TS ... forward end
"""
import sys, itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_xyz(path):
    lines = open(path).read().split('\n')
    frames, el, i = [], None, 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        rows = [l.split() for l in lines[i + 2:i + 2 + n]]
        el = [r[0] for r in rows]
        frames.append(np.array([[float(x) for x in r[1:4]] for r in rows]))
        i += n + 2
    return el, np.array(frames)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    base = sys.argv[1]
    npair_show = int(sys.argv[sys.argv.index('-n') + 1]) if '-n' in sys.argv else 3
    dat = np.loadtxt(base + '.irc.dat')
    s, erel = dat[:, 0], dat[:, 2]
    el, X = read_xyz(base + '.irc.xyz')
    if len(X) != len(s):
        sys.exit('%s.irc.xyz has %d frames, %s.irc.dat has %d rows' % (base, len(X), base, len(s)))

    pairs = list(itertools.combinations(range(len(el)), 2))
    D = np.array([[np.linalg.norm(x[i] - x[j]) for i, j in pairs] for x in X])
    # the pairs whose distance changes most between the two ends of the path
    change = np.abs(D[-1] - D[0])
    # two hydrogens count only when they come within bonding range somewhere on the path;
    # otherwise the largest changes are non-bonded H...H distances of the spectator groups
    for k, (i, j) in enumerate(pairs):
        if el[i] == 'H' and el[j] == 'H' and D[:, k].min() > 1.3:
            change[k] = 0.0
    order = np.argsort(-change)[:npair_show]

    fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True, gridspec_kw={'height_ratios': [3, 2]})
    ax[0].plot(s, erel, 'o-', ms=3, color='#1f77b4')
    ax[0].axvline(0, color='gray', lw=0.6, ls='--')
    ax[0].set_ylabel('E - E(TS)  (kcal/mol)')
    ax[0].set_title(base.split('/')[-1] + ' IRC')
    for k in order:
        i, j = pairs[k]
        ax[1].plot(s, D[:, k], label='%s%d-%s%d' % (el[i], i + 1, el[j], j + 1))
    ax[1].set_xlabel(r's  (amu$^{1/2}$ bohr)')
    ax[1].set_ylabel('distance (A)')
    ax[1].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.25)
    plt.tight_layout()
    out = base + '.irc.png'
    plt.savefig(out, dpi=140)
    print('written', out)


if __name__ == '__main__':
    main()
