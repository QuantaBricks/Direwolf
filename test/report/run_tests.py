#!/usr/bin/env python3
# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: AGPL-3.0-or-later

"""EngineUp report / driver regression: what the calculation hands back in opts%rep and what the driver
puts into <base>.json, plus the driver behaviours the report work touched.
Usage: python3 run_tests.py [path/to/Direwolf] [--slow] [--update]   (run from anywhere; ~2 min, --slow ~6 min more)
  checks    JSON content against the .out of the same run (energies, charges, TDDFT tables, frequencies) and
            against golden values in refs.json (made by --update from the current binary) and pyscf (make_refs.py)
  --slow    also: an interrupted TDDFT call leaves the result energy in the JSON (opts%progress), and an IRC
            with stability = 1 runs the stability analysis
  --update  rewrite the golden part of refs.json from this binary's results (review the diff before committing)"""
import json, os, re, signal, subprocess, sys, tempfile, time
here = os.path.dirname(os.path.abspath(__file__))
args = [a for a in sys.argv[1:] if not a.startswith('--')]
slow = '--slow' in sys.argv
update = '--update' in sys.argv
exe = os.path.abspath(args[0]) if args else os.path.join(here, '..', '..', 'Direwolf')
refs_path = os.path.join(here, 'refs.json')
refs = json.load(open(refs_path)) if os.path.exists(refs_path) else {}
gold = refs.get('golden', {})
new_gold = {}

W = "O 0.0 0.0 0.1173\nH 0.0 0.7572 -0.4692\nH 0.0 -0.7572 -0.4692"
F = "C 0.0 0.0 0.0\nO 0.0 0.0 1.22\nH 0.0 0.94 -0.58\nH 0.1 -0.94 -0.58"
OH = "O 0.0 0.0 0.0\nH 0.0 0.0 0.97"
H2 = "H 0.0 0.0 0.0\nH 0.0 0.0 2.0"
CBD = ("C 0.720000 0.720000 0.0\nC -0.720000 0.720000 0.0\nC -0.720000 -0.720000 0.0\nC 0.720000 -0.720000 0.0\n"
       "H 1.483651 1.483651 0.0\nH -1.483651 1.483651 0.0\nH -1.483651 -1.483651 0.0\nH 1.483651 -1.483651 0.0")

def inp(atoms, mult=1, extra='', J='exact', K='exact', blocks=''):
    return ("&molecule\n imult=%d\n icharge=0\n functional='PBE0'\n baselabel='def2-svp'\n verbose=1\n J='%s'\n K='%s'\n%s&end\n%s"
            "&atoms\n%s\n&end\n" % (mult, J, K, extra, blocks, atoms))

bad = 0
def check(name, ok, msg=''):
    global bad
    bad += not ok
    print('%s %-34s %s' % ('ok  ' if ok else 'FAIL', name, msg))

def run(d, name, text, env=None):
    f = os.path.join(d, name + '.inp')
    open(f, 'w').write(text)
    e = dict(os.environ, **(env or {}))
    p = subprocess.run([exe, f, os.path.join(d, name + '.out')], cwd=d, env=e, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = open(os.path.join(d, name + '.out')).read() if os.path.exists(os.path.join(d, name + '.out')) else ''
    try:
        js = json.load(open(os.path.join(d, name + '.json')))
    except Exception:
        js = None
    return p.returncode, out, js

def out_energy(out):
    m = re.search(r'Total Energy \(Hartree\)\s*=\s*(-?\d+\.\d+)', out) or re.search(r'FINAL SINGLE POINT ENERGY \(Eh\):\s*(-?\d+\.\d+)', out)
    return float(m.group(1)) if m else None

def close(a, b, tol):
    return a is not None and b is not None and abs(a - b) < tol

def golden(name, value, tol, label):
    new_gold[name] = value
    if update:
        check(name, True, '%s = %r (recorded)' % (label, value)); return
    if name not in gold:
        check(name, False, 'no golden value (run --update)'); return
    check(name, close(value, gold[name], tol), '%s %.10g vs %.10g (tol %.0e)' % (label, value, gold[name], tol))

with tempfile.TemporaryDirectory() as d:
    # ---- single point: properties and result ----
    rc, out, js = run(d, 'sp', inp(F, extra=' calc_force=.true.\n'))
    ok = rc == 0 and js is not None and set(js) == {'meta', 'properties', 'result'}
    check('sp keys', ok, str(sorted(js)) if js else 'no JSON')
    if js:
        check('sp energy = .out', close(js['result']['energy'], out_energy(out), 5e-9), '%.10f' % js['result']['energy'])
        q = js['properties']['charges']
        check('sp charges neutral', abs(sum(q['mulliken'])) < 1e-6 and abs(sum(q['lowdin'])) < 1e-6 and len(js['result']['forces']) == 4)
        golden('sp_energy', js['result']['energy'], 1e-7, 'E')

    # ---- TDDFT: closed-shell singlets and triplets, open-shell reference ----
    def td_table(out, header):
        i = out.index(header)
        return [float(r[1]) for r in re.findall(r'^\s+(\d+)\s+(-?\d+\.\d{4})\s+-?\d+\.\d{4}\s+\S+\s+\S+', out[i:], re.M)[:3]]
    rc, out, js = run(d, 'td_rks', inp(F, extra=' tddft_nroots=3\n tddft_triplets=.true.\n calc_force=.false.\n'))
    ok = rc == 0 and js is not None and {'tddft', 'tddft_triplet', 'result', 'properties'} <= set(js)
    check('tddft rks keys', ok)
    if ok:
        es = [s['energy_ev'] for s in js['tddft']['states']]; et = [s['energy_ev'] for s in js['tddft_triplet']['states']]
        check('tddft singlets = .out', len(es) == 3 and all(close(a, b, 6e-5) for a, b in zip(es, td_table(out, 'Singlet excited states'))), str([round(x, 4) for x in es]))
        check('tddft triplets = .out', len(et) == 3 and all(close(a, b, 6e-5) for a, b in zip(et, td_table(out, 'Triplet excited states'))), str([round(x, 4) for x in et]))
        check('tddft transitions listed', all(len(s['transitions']) > 0 for s in js['tddft']['states']))
        golden('td_rks_s1', es[0], 1e-5, 'S1 (eV)'); golden('td_rks_t1', et[0], 1e-5, 'T1 (eV)')
    rc, out, js = run(d, 'td_uks', inp(OH, mult=2, extra=' tddft_nroots=3\n calc_force=.false.\n'))
    ok = rc == 0 and js is not None and 'tddft' in js and js['tddft'].get('unrestricted')
    check('tddft uks keys', bool(ok))
    if ok:
        check('tddft uks spins', all('spin' in t for s in js['tddft']['states'] for t in s['transitions']))
        golden('td_uks_s2', js['tddft']['states'][1]['energy_ev'], 1e-5, 'state 2 (eV)')

    # ---- polarizability, static and dynamic ----
    rc, out, js = run(d, 'polar', inp(W, extra=" runtype='polar'\n polar_omega='0.05, 0.1'\n calc_force=.false.\n"))
    ok = rc == 0 and js is not None and 'polarizability' in js
    check('polar keys', bool(ok))
    if ok:
        p = js['polarizability']; t = p['tensor_au']
        check('polar tensor symmetric', all(close(t[i][j], t[j][i], 1e-7) for i in range(3) for j in range(3)))
        check('polar dynamic records', len(p.get('dynamic', [])) == 2)
        check('polar energy = .out', close(js['result']['energy'], out_energy(out), 5e-9) if out_energy(out) else True)
        golden('polar_iso', p['isotropic_au'], 1e-5, 'alpha_iso (a.u.)'); golden('polar_iso_w05', p['dynamic'][0]['isotropic_au'], 1e-5, 'alpha_iso(0.05)')

    # ---- frequencies: the process ends in the driver (opts%done), JSON from the report ----
    rc, out, js = run(d, 'freq', inp(W, extra=' freq=.true.\n calc_force=.false.\n'))
    ok = rc == 0 and js is not None and {'freq', 'thermo', 'dipole_derivatives', 'result'} <= set(js)
    check('freq keys', bool(ok), str(sorted(js)) if js else 'no JSON')
    if ok:
        fr = [m['frequency'] for m in js['freq']]
        fo = [float(x) for x in re.search(r'Frequency \(cm-1\)\s+(.*)', out).group(1).split()[:3]]
        check('freq = .out', len(fr) == 3 and all(close(a, b, 6e-5) for a, b in zip(fr, fo)), str([round(x, 2) for x in fr]))
        check('freq energy = .out', close(js['result']['energy'], out_energy(out), 5e-9))
        check('freq no results block', 'Hcore trace' not in out and '[THERMOEND]' in out)
        check('freq dipole derivatives', len(js['dipole_derivatives']) == 3)
        golden('freq_1', fr[0], 0.02, 'lowest (cm-1)')
    rc, out, js = run(d, 'oh_freq', inp(OH, mult=2, extra=' freq=.true.\n calc_force=.false.\n'))
    check('freq open shell', rc == 0 and js is not None and len(js.get('freq', [])) == 1 and 'thermo' in js)
    if js and js.get('freq'):
        golden('oh_freq', js['freq'][0]['frequency'], 0.02, 'OH (cm-1)')

    # ---- stability = 2 reaches the optimizer; the Hessian of the unrestricted solution is the open-shell one ----
    rc, out, js = run(d, 'h2_opt', inp(H2, extra=" runtype='opt, freq'\n stability=2\n"))
    ok = rc == 0 and js is not None and 'opt' in js and 'freq' in js
    check('opt+freq keys', bool(ok))
    if ok:
        c1 = js['opt']['cycles'][0]['energy']
        check('stability=2 in opt: step 1 unrestricted', c1 < -1.0100, 'E1 %.8f (restricted -0.98474, unrestricted -1.01118)' % c1)
        z = [g[2] for g in js['result']['geometry']]
        check('opt converged, H2 bond 0.7607', js['opt'].get('converged') and abs((z[1] - z[0]) - 0.7607) < 2e-3, '%.4f' % (z[1] - z[0]))
        golden('h2_opt_energy', js['result']['energy'], 1e-6, 'E')
    rc, out, js = run(d, 'cbd', inp(CBD, J='RI', K='cosx', extra=" runtype='freq'\n stability=2\n freq=.true.\n"))
    ok = rc == 0 and js is not None and len(js.get('freq', [])) == 18
    check('cbd uks hessian: 18 modes', bool(ok))
    if ok:
        fr = [m['frequency'] for m in js['freq']]
        check('cbd uks hessian: all real', min(fr) > 0.0, 'lowest %.1f cm-1 (closed-shell code gave -3102)' % min(fr))
        if 'cbd_uks_freq' in refs:
            err = max(abs(a - b) for a, b in zip(sorted(fr), sorted(refs['cbd_uks_freq'])))
            check('cbd uks hessian vs pyscf', err < 3.0, 'max|df| %.1f cm-1 (pyscf, grid level 5; tol 3)' % err)
        golden('cbd_lowest', min(fr), 0.5, 'lowest (cm-1)')

    # ---- an explicit def2 ECP name is the one def2 ECP (it used to run all-electron without a word) ----
    CH3I = "C 0.0 0.0 0.0\nI 0.0 0.0 2.14\nH 1.03 0.0 -0.37\nH -0.515 0.892 -0.37\nH -0.515 -0.892 -0.37"
    e_ecp = {}
    for tag, extra in (('default', ' calc_force=.false.\n'), ('explicit', " ecplabel='def2-svp'\n calc_force=.false.\n")):
        rc, out, js = run(d, 'ecp_' + tag, inp(CH3I, J='RI', K='cosx', extra=extra))
        e_ecp[tag] = js['result']['energy'] if js else None
    check('ecplabel=def2-svp is the def2 ECP', e_ecp['default'] is not None and close(e_ecp['default'], e_ecp['explicit'], 1e-9)
          and e_ecp['default'] > -400.0, 'E %s vs %s (all-electron would be -2039)' % (e_ecp['default'], e_ecp['explicit']))

    # ---- Molden file: the [Atoms] charge column is the effective nuclear charge (Z minus ECP core electrons) ----
    rc, out, js = run(d, 'molden_ecp', inp(CH3I, J='RI', K='cosx', extra=' calc_force=.false.\n',
                                           blocks="&molden\n molden_write=.true.\n molden_file='molden_ecp.molden'\n&end\n"))
    zcol = {}
    mf = os.path.join(d, 'molden_ecp.molden')
    if os.path.exists(mf):
        sec = False
        for ln in open(mf):
            if ln.startswith('['): sec = ln.startswith('[Atoms]'); continue
            t = ln.split()
            if sec and len(t) >= 6: zcol[t[0]] = int(t[2])
    check('molden: effective charge column', zcol == {'C': 6, 'I': 25, 'H': 1}, str(zcol))

    # ---- cube files: density integrates to the electron count, the potential is the exact one without RI-J ----
    for tag, J, cu in (('exact', 'exact', 'den, esp'), ('ri', 'RI', 'den, esp'), ('ri_exact', 'RI', 'den, esp_exact')):
        K = 'exact' if J == 'exact' else 'cosx'
        rc, out, js = run(d, 'cube_' + tag, inp(W, J=J, K=K, extra=" cube='%s'\n calc_force=.false.\n" % cu))
        vals = {}
        for kind in ('den', 'esp'):
            f = os.path.join(d, 'cube_%s.%s.cube' % (tag, kind))
            if not os.path.exists(f): continue
            L = open(f).read().split('\n')
            na = int(L[2].split()[0]); h = float(L[3].split()[1])
            vals[kind] = (h, [float(x) for x in ' '.join(L[6 + na:]).split()])
        ok = set(vals) == {'den', 'esp'}
        if ok:
            ne = sum(vals['den'][1]) * vals['den'][0] ** 3
            check('cube density integral (%s)' % tag, abs(ne - 10.0) < 0.05, 'N = %.4f' % ne)
            check('cube potential range (%s)' % tag, min(vals['esp'][1]) < 0.0 < max(vals['esp'][1]) and max(vals['esp'][1]) < 100.0,
                  '%.4f .. %.2f' % (min(vals['esp'][1]), max(vals['esp'][1])))
        else:
            check('cube files (%s)' % tag, False, 'missing')
        check('cube potential source (%s)' % tag, ('RI-J fitted density' in out) == (tag == 'ri'), 'fitted-density note %s' % ('RI-J fitted density' in out))

    # ---- TDDFT hole / electron densities: normalisations, the transition density's dipole, the measures in the JSON ----
    rc, out, js = run(d, 'hole', inp(W, J='RI', K='cosx', extra=" tddft_nroots=2\n tddft_cube='diff:1, hole:1, elec:1, trans:1'\n calc_force=.false.\n"))
    def cube_read(f):
        L = open(f).read().split('\n')
        na = int(L[2].split()[0]); h = float(L[3].split()[1]); o = [float(x) for x in L[2].split()[1:4]]
        n = [int(L[3 + i].split()[0]) for i in range(3)]
        return h, o, n, [float(x) for x in ' '.join(L[6 + na:]).split()]
    try:
        files = {k: cube_read(os.path.join(d, 'hole.s1.%s.cube' % k)) for k in ('diff', 'hole', 'elec', 'trans')}
    except Exception:
        files = None
    if files and js and 'hole_electron' in js['tddft']['states'][0]:
        h = files['hole'][0]; w = h ** 3
        nh = sum(files['hole'][3]) * w; ne = sum(files['elec'][3]) * w; nd = sum(files['diff'][3]) * w
        check('hole/electron integrals', abs(nh - ne) < 5e-3 and 1.0 <= nh < 1.3 and abs(nd) < 5e-3, 'hole %.4f electron %.4f diff %.4f' % (nh, ne, nd))
        h, o, n, v = files['trans']
        mom = [0.0, 0.0, 0.0]
        for ix in range(n[0]):
            for iy in range(n[1]):
                base = (ix * n[1] + iy) * n[2]
                for iz in range(n[2]):
                    t = v[base + iz]
                    if t != 0.0:
                        mom[0] += t * (o[0] + h * ix); mom[1] += t * (o[1] + h * iy); mom[2] += t * (o[2] + h * iz)
        mom = [m * h ** 3 for m in mom]
        td = js['tddft']['states'][0]['transition_dipole_au']
        check('transition density moment = transition dipole', max(abs(a - b) for a, b in zip(mom, td)) < 3e-3, 'moment %s vs %s' % ([round(m, 4) for m in mom], [round(x, 4) for x in td]))
        he = js['tddft']['states'][0]['hole_electron']
        check('hole/electron measures', 0.0 < he['lambda'] < 1.0 and 0.0 <= he['d_ct_angstrom'] < 3.0 and 0.0 < he['q_ct'] < 1.5, str({k: round(x, 3) for k, x in he.items()}))
    else:
        check('hole/electron cubes', False, 'missing cube files or JSON fields')

    # ---- the same for an unrestricted reference (CH3 radical): sums over both spins, no sqrt(2) in the transition density ----
    CH3 = "C 0.0 0.0 0.0\nH 1.07 0.0 0.0\nH -0.535 0.9267 0.0\nH -0.535 -0.9267 0.0"
    rc, out, js = run(d, 'holeu', inp(CH3, mult=2, J='RI', K='cosx', extra=" tddft_nroots=3\n tddft_cube='hole:3, elec:3, trans:3'\n calc_force=.false.\n"))
    try:
        fu = {k: cube_read(os.path.join(d, 'holeu.s3.%s.cube' % k)) for k in ('hole', 'elec', 'trans')}
    except Exception:
        fu = None
    if fu and js and 'hole_electron' in js['tddft']['states'][2]:
        w = fu['hole'][0] ** 3
        nh = sum(fu['hole'][3]) * w; ne = sum(fu['elec'][3]) * w
        check('hole/electron integrals (unrestricted)', abs(nh - ne) < 5e-3 and 1.0 <= nh < 1.3, 'hole %.4f electron %.4f' % (nh, ne))
        h, o, n, v = fu['trans']
        mom = [0.0, 0.0, 0.0]
        for ix in range(n[0]):
            for iy in range(n[1]):
                base = (ix * n[1] + iy) * n[2]
                for iz in range(n[2]):
                    t = v[base + iz]
                    if t != 0.0:
                        mom[0] += t * (o[0] + h * ix); mom[1] += t * (o[1] + h * iy); mom[2] += t * (o[2] + h * iz)
        mom = [m * h ** 3 for m in mom]
        td = js['tddft']['states'][2]['transition_dipole_au']
        check('transition density moment = transition dipole (unrestricted)', max(abs(a - b) for a, b in zip(mom, td)) < 3e-3,
              'moment %s vs %s' % ([round(m, 4) for m in mom], [round(x, 4) for x in td]))
    else:
        check('hole/electron cubes (unrestricted)', False, 'missing cube files or JSON fields')

    # ---- NMR: GIAO shielding of water, Hartree-Fock / def2-SVP (reference values from a published program: O 347.70, H 31.18 ppm) ----
    rc, out, js = run(d, 'nmr_w', inp(W, J='RI', K='cosx', extra=" functional='HF'\n runtype='nmr'\n").replace("functional='PBE0'\n", ""))
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr water O/H shielding', abs(at[0]['isotropic'] - 347.70) < 0.05 and abs(at[1]['isotropic'] - 31.18) < 0.02
              and abs(at[1]['isotropic'] - at[2]['isotropic']) < 1e-4, 'O %.3f H %.3f %.3f' % (at[0]['isotropic'], at[1]['isotropic'], at[2]['isotropic']))
        tr = js['nmr']['atoms'][0]['tensor']
        check('nmr tensor from the JSON', abs(sum(tr[i][i] for i in range(3)) / 3 - at[0]['isotropic']) < 1e-9, 'trace/3 = isotropic')
    else:
        check('nmr water', False, 'no nmr record in the JSON')

    # ---- NMR with a density functional (PBE / def2-SVP; reference values from a published program: O 343.968, H 31.935 ppm) ----
    rc, out, js = run(d, 'nmr_w_pbe', inp(W, J='RI', K='cosx', extra=" functional='PBE'\n runtype='nmr'\n").replace("functional='PBE0'\n", ""))
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr water PBE shielding', abs(at[0]['isotropic'] - 343.968) < 0.02 and abs(at[1]['isotropic'] - 31.935) < 0.01,
              'O %.3f H %.3f' % (at[0]['isotropic'], at[1]['isotropic']))
    else:
        check('nmr water PBE', False, 'no nmr record in the JSON')

    # ---- NMR with a range-separated hybrid (CAM-B3LYP / def2-SVP; reference values from a published program: O 349.360, H 31.533 ppm) ----
    rc, out, js = run(d, 'nmr_w_cam', inp(W, J='RI', K='cosx', extra=" functional='CAM-B3LYP'\n runtype='nmr'\n").replace("functional='PBE0'\n", ""))
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr water CAM-B3LYP shielding', abs(at[0]['isotropic'] - 349.360) < 0.02 and abs(at[1]['isotropic'] - 31.533) < 0.01,
              'O %.3f H %.3f' % (at[0]['isotropic'], at[1]['isotropic']))
    else:
        check('nmr water CAM-B3LYP', False, 'no nmr record in the JSON')

    # ---- NMR with a meta-GGA (TPSS / def2-SVP; reference values from a published program: O 341.324, H 31.994 ppm) ----
    rc, out, js = run(d, 'nmr_w_tpss', inp(W, J='RI', K='cosx', extra=" functional='TPSS'\n runtype='nmr'\n").replace("functional='PBE0'\n", ""))
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr water TPSS shielding', abs(at[0]['isotropic'] - 341.324) < 0.03 and abs(at[1]['isotropic'] - 31.994) < 0.01,
              'O %.3f H %.3f' % (at[0]['isotropic'], at[1]['isotropic']))
    else:
        check('nmr water TPSS', False, 'no nmr record in the JSON')

    # ---- NMR environment: point charges, C-PCM water, an ECP (reference values from a published program, def2-SVP) ----
    pcblk = "&pointcharges\n npc=2\n pc_q=0.4,-0.4\n pc_x=0.0,2.0\n pc_y=0.0,1.0\n pc_z=3.0,-2.0\n&end\n"
    rc, out, js = run(d, 'nmr_pc', inp(W, extra=" functional='PBE'\n runtype='nmr'\n").replace("functional='PBE0'\n", "") + pcblk)
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr water with point charges', abs(at[0]['isotropic'] - 348.037) < 0.02 and abs(at[1]['isotropic'] - 30.884) < 0.01
              and abs(at[2]['isotropic'] - 31.605) < 0.01, 'O %.3f H %.3f %.3f' % (at[0]['isotropic'], at[1]['isotropic'], at[2]['isotropic']))
    else:
        check('nmr water with point charges', False, 'no nmr record in the JSON')
    rc, out, js = run(d, 'nmr_cpcm', inp(W, extra=" functional='PBE'\n runtype='nmr'\n").replace("functional='PBE0'\n", "")
                      + "&cosmo\n solvent='water'\n convention='cpcm'\n&end\n")
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr water in C-PCM water', abs(at[0]['isotropic'] - 352.187) < 0.02 and abs(at[1]['isotropic'] - 31.112) < 0.01,
              'O %.3f H %.3f' % (at[0]['isotropic'], at[1]['isotropic']))
    else:
        check('nmr water in C-PCM water', False, 'no nmr record in the JSON')
    CH3I = "C 0.0 0.0 0.0\nI 0.0 0.0 2.14\nH 1.03 0.0 -0.36\nH -0.515 0.892 -0.36\nH -0.515 -0.892 -0.36"
    rc, out, js = run(d, 'nmr_ecp', inp(CH3I, extra=" functional='HF'\n runtype='nmr'\n").replace("functional='PBE0'\n", ""))
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr CH3I with the iodine ECP', abs(at[0]['isotropic'] - 195.199) < 0.01 and abs(at[1]['isotropic'] - 460.660) < 0.02
              and abs(at[2]['isotropic'] - 30.155) < 0.01, 'C %.3f I %.3f H %.3f' % (at[0]['isotropic'], at[1]['isotropic'], at[2]['isotropic']))
    else:
        check('nmr CH3I with the iodine ECP', False, 'no nmr record in the JSON')

    # ---- NMR grid: CH3+ CAM-B3LYP carbon, exact integrals; the SCF grid alone leaves 0.019 ppm (reference: a published program, and Direwolf on a 307k-point grid, -203.817) ----
    CH3P = "C 0.0 0.0 0.0\nH 1.079 0.0 0.0\nH -0.5395 0.934441 0.0\nH -0.5395 -0.934441 0.05"
    rc, out, js = run(d, 'nmr_ch3p', inp(CH3P, extra=" functional='CAM-B3LYP'\n runtype='nmr'\n icharge=1\n").replace("functional='PBE0'\n", "").replace("icharge=0\n", ""))
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr CH3+ CAM-B3LYP (NMR grid)', abs(at[0]['isotropic'] + 203.817) < 0.006, 'C %.3f' % at[0]['isotropic'])
    else:
        check('nmr CH3+ CAM-B3LYP (NMR grid)', False, 'no nmr record in the JSON')
    # ---- open-shell NMR: methyl radical (doublet), def2-SVP; reference: UHF C 127.916, H 28.094 / 28.074; UKS PBE0 RIJCOSX C 108.833 ----
    CH3R = "C 0.0 0.0 0.0\nH 1.079 0.0 0.0\nH -0.5395 0.934441 0.0\nH -0.5395 -0.934441 0.05"
    rc, out, js = run(d, 'nmr_uhf', inp(CH3R, mult=2, extra=" functional='HF'\n runtype='nmr'\n").replace("functional='PBE0'\n", ""))
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr methyl radical UHF', abs(at[0]['isotropic'] - 127.916) < 0.01 and abs(at[1]['isotropic'] - 28.094) < 0.002
              and abs(at[3]['isotropic'] - 28.074) < 0.002, 'C %.3f H %.3f %.3f' % (at[0]['isotropic'], at[1]['isotropic'], at[3]['isotropic']))
    else:
        check('nmr methyl radical UHF', False, 'no nmr record in the JSON')
    rc, out, js = run(d, 'nmr_uks', inp(CH3R, mult=2, J='RI', K='cosx', extra=" runtype='nmr'\n"))
    if js and 'nmr' in js:
        at = js['nmr']['atoms']
        check('nmr methyl radical UKS PBE0', abs(at[0]['isotropic'] - 108.833) < 0.01 and abs(at[1]['isotropic'] - 27.426) < 0.002,
              'C %.3f H %.3f' % (at[0]['isotropic'], at[1]['isotropic']))
    else:
        check('nmr methyl radical UKS PBE0', False, 'no nmr record in the JSON')

    # ---- solvent convention: C-PCM (f = (eps-1)/eps, H radius 1.10) solvates water 8-18 % more strongly than the default ----
    esol = {}
    for tag, ex in (('gas', ''), ('cosmo', "&cosmo\n solvent='water'\n&end\n"), ('cpcm', "&cosmo\n solvent='water'\n convention='cpcm'\n&end\n")):
        rc, out, js = run(d, 'conv_' + tag, inp(W, J='RI', K='cosx', extra=' calc_force=.false.\n', blocks=ex))
        esol[tag] = out_energy(out)
    ok = all(v is not None for v in esol.values())
    ratio = (esol['cpcm'] - esol['gas']) / (esol['cosmo'] - esol['gas']) if ok else None
    check('cosmo convention cpcm', ok and 1.08 < ratio < 1.20, 'ratio of solvation energies %s' % ('%.3f' % ratio if ok else 'n/a'))

    # ---- thermochemistry conditions: a 1 mol/L standard state raises G by RT ln(24.4654) = 1.894 kcal/mol at 298.15 K ----
    gcorr = {}
    for tag, ex in (('gas', ''), ('conc', ' freq_conc=1.0\n')):
        rc, out, js = run(d, 'thermo_' + tag, inp(W, J='RI', K='cosx', extra=" runtype='freq'\n" + ex))
        m = re.search(r'Thermal correction to G =\s*(-?\d+\.\d+)', out)
        gcorr[tag] = float(m.group(1)) if m else None
    ok = gcorr['gas'] is not None and gcorr['conc'] is not None
    check('thermo: 1 mol/L standard state', ok and abs((gcorr['conc'] - gcorr['gas']) * 627.509 - 1.894) < 0.01,
          'dG %.3f kcal/mol' % ((gcorr['conc'] - gcorr['gas']) * 627.509) if ok else 'no THERMO block')

    # ---- optimization: TDDFT and cube files belong to the last frame only ----
    WOPT = "O 0.0 0.0 0.15\nH 0.0 0.77 -0.48\nH 0.0 -0.77 -0.48"
    rc, out, js = run(d, 'opt_last', inp(WOPT, J='RI', K='cosx', extra=" runtype='opt'\n tddft_nroots=2\n cube='den'\n"))
    check('opt + TDDFT + cube: one TDDFT block', rc == 0 and out.count('[TDDFT]') == 1 and '[LAST FRAME]' in out,
          'blocks %d' % out.count('[TDDFT]'))
    mf = os.path.join(d, 'opt_last.den.cube')
    if os.path.exists(mf):
        zc = float(open(mf).read().split('\n')[6].split()[4])
        xyz = [l.split() for l in open(os.path.join(d, 'opt_last.opt.xyz')).read().strip().split('\n')]
        zopt = [float(l[3]) for l in xyz[-3:] if l[0] == 'O']
        check('opt cube is the last geometry', bool(zopt) and abs(zc * 0.529177210903 - zopt[0]) < 1e-4, '%.5f vs %s' % (zc * 0.529177210903, zopt))
    else:
        check('opt cube is the last geometry', False, 'no cube file')

    # ---- cube = '..., all' in an optimization: one set of files per cycle and the last frame ----
    rc, out, js = run(d, 'opt_all', inp(WOPT, J='RI', K='cosx', extra=" runtype='opt'\n cube='den, all'\n"))
    cyc_files = sorted(f for f in os.listdir(d) if re.match(r'opt_all\.c\d+\.den\.cube$', f))
    ncyc = len(re.findall(r'\[OPTCYCLE \d+\]', out))
    check('cube all in opt', rc == 0 and len(cyc_files) == ncyc >= 2 and os.path.exists(os.path.join(d, 'opt_all.den.cube')),
          '%d cycle files for %d cycles' % (len(cyc_files), ncyc))

    # ---- slow: progress callback and IRC ----
    if slow:
        caf = open(os.path.join(here, '..', '..', 'examples', 'drug_molecules', 'caffeine_wb97x3c.inp')).read()
        atoms = caf.split('&atoms')[1].split('&end')[0].strip()
        f = os.path.join(d, 'caf.inp')
        open(f, 'w').write(inp(atoms, J='RI', K='cosx', extra=' tddft_nroots=6\n calc_force=.false.\n n_threads=8\n'))
        p = subprocess.Popen([exe, f, os.path.join(d, 'caf.out')], cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        jp = os.path.join(d, 'caf.json'); t0 = time.time(); got = False
        while time.time() - t0 < 600 and p.poll() is None:
            try:
                if '"result"' in open(jp).read(): got = True; break
            except OSError:
                pass
            time.sleep(1)
        time.sleep(2)
        alive = p.poll() is None
        p.kill(); p.wait()
        js = json.load(open(jp)) if os.path.exists(jp) else {}
        check('progress: result energy before TDDFT ends', got and alive and 'result' in js and 'tddft' not in js,
              'E %.6f' % js['result']['energy'] if js.get('result') else 'no result')
        irc = open(os.path.join(here, '..', '..', 'examples', 'irc', 'hcn_hnc_irc.inp')).read().replace(" K          = 'cosx'", " K          = 'cosx'\n stability  = 1")
        rc, out, js = run(d, 'irc', irc)
        check('irc stability reaches the driver', rc == 0 and out.count('STABILITY') > 10 and js is not None and len(js['irc']['points']) > 10,
              '%d STABILITY lines, %d points' % (out.count('STABILITY'), len(js['irc']['points']) if js else 0))

if update:
    refs['golden'] = new_gold
    json.dump(refs, open(refs_path, 'w'), indent=1)
    print('refs.json golden values rewritten')
sys.exit(1 if bad else 0)
