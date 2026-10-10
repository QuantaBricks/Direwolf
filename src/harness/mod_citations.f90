! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Literature of the methods a Direwolf run actually used.

module mod_citations
use mod_engine_input_types, only: engine_input_t
use mod_xc, only: xc_write_references, xc_uses_vv10
implicit none
private
public :: write_citations

integer, parameter :: MAXREF = 80
type :: ref_t
   character(len=16)  :: key = ''
   character(len=60)  :: cat = ''
   character(len=300) :: authors = ''
   character(len=300) :: title = ''
   character(len=60)  :: journal = ''
   character(len=8)   :: volume = ''
   character(len=16)  :: pages = ''
   character(len=6)   :: year = ''
   logical            :: used = .false.
end type ref_t
type(ref_t), save :: tab(MAXREF)
integer, save :: nref = 0

contains

subroutine add(key, cat, authors, title, journal, volume, pages, year)
character(len=*), intent(in) :: key, cat, authors, title, journal, volume, pages, year
nref = nref + 1
tab(nref)%key = key ; tab(nref)%cat = cat ; tab(nref)%authors = authors
tab(nref)%title = title ; tab(nref)%journal = journal
tab(nref)%volume = volume ; tab(nref)%pages = pages ; tab(nref)%year = year
tab(nref)%used = .false.
end subroutine add

subroutine build_table()
if (nref > 0) return
call add('libcint', 'Integrals', 'Sun, Q.', &
     'Libcint: An efficient general integral library for Gaussian basis functions', &
     'J. Comput. Chem.', '36', '1664', '2015')
call add('becke88g', 'Numerical integration grid', 'Becke, A. D.', &
     'A multicenter numerical integration scheme for polyatomic molecules', 'J. Chem. Phys.', '88', '2547', '1988')
call add('lebedev', 'Numerical integration grid', 'Lebedev, V. I. and Laikov, D. N.', &
     'A quadrature formula for the sphere of the 131st algebraic order of accuracy', 'Dokl. Math.', '59', '477', '1999')
call add('treutler', 'Numerical integration grid', 'Treutler, O. and Ahlrichs, R.', &
     'Efficient molecular numerical integration schemes', 'J. Chem. Phys.', '102', '346', '1995')
call add('cosx1', 'Semi-numerical exchange (COSX)', 'Neese, F. and Wennmohs, F. and Hansen, A. and Becker, U.', &
     'Efficient, approximate and parallel Hartree-Fock and hybrid DFT calculations. A chain-of-spheres '// &
     'algorithm for the Hartree-Fock exchange', 'Chem. Phys.', '356', '98', '2009')
call add('cosx2', 'Semi-numerical exchange (COSX)', 'Izsák, R. and Neese, F.', &
     'An overlap fitted chain of spheres exchange method', 'J. Chem. Phys.', '135', '144105', '2011')
call add('cosx3', 'Semi-numerical exchange (COSX)', 'Helmich-Paris, B. and de Souza, B. and Neese, F. and Izsák, R.', &
     'An improved chain of spheres for exchange algorithm', 'J. Chem. Phys.', '155', '104109', '2021')
call add('rij1', 'Density fitting (RI-J)', 'Whitten, J. L.', 'Coulombic potential energy integrals and approximations', &
     'J. Chem. Phys.', '58', '4496', '1973')
call add('rij2', 'Density fitting (RI-J)', 'Dunlap, B. I. and Connolly, J. W. D. and Sabin, J. R.', &
     'On some approximations in applications of Xα theory', 'J. Chem. Phys.', '71', '3396', '1979')
call add('rij3', 'Density fitting (RI-J)', 'Vahtras, O. and Almlöf, J. and Feyereisen, M. W.', &
     'Integral approximations for LCAO-SCF calculations', 'Chem. Phys. Lett.', '213', '514', '1993')
call add('rij4', 'Density fitting (RI-J)', 'Eichkorn, K. and Treutler, O. and Öhm, H. and Häser, M. and Ahlrichs, R.', &
     'Auxiliary basis sets to approximate Coulomb potentials', 'Chem. Phys. Lett.', '240', '283', '1995')
call add('rik', 'Density fitting (RI-K)', 'Weigend, F.', &
     'A fully direct RI-HF algorithm: Implementation, optimised auxiliary basis sets, demonstration of accuracy '// &
     'and efficiency', 'Phys. Chem. Chem. Phys.', '4', '4285', '2002')
call add('def2', 'Basis sets', 'Weigend, F. and Ahlrichs, R.', &
     'Balanced basis sets of split valence, triple zeta valence and quadruple zeta valence quality for H to Rn: '// &
     'Design and assessment of accuracy', 'Phys. Chem. Chem. Phys.', '7', '3297', '2005')
call add('def2j', 'Basis sets', 'Weigend, F.', 'Accurate Coulomb-fitting basis sets for H to Rn', &
     'Phys. Chem. Chem. Phys.', '8', '1057', '2006')
call add('def2jk', 'Basis sets', 'Weigend, F.', 'Hartree-Fock exchange fitting basis sets for H to Rn', &
     'J. Comput. Chem.', '29', '167', '2008')
call add('d2', 'Dispersion', 'Grimme, S.', &
     'Semiempirical GGA-type density functional constructed with a long-range dispersion correction', &
     'J. Comput. Chem.', '27', '1787', '2006')
call add('d3', 'Dispersion', 'Grimme, S. and Antony, J. and Ehrlich, S. and Krieg, H.', &
     'A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) '// &
     'for the 94 elements H-Pu', 'J. Chem. Phys.', '132', '154104', '2010')
call add('d3bj', 'Dispersion', 'Grimme, S. and Ehrlich, S. and Goerigk, L.', &
     'Effect of the damping function in dispersion corrected density functional theory', &
     'J. Comput. Chem.', '32', '1456', '2011')
call add('d4a', 'Dispersion', 'Caldeweyher, E. and Bannwarth, C. and Grimme, S.', &
     'Extension of the D3 dispersion coefficient model', 'J. Chem. Phys.', '147', '034112', '2017')
call add('d4b', 'Dispersion', 'Caldeweyher, E. and Ehlert, S. and Hansen, A. and Neugebauer, H. and Spicher, S. '// &
     'and Bannwarth, C. and Grimme, S.', 'A generally applicable atomic-charge dependent London dispersion correction', &
     'J. Chem. Phys.', '150', '154122', '2019')
call add('gcp', 'Dispersion', 'Kruse, H. and Grimme, S.', &
     'A geometrical correction for the inter- and intra-molecular basis set superposition error in Hartree-Fock '// &
     'and density functional theory calculations for large systems', 'J. Chem. Phys.', '136', '154101', '2012')
call add('wb97xd', 'Dispersion', 'Chai, J.-D. and Head-Gordon, M.', &
     'Long-range corrected hybrid density functionals with damped atom-atom dispersion corrections', &
     'Phys. Chem. Chem. Phys.', '10', '6615', '2008')
call add('vv10', 'Non-local correlation', 'Vydrov, O. A. and Van Voorhis, T.', &
     'Nonlocal van der Waals density functional: The simpler the better', 'J. Chem. Phys.', '133', '244103', '2010')
call add('wb97mv', 'Non-local correlation', 'Mardirossian, N. and Head-Gordon, M.', &
     'ωB97M-V: A combinatorially optimized, range-separated hybrid, meta-GGA density functional with VV10 '// &
     'nonlocal correlation', 'J. Chem. Phys.', '144', '214110', '2016')
call add('wb97xv', 'Non-local correlation', 'Mardirossian, N. and Head-Gordon, M.', &
     'ωB97X-V: A 10-parameter, range-separated hybrid, generalized gradient approximation density functional with '// &
     'nonlocal correlation, designed by a survival-of-the-fittest strategy', 'Phys. Chem. Chem. Phys.', '16', '9904', '2014')
call add('b97mv', 'Non-local correlation', 'Mardirossian, N. and Head-Gordon, M.', &
     'Mapping the genome of meta-generalized gradient approximation density functionals: The search for B97M-V', &
     'J. Chem. Phys.', '142', '074111', '2015')
call add('pbeh3c', 'Composite -3c methods', 'Grimme, S. and Brandenburg, J. G. and Bannwarth, C. and Hansen, A.', &
     'Consistent structures and interactions by density functional theory with small atomic orbital basis sets', &
     'J. Chem. Phys.', '143', '054107', '2015')
call add('b973c', 'Composite -3c methods', 'Brandenburg, J. G. and Bannwarth, C. and Hansen, A. and Grimme, S.', &
     'B97-3c: A revised low-cost variant of the B97-D density functional method', 'J. Chem. Phys.', '148', '064104', '2018')
call add('r2scan3c', 'Composite -3c methods', 'Grimme, S. and Hansen, A. and Ehlert, S. and Mewes, J.-M.', &
     'r2SCAN-3c: A "Swiss army knife" composite electronic-structure method', 'J. Chem. Phys.', '154', '064103', '2021')
call add('wb97x3c', 'Composite -3c methods', 'Müller, M. and Hansen, A. and Grimme, S.', &
     'ωB97X-3c: A composite range-separated hybrid DFT method with a molecule-optimized polarized valence '// &
     'double-ζ basis set', 'J. Chem. Phys.', '158', '014103', '2023')
call add('cosmo', 'Solvation', 'Klamt, A. and Schüürmann, G.', &
     'COSMO: a new approach to dielectric screening in solvents with explicit expressions for the screening '// &
     'energy and its gradient', 'J. Chem. Soc., Perkin Trans. 2', '', '799', '1993')
call add('smd', 'Solvation', 'Marenich, A. V. and Cramer, C. J. and Truhlar, D. G.', &
     'Universal solvation model based on solute electron density and on a continuum model of the solvent defined '// &
     'by the bulk dielectric constant and atomic surface tensions', 'J. Phys. Chem. B', '113', '6378', '2009')
call add('diis', 'SCF convergence', 'Pulay, P.', 'Convergence acceleration of iterative sequences. The case of SCF iteration', &
     'Chem. Phys. Lett.', '73', '393', '1980')
call add('adiis', 'SCF convergence', 'Hu, X. and Yang, W.', &
     'Accelerating self-consistent field convergence with the augmented Roothaan-Hall energy function', &
     'J. Chem. Phys.', '132', '054109', '2010')
call add('rfo', 'Geometry optimization', 'Banerjee, A. and Adams, N. and Simons, J. and Shepard, R.', &
     'Search for stationary points on surfaces', 'J. Phys. Chem.', '89', '52', '1985')
call add('ric', 'Geometry optimization', 'Peng, C. and Ayala, P. Y. and Schlegel, H. B. and Frisch, M. J.', &
     'Using redundant internal coordinates to optimize equilibrium geometries and transition states', &
     'J. Comput. Chem.', '17', '49', '1996')
call add('lindh', 'Geometry optimization', 'Lindh, R. and Bernhardsson, A. and Karlström, G. and Malmqvist, P.-Å.', &
     'On the use of a Hessian model function in molecular geometry optimizations', 'Chem. Phys. Lett.', '241', '423', '1995')
call add('bofill', 'Transition state search', 'Bofill, J. M.', &
     'Updated Hessian matrix and the restricted step method for locating transition structures', &
     'J. Comput. Chem.', '15', '1', '1994')
call add('fukui', 'Reaction path (IRC)', 'Fukui, K.', 'The path of chemical reactions - the IRC approach', &
     'Acc. Chem. Res.', '14', '363', '1981')
call add('lqa', 'Reaction path (IRC)', 'Page, M. and McIver, J. W.', 'On evaluating the reaction path Hamiltonian', &
     'J. Chem. Phys.', '88', '922', '1988')
call add('hpc1', 'Reaction path (IRC)', 'Hratchian, H. P. and Schlegel, H. B.', &
     'Accurate reaction paths using a Hessian based predictor-corrector integrator', 'J. Chem. Phys.', '120', '9918', '2004')
call add('hpc2', 'Reaction path (IRC)', 'Hratchian, H. P. and Schlegel, H. B.', &
     'Using Hessian updating to increase the efficiency of a Hessian based predictor-corrector reaction path '// &
     'following method', 'J. Chem. Theory Comput.', '1', '61', '2005')
call add('mayer', 'Bond orders', 'Mayer, I.', 'Charge, bond order and valence in the AB initio SCF theory', &
     'Chem. Phys. Lett.', '97', '270', '1983')
call add('sad', 'Initial guess', 'van Lenthe, J. H. and Zwaans, R. and van Dam, H. J. J. and Guest, M. F.', &
     'Starting SCF calculations by superposition of atomic densities', 'J. Comput. Chem.', '27', '926', '2006')
call add('pyscf1', 'Effective core potential integrals (code from PySCF)', 'Sun, Q. and Berkelbach, T. C. and Blunt, N. S. '// &
     'and Booth, G. H. and Guo, S. and Li, Z. and Liu, J. and McClain, J. D. and Sayfutyarova, E. R. and Sharma, S. '// &
     'and Wouters, S. and Chan, G. K.-L.', &
     'PySCF: the Python-based simulations of chemistry framework', 'WIREs Comput. Mol. Sci.', '8', 'e1340', '2018')
call add('pyscf2', 'Effective core potential integrals (code from PySCF)', 'Sun, Q. and Zhang, X. and Banerjee, S. '// &
     'and Bao, P. and Barbry, M. and Blunt, N. S. and Bogdanov, N. A. and Booth, G. H. and Chen, J. and Cui, Z.-H. '// &
     'and others', 'Recent developments in the PySCF program package', 'J. Chem. Phys.', '153', '024109', '2020')
call add('hess', 'Analytic Hessian', 'Pople, J. A. and Krishnan, R. and Schlegel, H. B. and Binkley, J. S.', &
     'Derivative studies in Hartree-Fock and Møller-Plesset theories', 'Int. J. Quantum Chem. Symp.', '13', '225', '1979')
call add('mulliken', 'Population analysis', 'Mulliken, R. S.', &
     'Electronic population analysis on LCAO-MO molecular wave functions. I', 'J. Chem. Phys.', '23', '1833', '1955')
call add('lowdin', 'Population analysis', 'Löwdin, P.-O.', &
     'On the non-orthogonality problem connected with the use of atomic wave functions in the theory of molecules '// &
     'and crystals', 'J. Chem. Phys.', '18', '365', '1950')
call add('resp1', 'Population analysis', 'Bayly, C. I. and Cieplak, P. and Cornell, W. and Kollman, P. A.', &
     'A well-behaved electrostatic potential based method using charge restraints for deriving atomic charges: '// &
     'the RESP model', 'J. Phys. Chem.', '97', '10269', '1993')
call add('resp2', 'Population analysis', 'Singh, U. C. and Kollman, P. A.', &
     'An approach to computing electrostatic charges for molecules', 'J. Comput. Chem.', '5', '129', '1984')
call add('stab', 'Wavefunction stability', 'Seeger, R. and Pople, J. A.', &
     'Self-consistent molecular orbital methods. XVIII. Constraints and stability in Hartree-Fock theory', &
     'J. Chem. Phys.', '66', '3045', '1977')
end subroutine build_table

subroutine use_ref(key)
character(len=*), intent(in) :: key
integer :: i
do i = 1, nref
   if (trim(tab(i)%key) == key) then
      tab(i)%used = .true.
      return
   endif
enddo
end subroutine use_ref

function upper(s) result(r)
character(len=*), intent(in) :: s
character(len=len(s)) :: r
integer :: i, c
r = s
do i = 1, len(s)
   c = ichar(s(i:i))
   if (c >= ichar('a') .and. c <= ichar('z')) r(i:i) = char(c - 32)
enddo
end function upper

function lower(s) result(r)
character(len=*), intent(in) :: s
character(len=len(s)) :: r
integer :: i, c
r = s
do i = 1, len(s)
   c = ichar(s(i:i))
   if (c >= ichar('A') .and. c <= ichar('Z')) r(i:i) = char(c + 32)
enddo
end function lower

logical function uses_ecp(spec, bb)
type(engine_input_t), intent(in) :: spec
character(len=*), intent(in) :: bb
integer :: j
uses_ecp = .false.
if (len_trim(spec%ecplabel) > 0) uses_ecp = .true.
do j = 1, spec%ncenters
   if (len_trim(spec%atom_ecp(j)) > 0) uses_ecp = .true.
   if (index(bb, 'def2') > 0 .and. spec%atomchg(j) >= 37) uses_ecp = .true.
enddo
end function uses_ecp

subroutine select_refs(spec, is_hf)
type(engine_input_t), intent(in) :: spec
logical, intent(out) :: is_hf
character(len=64) :: fu, jj, kk, bb, aux
logical :: has_d3, has_d4, has_3c
integer :: i

fu = upper(adjustl(spec%functional))
jj = lower(adjustl(spec%J))
kk = lower(adjustl(spec%K))
bb = lower(adjustl(spec%baselabel))
aux = lower(adjustl(spec%ri_aux_basis))
is_hf = (trim(fu) == 'HF')

call use_ref('libcint')
call use_ref('diis')
call use_ref('adiis')
call use_ref('sad')
if (.not. (spec%opt_run .or. spec%irc_run)) then
   call use_ref('mulliken') ; call use_ref('lowdin')
endif
if (spec%resp_charges_on) then
   call use_ref('resp1') ; call use_ref('resp2')
endif
if (spec%stability_mode > 0) call use_ref('stab')
if (spec%freq .or. spec%irc_run .or. (spec%opt_run .and. spec%opt_ts)) call use_ref('hess')
if (uses_ecp(spec, bb)) then
   call use_ref('pyscf1') ; call use_ref('pyscf2')
endif
if (.not. is_hf) then
   call use_ref('becke88g') ; call use_ref('lebedev') ; call use_ref('treutler')
endif

if (trim(kk) == 'cosx') then
   call use_ref('cosx1') ; call use_ref('cosx2') ; call use_ref('cosx3')
endif
if (trim(jj) == 'ri') then
   call use_ref('rij1') ; call use_ref('rij2') ; call use_ref('rij3') ; call use_ref('rij4')
endif
if (trim(kk) == 'ri') call use_ref('rik')
if (len_trim(aux) == 0 .or. index(aux, 'def2') > 0) then
   if (index(bb, 'def2') > 0) then
      if (trim(jj) == 'ri' .or. trim(kk) == 'cosx') call use_ref('def2j')
      if (trim(kk) == 'ri') call use_ref('def2jk')
   endif
endif
if (index(bb, 'def2') > 0) call use_ref('def2')

has_3c = index(fu, '3C') > 0
has_d3 = index(fu, 'D3') > 0
has_d4 = index(fu, 'D4') > 0
if (has_d3) call use_ref('d3')
if (index(fu, 'D3BJ') > 0 .or. index(fu, 'BJ') > 0) then
   call use_ref('d3') ; call use_ref('d3bj')
endif
if (has_d4) then
   call use_ref('d4a') ; call use_ref('d4b')
endif
if (index(fu, '-D2') > 0 .or. index(fu, '_D2') > 0) call use_ref('d2')
if (index(fu, 'WB97X-D') > 0 .and. .not. has_d3 .and. .not. has_d4) call use_ref('wb97xd')
if (has_3c) then
   call use_ref('gcp')
   if (index(fu, 'B97-3C') > 0 .and. index(fu, 'WB97') == 0) then
      call use_ref('b973c') ; call use_ref('d3') ; call use_ref('d3bj')
   else if (index(fu, 'R2SCAN-3C') > 0) then
      call use_ref('r2scan3c') ; call use_ref('d4a') ; call use_ref('d4b')
   else if (index(fu, 'WB97X-3C') > 0) then
      call use_ref('wb97x3c') ; call use_ref('d4a') ; call use_ref('d4b')
   else if (index(fu, 'PBEH-3C') > 0) then
      call use_ref('pbeh3c') ; call use_ref('d3') ; call use_ref('d3bj')
   endif
endif
if (.not. is_hf) then
   if (xc_uses_vv10()) call use_ref('vv10')
endif
if (index(fu, 'WB97M-V') > 0) call use_ref('wb97mv')
if (index(fu, 'WB97X-V') > 0) call use_ref('wb97xv')
if (index(fu, 'B97M-V') > 0 .and. index(fu, 'WB97M-V') == 0) call use_ref('b97mv')

if (spec%cosmo_on) call use_ref('cosmo')
if (spec%cosmo_on .and. spec%cosmo_smd) call use_ref('smd')

if (spec%opt_run) then
   call use_ref('rfo') ; call use_ref('lindh')
   if (trim(spec%opt_coord) == 'ric' .or. trim(spec%opt_coord) == 'internal' .or. &
       trim(spec%opt_coord) == 'redundant') call use_ref('ric')
   if (spec%opt_ts) then
      call use_ref('bofill')
   endif
endif
if (spec%irc_run) then
   call use_ref('bofill') ; call use_ref('fukui') ; call use_ref('lqa')
   call use_ref('hpc1') ; call use_ref('hpc2')
   if (spec%irc_stop(2) > 0.0d0) call use_ref('mayer')
endif
do i = 1, 0
enddo
end subroutine select_refs

function plain(r) result(str)
type(ref_t), intent(in) :: r
character(len=900) :: str
character(len=300) :: au
integer :: p
au = r%authors
do
   p = index(au, ' and ')
   if (p == 0) exit
   au = au(1:p-1)//';'//au(p+4:)
enddo
if (au(len_trim(au):len_trim(au)) == '.') then
   str = trim(au)//' '//trim(r%title)//'. '//trim(r%journal)
else
   str = trim(au)//'. '//trim(r%title)//'. '//trim(r%journal)
endif
if (len_trim(r%volume) > 0) str = trim(str)//' '//trim(r%volume)
str = trim(str)//', '//trim(r%pages)//' ('//trim(r%year)//').'
end function plain

subroutine write_bib_entry(u, r)
integer, intent(in) :: u
type(ref_t), intent(in) :: r
write(u,'(A)') '@article{'//trim(r%key)//','
write(u,'(A)') '  author = {'//trim(r%authors)//'},'
write(u,'(A)') '  title = {{'//trim(r%title)//'}},'
write(u,'(A)') '  journaltitle = {'//trim(r%journal)//'},'
write(u,'(A)') '  date = {'//trim(r%year)//'},'
if (len_trim(r%volume) > 0) write(u,'(A)') '  volume = {'//trim(r%volume)//'},'
write(u,'(A)') '  pages = {'//trim(r%pages)//'},'
write(u,'(A)') '}'
write(u,'(A)') ''
end subroutine write_bib_entry

subroutine write_citations(spec, u, bibfile)
type(engine_input_t), intent(in) :: spec
integer, intent(in) :: u
character(len=*), intent(in), optional :: bibfile
logical :: is_hf, wrote_cat, anyb
integer :: i, j, ub, ios
character(len=60) :: lastcat
character(len=64) :: bb

call build_table()
do i = 1, nref
   tab(i)%used = .false.
enddo
call select_refs(spec, is_hf)

ub = 0
if (present(bibfile)) then
   if (len_trim(bibfile) > 0) then
      open(newunit=ub, file=trim(bibfile), status='replace', action='write', iostat=ios)
      if (ios /= 0) ub = 0
   endif
endif

write(u,'(A)') '[CITATIONS]'
write(u,'(A)') '  Direwolf (QuantaBricks), https://www.quantabricks.xyz'
write(u,'(A)') '  Please cite the methods used in this calculation:'
if (ub /= 0) then
   write(ub,'(A)') '% Literature of the methods used in this Direwolf calculation (BibLaTeX).'
   write(ub,'(A)') '@software{direwolf,'
   write(ub,'(A)') '  author = {Ma, Yizhou and Chen, Xin},'
   write(ub,'(A)') '  title = {Direwolf},'
   write(ub,'(A)') '  organization = {QuantaBricks},'
   write(ub,'(A)') '  url = {https://www.quantabricks.xyz},'
   write(ub,'(A)') '  date = {2026},'
   write(ub,'(A)') '}'
   write(ub,'(A)') ''
endif

lastcat = ''
do i = 1, nref
   if (.not. tab(i)%used) cycle
   if (trim(tab(i)%cat) /= trim(lastcat)) then
      write(u,'(A)') ''
      write(u,'(A)') '  '//trim(tab(i)%cat)
      lastcat = tab(i)%cat
   endif
   write(u,'(A)') '    '//trim(plain(tab(i)))
   if (ub /= 0) call write_bib_entry(ub, tab(i))
enddo

bb = lower(adjustl(spec%baselabel))
if (len_trim(bb) > 0 .and. index(bb, 'def2') == 0) then
   write(u,'(A)') ''
   write(u,'(A)') '  Basis sets'
   write(u,'(A)') '    '//trim(spec%baselabel)//': cite the original paper of the basis set; see its Basis Set Exchange entry'
endif
if (index(bb, 'def2') > 0) then
   anyb = .false.
   do j = 1, spec%ncenters
      if (spec%atomchg(j) >= 37) anyb = .true.
   enddo
   if (anyb) then
      write(u,'(A)') ''
      write(u,'(A)') '  Effective core potentials'
      write(u,'(A)') '    def2 ECPs (Z >= 37): cite the Basis Set Exchange entry of the def2 ECPs'
   endif
endif

if (.not. is_hf) then
   write(u,'(A)') ''
   write(u,'(A)') '  Exchange-correlation functionals (library and functional references from Libxc)'
   call xc_write_references(u, ub)
endif
write(u,'(A)') ''
write(u,'(A)') '  Libraries: libcint (above), OpenBLAS (https://www.openblas.net)'
if (index(upper(spec%functional), 'D3') > 0 .or. index(upper(spec%functional), '3C') > 0) &
   write(u,'(A)') '    simple-dftd3 (https://github.com/dftd3/simple-dftd3)'
if (index(upper(spec%functional), 'D4') > 0 .or. index(upper(spec%functional), 'R2SCAN-3C') > 0 .or. &
    index(upper(spec%functional), 'WB97X-3C') > 0) write(u,'(A)') '    dftd4 (https://github.com/dftd4/dftd4)'
if (index(upper(spec%functional), '3C') > 0) write(u,'(A)') '    gCP (https://github.com/grimme-lab/gcp)'
if (ub /= 0) then
   write(u,'(A)') '  BibLaTeX file: '//trim(bibfile)
   close(ub)
endif
write(u,'(A)') '[CITATIONSEND]'
write(u,'(A)') ''
end subroutine write_citations

end module mod_citations
