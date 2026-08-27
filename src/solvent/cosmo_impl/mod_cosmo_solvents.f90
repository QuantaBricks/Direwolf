! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_cosmo_solvents: named-solvent lookup tables (plain epsilon, and

module mod_cosmo_solvents
implicit none

contains

function cosmo_solvent_epsilon(name) result(eps)
implicit none
character(len=*),intent(in) :: name
real(8) :: eps
character(len=32) :: key
integer :: i

key = trim(adjustl(name))
do i = 1,len_trim(key)
   if (key(i:i) >= 'A' .and. key(i:i) <= 'Z') key(i:i) = achar(iachar(key(i:i))+32)
enddo

select case (trim(key))
case ('water');               eps = 78.39d0
case ('methanol');             eps = 32.61d0
case ('ethanol');               eps = 24.85d0
case ('isopropanol','2-propanol','propan-2-ol'); eps = 19.26d0
case ('acetone');               eps = 20.49d0
case ('acetonitrile');          eps = 35.69d0
case ('dmso','dimethylsulfoxide'); eps = 46.83d0
case ('dmf','dimethylformamide'); eps = 36.71d0
case ('thf','tetrahydrofuran'); eps = 7.43d0
case ('dichloromethane','dcm','methylenechloride'); eps = 8.93d0
case ('chloroform');            eps = 4.71d0
case ('carbontetrachloride','ccl4'); eps = 2.23d0
case ('toluene');                eps = 2.37d0
case ('benzene');                eps = 2.27d0
case ('hexane','n-hexane');      eps = 1.88d0
case ('cyclohexane');            eps = 2.02d0
case ('dioxane','1,4-dioxane');  eps = 2.21d0
case ('diethylether','ether');   eps = 4.24d0
case ('ethylacetate');           eps = 5.99d0
case ('nitromethane');           eps = 35.87d0
case ('pyridine');               eps = 12.98d0
case ('formicacid');             eps = 51.10d0
case ('aceticacid');             eps = 6.25d0
case ('aniline');                eps = 6.89d0
case ('phenol');                 eps = 8.00d0
case ('nitrobenzene');           eps = 34.82d0
case ('ammonia');                eps = 22.40d0
case ('formamide');              eps = 108.94d0
case default
   print '("COSMO: unknown cosmo_solvent=''",A,"'' - not in the built-in solvent table.")', trim(name)
   print '("       Either fix the name or set cosmo_epsilon directly instead of cosmo_solvent.")'
   stop 1
end select
end function cosmo_solvent_epsilon

subroutine cosmo_solvent_smd_params(name, n, alpha, beta, gamma, phi, psi)
implicit none
character(len=*),intent(in) :: name
real(8),intent(out) :: n, alpha, beta, gamma, phi, psi
character(len=32) :: key
integer :: i

key = trim(adjustl(name))
do i = 1,len_trim(key)
   if (key(i:i) >= 'A' .and. key(i:i) <= 'Z') key(i:i) = achar(iachar(key(i:i))+32)
enddo

select case (trim(key))
case ('water');            n=1.3328d0; alpha=0.82d0; beta=0.35d0; gamma=-1.0d0;  phi=-1.0d0;  psi=-1.0d0
case ('methanol');         n=1.3288d0; alpha=0.43d0; beta=0.47d0; gamma=31.77d0; phi=0.0d0;   psi=0.0d0
case ('ethanol');          n=1.3611d0; alpha=0.37d0; beta=0.48d0; gamma=31.62d0; phi=0.0d0;   psi=0.0d0
case ('isopropanol','2-propanol','propan-2-ol'); &
                           n=1.3776d0; alpha=0.33d0; beta=0.56d0; gamma=30.13d0; phi=0.0d0;   psi=0.0d0
case ('acetone');          n=1.3588d0; alpha=0.04d0; beta=0.49d0; gamma=33.77d0; phi=0.0d0;   psi=0.0d0
case ('acetonitrile');     n=1.3442d0; alpha=0.07d0; beta=0.32d0; gamma=41.25d0; phi=0.0d0;   psi=0.0d0
case ('dmso','dimethylsulfoxide'); &
                           n=1.4783d0; alpha=0.0d0;  beta=0.88d0; gamma=61.78d0; phi=0.0d0;   psi=0.0d0
case ('dmf','dimethylformamide'); &
                           n=1.4305d0; alpha=0.0d0;  beta=0.74d0; gamma=49.56d0; phi=0.0d0;   psi=0.0d0
case ('thf','tetrahydrofuran'); &
                           n=1.4050d0; alpha=0.0d0;  beta=0.48d0; gamma=39.44d0; phi=0.0d0;   psi=0.0d0
case ('dichloromethane','dcm','methylenechloride'); &
                           n=1.4242d0; alpha=0.10d0; beta=0.05d0; gamma=39.15d0; phi=0.0d0;   psi=0.667d0
case ('chloroform');       n=1.4459d0; alpha=0.15d0; beta=0.02d0; gamma=38.39d0; phi=0.0d0;   psi=0.75d0
case ('carbontetrachloride','ccl4'); &
                           n=1.4601d0; alpha=0.0d0;  beta=0.0d0;  gamma=38.04d0; phi=0.0d0;   psi=0.8d0
case ('toluene');          n=1.4961d0; alpha=0.0d0;  beta=0.14d0; gamma=40.20d0; phi=0.857d0; psi=0.0d0
case ('benzene');          n=1.5011d0; alpha=0.0d0;  beta=0.14d0; gamma=40.62d0; phi=1.0d0;   psi=0.0d0
case ('hexane','n-hexane'); &
                           n=1.3749d0; alpha=0.0d0;  beta=0.0d0;  gamma=25.75d0; phi=0.0d0;   psi=0.0d0
case ('cyclohexane');      n=1.4266d0; alpha=0.0d0;  beta=0.0d0;  gamma=35.48d0; phi=0.0d0;   psi=0.0d0
case ('dioxane','1,4-dioxane'); &
                           n=1.4224d0; alpha=0.0d0;  beta=0.64d0; gamma=47.14d0; phi=0.0d0;   psi=0.0d0
case ('diethylether','ether'); &
                           n=1.3526d0; alpha=0.0d0;  beta=0.41d0; gamma=23.96d0; phi=0.0d0;   psi=0.0d0
case ('ethylacetate');     n=1.3723d0; alpha=0.0d0;  beta=0.45d0; gamma=33.67d0; phi=0.0d0;   psi=0.0d0
case ('nitromethane');     n=1.3817d0; alpha=0.06d0; beta=0.31d0; gamma=52.58d0; phi=0.0d0;   psi=0.0d0
case ('pyridine');         n=1.5095d0; alpha=0.0d0;  beta=0.52d0; gamma=52.62d0; phi=0.833d0; psi=0.0d0
case ('formicacid');       n=1.3714d0; alpha=0.75d0; beta=0.38d0; gamma=53.44d0; phi=0.0d0;   psi=0.0d0
case ('aceticacid');       n=1.3720d0; alpha=0.61d0; beta=0.44d0; gamma=39.01d0; phi=0.0d0;   psi=0.0d0
case ('aniline');          n=1.5863d0; alpha=0.26d0; beta=0.41d0; gamma=60.62d0; phi=0.857d0; psi=0.0d0
case ('nitrobenzene');     n=1.5562d0; alpha=0.0d0;  beta=0.28d0; gamma=57.54d0; phi=0.667d0; psi=0.0d0
case ('formamide');        n=1.4472d0; alpha=0.62d0; beta=0.60d0; gamma=82.08d0; phi=0.0d0;   psi=0.0d0
case ('ammonia','phenol')
   print '("COSMO: cosmo_solvent=''",A,"'' has an epsilon entry but no SMD descriptors ", &
          &"in the built-in table (not in pyscf''s SMD solvent_db either) - SMD is not ", &
          &"available for this solvent, only plain CPCM (cosmo_smd=.false.).")', trim(name)
   stop 1
case default
   print '("COSMO: unknown cosmo_solvent=''",A,"'' - not in the built-in SMD solvent table.")', trim(name)
   stop 1
end select
end subroutine cosmo_solvent_smd_params

end module mod_cosmo_solvents
