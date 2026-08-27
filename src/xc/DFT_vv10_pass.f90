! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! VV10 nonlocal-correlation Fock/energy passes - split out of DFT.f90


subroutine DFT_vv10_nlc_pass(Fxc_a, Fxc_b, Exc)
    use MOL_info
    use mod_vv10, only: vv10_evaluate, vv10_nlc_grid_build, vv10_nlc_basis_build, &
                         NLC_CHUNK, vv10_nlc_chunk_info, vv10_nlc_chunk_val
    use mod_profile, only: prof_start, prof_stop
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
    implicit none
    interface
       subroutine openblas_set_num_threads(num_threads) bind(c, name="openblas_set_num_threads")
       use iso_c_binding, only: c_int
       integer(c_int),value :: num_threads
       end subroutine openblas_set_num_threads
    end interface
    real(8),intent(inout) :: Fxc_a(nconts,nconts), Fxc_b(nconts,nconts), Exc
    real(8),allocatable :: coor(:,:), weight(:)
    real(8),allocatable :: rho(:), grad(:,:), sig(:), eps(:), Fn(:), Fg(:)
    real(8),allocatable :: Fxc_a_s(:,:)
    real(8) :: Enl, tmp
    integer :: np, i, mu, nu
    integer :: ic, cstart, cn, nsig

    call vv10_nlc_grid_build(np, coor, weight)
    allocate(rho(np), grad(3,np), sig(np))
    allocate(eps(np), Fn(np), Fg(np))

    call prof_start("vv10_gtoeval")
    call vv10_nlc_basis_build(np)
    call prof_stop("vv10_gtoeval")

    call openblas_set_num_threads(1)

    call prof_start("vv10_density_gemm")
    !$omp parallel private(ic,cstart,cn,nsig,i)
    block
      real(8),allocatable :: val0_chunk(:,:), val1_chunk(:,:,:), wtot_chunk(:,:)
      real(8),allocatable :: Pa_s(:,:), Pb_s(:,:)
      integer,allocatable :: sig_idx(:)
      real(4),pointer :: val0_ptr(:,:), val1_ptr(:,:,:)
      integer :: j
      allocate(val0_chunk(nconts,NLC_CHUNK), val1_chunk(nconts,3,NLC_CHUNK), wtot_chunk(nconts,NLC_CHUNK))
      allocate(Pa_s(nconts,nconts), Pb_s(nconts,nconts))
      !$omp do schedule(dynamic)
      do ic = 1,(np+NLC_CHUNK-1)/NLC_CHUNK
         call vv10_nlc_chunk_info(ic, cstart, cn, nsig, sig_idx)
         call vv10_nlc_chunk_val(ic, val0_ptr, val1_ptr)
         val0_chunk(1:nsig,1:cn) = real(val0_ptr(1:nsig,1:cn),8)
         val1_chunk(1:nsig,:,1:cn) = real(val1_ptr(1:nsig,:,1:cn),8)
         do j = 1,nsig
            Pa_s(1:nsig,j) = Pa(sig_idx(1:nsig),sig_idx(j))
         enddo
         if (Multi .ne. 1) then
            do j = 1,nsig
               Pb_s(1:nsig,j) = Pb(sig_idx(1:nsig),sig_idx(j))
            enddo
            call dgemm('N','N',nsig,cn,nsig,1.0d0,Pa_s,nconts,val0_chunk,nconts,0.0d0,wtot_chunk,nconts)
            call dgemm('N','N',nsig,cn,nsig,1.0d0,Pb_s,nconts,val0_chunk,nconts,1.0d0,wtot_chunk,nconts)
         else
            call dgemm('N','N',nsig,cn,nsig,2.0d0,Pa_s,nconts,val0_chunk,nconts,0.0d0,wtot_chunk,nconts)
         endif
         do i = 1,cn
            rho(cstart+i-1) = dot_product(val0_chunk(1:nsig,i), wtot_chunk(1:nsig,i))
            grad(1,cstart+i-1) = 2.0d0*dot_product(val1_chunk(1:nsig,1,i), wtot_chunk(1:nsig,i))
            grad(2,cstart+i-1) = 2.0d0*dot_product(val1_chunk(1:nsig,2,i), wtot_chunk(1:nsig,i))
            grad(3,cstart+i-1) = 2.0d0*dot_product(val1_chunk(1:nsig,3,i), wtot_chunk(1:nsig,i))
            sig(cstart+i-1) = grad(1,cstart+i-1)**2 + grad(2,cstart+i-1)**2 + grad(3,cstart+i-1)**2
         enddo
      enddo
      !$omp end do
      deallocate(val0_chunk, val1_chunk, wtot_chunk, Pa_s, Pb_s)
    end block
    !$omp end parallel
    call prof_stop("vv10_density_gemm")

    call prof_start("vv10_kernel_sum")
    call vv10_evaluate(np, rho, sig, weight, coor, Enl, eps, Fn, Fg)
    Exc = Exc + Enl
    call prof_stop("vv10_kernel_sum")

    call prof_start("vv10_fock_gemm")
    call openblas_set_num_threads(1)
    allocate(Fxc_a_s(nconts,nconts))
    Fxc_a_s = 0.0d0
    !$omp parallel private(ic,cstart,cn,nsig,i)
    block
      real(8),allocatable :: val0_chunk(:,:), val1_chunk(:,:,:)
      real(8),allocatable :: Gmat_chunk(:,:), Bmat_chunk(:,:), Cmat_chunk(:,:)
      real(8),allocatable :: Fxc_a_s_local(:,:), Fxc_compact(:,:)
      integer,allocatable :: sig_idx(:)
      real(4),pointer :: val0_ptr(:,:), val1_ptr(:,:,:)
      integer :: a, b
      allocate(val0_chunk(nconts,NLC_CHUNK), val1_chunk(nconts,3,NLC_CHUNK))
      allocate(Gmat_chunk(nconts,NLC_CHUNK), Bmat_chunk(nconts,NLC_CHUNK), Cmat_chunk(nconts,NLC_CHUNK))
      allocate(Fxc_a_s_local(nconts,nconts), Fxc_compact(nconts,nconts))
      Fxc_a_s_local = 0.0d0
      !$omp do schedule(static)
      do ic = 1,(np+NLC_CHUNK-1)/NLC_CHUNK
         call vv10_nlc_chunk_info(ic, cstart, cn, nsig, sig_idx)
         call vv10_nlc_chunk_val(ic, val0_ptr, val1_ptr)
         val0_chunk(1:nsig,1:cn) = real(val0_ptr(1:nsig,1:cn),8)
         val1_chunk(1:nsig,:,1:cn) = real(val1_ptr(1:nsig,:,1:cn),8)
         do i = 1,cn
            Gmat_chunk(1:nsig,i) = grad(1,cstart+i-1)*val1_chunk(1:nsig,1,i) + grad(2,cstart+i-1)*val1_chunk(1:nsig,2,i) &
                             + grad(3,cstart+i-1)*val1_chunk(1:nsig,3,i)
            Bmat_chunk(1:nsig,i) = weight(cstart+i-1)*Fn(cstart+i-1)*val0_chunk(1:nsig,i) &
                             + weight(cstart+i-1)*2.0d0*Fg(cstart+i-1)*Gmat_chunk(1:nsig,i)
            Cmat_chunk(1:nsig,i) = weight(cstart+i-1)*2.0d0*Fg(cstart+i-1)*val0_chunk(1:nsig,i)
         enddo
         Fxc_compact(1:nsig,1:nsig) = 0.0d0
         call dgemm('N','T',nsig,nsig,cn,1.0d0,Bmat_chunk,nconts,val0_chunk,nconts,1.0d0,Fxc_compact,nconts)
         call dgemm('N','T',nsig,nsig,cn,1.0d0,Cmat_chunk,nconts,Gmat_chunk,nconts,1.0d0,Fxc_compact,nconts)
         do b = 1,nsig
            do a = 1,nsig
               Fxc_a_s_local(sig_idx(a),sig_idx(b)) = Fxc_a_s_local(sig_idx(a),sig_idx(b)) + Fxc_compact(a,b)
            enddo
         enddo
      enddo
      !$omp end do
      block
        integer :: merge_tid, merge_nthreads
        merge_nthreads = omp_get_num_threads()
        do merge_tid = 0, merge_nthreads-1
           !$omp barrier
           if (omp_get_thread_num() .eq. merge_tid) then
              Fxc_a_s = Fxc_a_s + Fxc_a_s_local
           endif
        enddo
      end block
      deallocate(val0_chunk, val1_chunk, Gmat_chunk, Bmat_chunk, Cmat_chunk, Fxc_a_s_local, Fxc_compact)
    end block
    !$omp end parallel
    call prof_stop("vv10_fock_gemm")

    do mu = 1,nconts
       do nu = mu+1,nconts
          tmp = 0.5d0*(Fxc_a_s(mu,nu) + Fxc_a_s(nu,mu))
          Fxc_a_s(mu,nu) = tmp
          Fxc_a_s(nu,mu) = tmp
       enddo
    enddo

    Fxc_a = Fxc_a + Fxc_a_s
    Fxc_b = Fxc_b + Fxc_a_s
    deallocate(coor,weight,rho,grad,sig,eps,Fn,Fg)
    deallocate(Fxc_a_s)
end subroutine DFT_vv10_nlc_pass

subroutine DFT_vv10_postscf_aca(Enl_out)
    use MOL_info
    use mod_vv10, only: vv10_evaluate, vv10_nlc_grid_build, vv10_nlc_basis_build, NLC_CHUNK, &
                         vv10_nlc_chunk_info, vv10_nlc_chunk_val, vv10_omega_kappa, &
                         vv10_beta_val, vv10_report_set, vv10_rcut_val
    use mod_vv10_aca, only: vv10_evaluate_aca
    use mod_profile, only: prof_start, prof_stop
    implicit none
    interface
       subroutine openblas_set_num_threads(num_threads) bind(c, name="openblas_set_num_threads")
       use iso_c_binding, only: c_int
       integer(c_int),value :: num_threads
       end subroutine openblas_set_num_threads
    end interface
    real(8),intent(out) :: Enl_out
    real(8),allocatable :: coor(:,:), weight(:)
    real(8),allocatable :: rho(:), grad(:,:), sig(:)
    real(8),allocatable :: omega0(:), kappa(:), A_out(:)
    integer,allocatable :: act(:)
    integer :: np, n_active, i
    integer :: ic, cstart, cn, nsig
    integer,parameter :: ACA_CROSSOVER_N = 200000

    call vv10_nlc_grid_build(np, coor, weight)
    allocate(rho(np), grad(3,np), sig(np))

    call prof_start("vv10_gtoeval")
    call vv10_nlc_basis_build(np)
    call prof_stop("vv10_gtoeval")

    call openblas_set_num_threads(1)
    call prof_start("vv10_density_gemm")
    !$omp parallel private(ic,cstart,cn,nsig,i)
    block
      real(8),allocatable :: val0_chunk(:,:), val1_chunk(:,:,:), wtot_chunk(:,:)
      real(8),allocatable :: Pa_s(:,:), Pb_s(:,:)
      integer,allocatable :: sig_idx(:)
      real(4),pointer :: val0_ptr(:,:), val1_ptr(:,:,:)
      integer :: j
      allocate(val0_chunk(nconts,NLC_CHUNK), val1_chunk(nconts,3,NLC_CHUNK), wtot_chunk(nconts,NLC_CHUNK))
      allocate(Pa_s(nconts,nconts), Pb_s(nconts,nconts))
      !$omp do schedule(dynamic)
      do ic = 1,(np+NLC_CHUNK-1)/NLC_CHUNK
         call vv10_nlc_chunk_info(ic, cstart, cn, nsig, sig_idx)
         call vv10_nlc_chunk_val(ic, val0_ptr, val1_ptr)
         val0_chunk(1:nsig,1:cn) = real(val0_ptr(1:nsig,1:cn),8)
         val1_chunk(1:nsig,:,1:cn) = real(val1_ptr(1:nsig,:,1:cn),8)
         do j = 1,nsig
            Pa_s(1:nsig,j) = Pa(sig_idx(1:nsig),sig_idx(j))
         enddo
         if (Multi .ne. 1) then
            do j = 1,nsig
               Pb_s(1:nsig,j) = Pb(sig_idx(1:nsig),sig_idx(j))
            enddo
            call dgemm('N','N',nsig,cn,nsig,1.0d0,Pa_s,nconts,val0_chunk,nconts,0.0d0,wtot_chunk,nconts)
            call dgemm('N','N',nsig,cn,nsig,1.0d0,Pb_s,nconts,val0_chunk,nconts,1.0d0,wtot_chunk,nconts)
         else
            call dgemm('N','N',nsig,cn,nsig,2.0d0,Pa_s,nconts,val0_chunk,nconts,0.0d0,wtot_chunk,nconts)
         endif
         do i = 1,cn
            rho(cstart+i-1) = dot_product(val0_chunk(1:nsig,i), wtot_chunk(1:nsig,i))
            grad(1,cstart+i-1) = 2.0d0*dot_product(val1_chunk(1:nsig,1,i), wtot_chunk(1:nsig,i))
            grad(2,cstart+i-1) = 2.0d0*dot_product(val1_chunk(1:nsig,2,i), wtot_chunk(1:nsig,i))
            grad(3,cstart+i-1) = 2.0d0*dot_product(val1_chunk(1:nsig,3,i), wtot_chunk(1:nsig,i))
            sig(cstart+i-1) = grad(1,cstart+i-1)**2 + grad(2,cstart+i-1)**2 + grad(3,cstart+i-1)**2
         enddo
      enddo
      !$omp end do
      deallocate(val0_chunk, val1_chunk, wtot_chunk, Pa_s, Pb_s)
    end block
    !$omp end parallel
    call prof_stop("vv10_density_gemm")

    allocate(omega0(np), kappa(np), act(np), A_out(np))
    call vv10_omega_kappa(np, rho, sig, omega0, kappa, n_active, act)

    call prof_start("vv10_kernel_sum")
    if (n_active .ge. ACA_CROSSOVER_N) then
       call vv10_evaluate_aca(np, rho, sig, weight, coor, omega0, kappa, act, n_active, &
                               vv10_beta_val(), Enl_out, A_out, &
                               leaf_size_in=512, eta_in=1.0d0, tol_in=1.0d-5, &
                               rcut_in=vv10_rcut_val())
    else
       block
          real(8),allocatable :: eps(:), Fn(:), Fg(:)
          allocate(eps(np), Fn(np), Fg(np))
          call vv10_evaluate(np, rho, sig, weight, coor, Enl_out, eps, Fn, Fg)
          deallocate(eps, Fn, Fg)
       end block
    endif
    call prof_stop("vv10_kernel_sum")

    call vv10_report_set(Enl_out)

    deallocate(rho, grad, sig, omega0, kappa, act, A_out, coor, weight)
end subroutine DFT_vv10_postscf_aca
