program exact
   use openacc
   use mpi
   implicit none
 ! 3D TDSE on 3 coupled PESs — multi-GPU version (MPI + OpenACC).
 ! Domain decomposed along igrid1 (one slab per GPU).
 ! Leapfrog (Störmer-Verlet) time integration.

 !--------------- grid size and domain --------------------------------
   integer, parameter :: ngrid1 = 1000, ngrid2 = 1000, ngrid3 = 1000
   real*8, parameter :: x1min = -4.0d0, x1max =  8.0d0
   real*8, parameter :: x2min = -4.0d0, x2max =  8.0d0
   real*8, parameter :: x3min = -4.0d0, x3max =  8.0d0

 !--------------- model parameters ------------------------------------
   real*8, parameter :: w1    = 0.25d0
   real*8, parameter :: w2    = 0.025d0
   real*8, parameter :: delta = 0.01d0
   real*8, parameter :: epsil = 0.00d0
   real*8, parameter :: c     = 0.025d0
   real*8, parameter :: pmass = 1.845d3

 !--------------- integration parameters ------------------------------
   real*8, parameter :: deltat   = 0.01d0
   real*8, parameter :: tstepmax = 3.1d4

 !--------------- initial conditions ----------------------------------
   real*8, parameter :: R1bar = -1.0d0, R2bar = 0.0d0, R3bar = 0.0d0
   real*8, parameter :: P1bar = 10.0d0, P2bar = 10.0d0, P3bar = 10.0d0

 !--------------- other parameters ------------------------------------
   real*8, parameter :: pi     = 3.14159265359d0
   real*8, parameter :: alpha1 = 6.0d0, alpha2 = 6.0d0, alpha3 = 6.0d0

 !--------------- wavefunction arrays (decomposed along igrid1) -------
 ! Local index -1,0 = left halos, 1..local_n1 = interior, local_n1+1,+2 = right halos.
 ! Two-cell halos support 4th-order finite-difference stencils in compute_momentum
 ! (via exchange_halo_wide). The propagator uses only the inner halos 0 and
 ! local_n1+1 and calls the narrow exchange_halo which transfers a single slab.
   real*8, allocatable :: wfr(:,:,:,:)
   real*8, allocatable :: wfi(:,:,:,:)

 !--------------- halo exchange buffers --------------------------------
 ! Last dim is the slab index: 1 = innermost (adjacent to interior),
 !                             2 = outermost (one cell deeper into neighbor).
   real*8, allocatable :: sbuf_l(:,:,:,:), sbuf_r(:,:,:,:)
   real*8, allocatable :: rbuf_l(:,:,:,:), rbuf_r(:,:,:,:)

 !--------------- workspace variables ---------------------------------
   real*8 :: deltaR1, deltaR2, deltaR3
   real*8 :: prefac1, prefac2, prefac3
   real*8 :: gauss1, gauss2, gauss3
   real*8 :: R1, R2, R3, R1disp, R2disp, R3disp
   real*8 :: PbarRdisp, cosPbarRdisp, sinPbarRdisp
   real*8 :: inv_ke1a, inv_ke2a, inv_ke3a    ! 4th-order: coeff. on psi(+/-1)
   real*8 :: inv_ke1b, inv_ke2b, inv_ke3b    ! 4th-order: coeff. on psi(+/-2)
   real*8 :: diag_ke, halfdeltat
   real*8 :: norm, norm_local, V1, V2, V3, cR2, dval, dens
   real*8 :: pop(3), pop_local(3)
   real*8 :: pop1_l, pop2_l, pop3_l
   real*8 :: start, finish
   integer :: igrid1, igrid2, igrid3, istate, iloc1
   integer*8 :: n

 !--------------- momentum computation variables -----------------------
   real*8 :: inv_2dr1, inv_2dr2, inv_2dr3

 !--------------- MPI / decomposition variables -----------------------
   integer :: myrank, nprocs, ierr
   integer :: left_rank, right_rank
   integer :: local_n1, global_start
   integer :: ninterior, base_n, remainder_n
   integer :: halo_size
   integer :: status_mpi(MPI_STATUS_SIZE)

 !--------------- for initial wf output --------------------------------
   real*8, allocatable :: rho1d_local(:), rho1d_global(:)

 ! =====================================================================
 ! MPI + GPU SETUP
 ! =====================================================================
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

   call acc_set_device_num(myrank, acc_device_nvidia)

 ! =====================================================================
 ! PRECOMPUTE CONSTANTS
 ! =====================================================================
   deltaR1    = (x1max - x1min) / dble(ngrid1)
   deltaR2    = (x2max - x2min) / dble(ngrid2)
   deltaR3    = (x3max - x3min) / dble(ngrid3)
   prefac1    = sqrt(sqrt(2.0d0 * alpha1 / pi))
   prefac2    = sqrt(sqrt(2.0d0 * alpha2 / pi))
   prefac3    = sqrt(sqrt(2.0d0 * alpha3 / pi))
   ! 4th-order central-difference Laplacian for T = -nabla^2/(2m):
   !   (nabla^2 psi)/(2m) = sum_n { 16 (psi(+/-1)) - (psi(+/-2)) - 30 psi } / (24 m dR_n^2)
   inv_ke1a   =  16.0d0 / (24.0d0 * deltaR1**2 * pmass)
   inv_ke2a   =  16.0d0 / (24.0d0 * deltaR2**2 * pmass)
   inv_ke3a   =  16.0d0 / (24.0d0 * deltaR3**2 * pmass)
   inv_ke1b   =  -1.0d0 / (24.0d0 * deltaR1**2 * pmass)
   inv_ke2b   =  -1.0d0 / (24.0d0 * deltaR2**2 * pmass)
   inv_ke3b   =  -1.0d0 / (24.0d0 * deltaR3**2 * pmass)
   diag_ke    = -30.0d0 / (24.0d0 * pmass) * &
                (1.0d0/deltaR1**2 + 1.0d0/deltaR2**2 + 1.0d0/deltaR3**2)
   halfdeltat = 0.5d0 * deltat
   inv_2dr1   = 0.5d0 / deltaR1
   inv_2dr2   = 0.5d0 / deltaR2
   inv_2dr3   = 0.5d0 / deltaR3

 ! =====================================================================
 ! DOMAIN DECOMPOSITION ALONG igrid1
 ! =====================================================================
   ninterior   = ngrid1 - 2          ! interior points 2..ngrid1-1
   base_n      = ninterior / nprocs
   remainder_n = mod(ninterior, nprocs)
   if (myrank < remainder_n) then
     local_n1     = base_n + 1
     global_start = 2 + myrank * (base_n + 1)
   else
     local_n1     = base_n
     global_start = 2 + remainder_n * (base_n + 1) + (myrank - remainder_n) * base_n
   end if
   ! This rank owns global igrid1 = global_start .. global_start + local_n1 - 1

   left_rank  = myrank - 1
   right_rank = myrank + 1
   if (left_rank  < 0)      left_rank  = MPI_PROC_NULL
   if (right_rank >= nprocs) right_rank = MPI_PROC_NULL

   halo_size = ngrid3 * ngrid2 * 3

   if (myrank == 0) then
     print*, "Running on", nprocs, "GPUs"
     print*, "Grid:", ngrid1, "x", ngrid2, "x", ngrid3
     print*, "Slabs per GPU:", local_n1
   end if

 ! =====================================================================
 ! ALLOCATE ARRAYS
 ! =====================================================================
   allocate(wfr(ngrid3, ngrid2, -1:local_n1+2, 3))
   allocate(wfi(ngrid3, ngrid2, -1:local_n1+2, 3))
   allocate(sbuf_l(ngrid3, ngrid2, 3, 2))
   allocate(sbuf_r(ngrid3, ngrid2, 3, 2))
   allocate(rbuf_l(ngrid3, ngrid2, 3, 2))
   allocate(rbuf_r(ngrid3, ngrid2, 3, 2))

   wfr    = 0.0d0;  wfi    = 0.0d0
   sbuf_l = 0.0d0;  sbuf_r = 0.0d0
   rbuf_l = 0.0d0;  rbuf_r = 0.0d0

 ! =====================================================================
 ! INITIALIZE WAVEFUNCTION (on host)
 ! =====================================================================
   do iloc1 = 1, local_n1
     igrid1 = global_start + iloc1 - 1
     R1 = x1min + (dble(igrid1)-0.5d0)*deltaR1
     do igrid2 = 2, ngrid2-1
       R2 = x2min + (dble(igrid2)-0.5d0)*deltaR2
       do igrid3 = 2, ngrid3-1
         R3 = x3min + (dble(igrid3)-0.5d0)*deltaR3
         R1disp       = R1 - R1bar
         R2disp       = R2 - R2bar
         R3disp       = R3 - R3bar
         gauss1       = exp(-alpha1 * R1disp**2)
         gauss2       = exp(-alpha2 * R2disp**2)
         gauss3       = exp(-alpha3 * R3disp**2)
         PbarRdisp    = P1bar*R1disp + P2bar*R2disp + P3bar*R3disp
         cosPbarRdisp = cos(PbarRdisp)
         sinPbarRdisp = sin(PbarRdisp)
         wfr(igrid3,igrid2,iloc1,1) = prefac1*prefac2*prefac3*gauss1*gauss2*gauss3*cosPbarRdisp
         wfi(igrid3,igrid2,iloc1,1) = prefac1*prefac2*prefac3*gauss1*gauss2*gauss3*sinPbarRdisp
       end do
     end do
   end do

 ! =====================================================================
 ! WRITE INITIAL WAVEFUNCTION (1D marginal along R1)
 ! =====================================================================
   allocate(rho1d_local(ngrid1), rho1d_global(ngrid1))
   rho1d_local = 0.0d0
   do iloc1 = 1, local_n1
     igrid1 = global_start + iloc1 - 1
     do igrid2 = 1, ngrid2
       do igrid3 = 1, ngrid3
         rho1d_local(igrid1) = rho1d_local(igrid1) + &
           wfr(igrid3,igrid2,iloc1,1)**2 + wfi(igrid3,igrid2,iloc1,1)**2
       end do
     end do
   end do
   call MPI_Reduce(rho1d_local, rho1d_global, ngrid1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   if (myrank == 0) then
     open(unit=21, file="wf_init_R1.dat", status="replace")
     do igrid1 = 1, ngrid1
       R1 = x1min + (dble(igrid1)-0.5d0)*deltaR1
       write(21,*) R1, rho1d_global(igrid1) * deltaR2 * deltaR3
     end do
     close(21)
   end if
   deallocate(rho1d_local, rho1d_global)

 ! =====================================================================
 ! INITIAL NORM & POPULATIONS
 ! =====================================================================
   norm_local = 0.0d0;  pop_local = 0.0d0
   do iloc1 = 1, local_n1
     do igrid2 = 2, ngrid2-1
       do igrid3 = 2, ngrid3-1
         do istate = 1, 3
           dens = wfr(igrid3,igrid2,iloc1,istate)**2 + wfi(igrid3,igrid2,iloc1,istate)**2
           norm_local       = norm_local       + dens
           pop_local(istate) = pop_local(istate) + dens
         end do
       end do
     end do
   end do
   call MPI_Allreduce(norm_local, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(pop_local,  pop,  3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   norm = norm * deltaR1 * deltaR2 * deltaR3
   pop  = pop  * deltaR1 * deltaR2 * deltaR3

   if (myrank == 0) then
     open(unit=30, file="norm.dat", status="replace")
     open(unit=31, file="pop.dat",  status="replace")
     open(unit=32, file="mom.dat",  status="replace")
     write(30,*) 0_8, norm
     write(31,*) 0_8, pop
   end if

 ! =====================================================================
 ! COPY DATA TO GPUs
 ! =====================================================================
   start = MPI_Wtime()
   !$acc enter data copyin(wfr, wfi, sbuf_l, sbuf_r, rbuf_l, rbuf_r)

   call compute_momentum(0_8)

   n = 1_8

 ! =====================================================================
 ! MAIN PROPAGATION LOOP
 ! =====================================================================
   do while (n < tstepmax)

     ! ----- Halo exchange for wfi (2-cell, for 4th-order stencil) -----
     call exchange_halo_wide(wfi)

     ! ----- Step 1: dwfrdt from wfi, advance wfr by dt/2 -----
     ! 4th-order Laplacian uses a 5-point stencil per direction; R1 reaches
     ! into the 2-cell halo, R2/R3 loops shrink to (3..ngrid-2) so the
     ! stencil stays in bounds. Cells at ig{2,3} = 2 and ngrid{2,3}-1 stay
     ! at their initial value (zero for a localized wavefunction) -- a
     ! 2-cell-thick Dirichlet wall instead of 1-cell.
     !$acc parallel loop collapse(3) present(wfr, wfi)
     do iloc1 = 1, local_n1
       do igrid2 = 3, ngrid2-2
         do igrid3 = 3, ngrid3-2
           igrid1 = global_start + iloc1 - 1
           R1  = x1min + (dble(igrid1)-0.5d0)*deltaR1
           R2  = x2min + (dble(igrid2)-0.5d0)*deltaR2
           V1  = -w1 * R1
           V2  =  w2 * R1 - delta
           V3  =  w2 * R1 - 2.0d0 * delta
           cR2 = c * R2

           dval = -( &
               (wfi(igrid3,igrid2,iloc1-1,1) + wfi(igrid3,igrid2,iloc1+1,1)) * inv_ke1a + &
               (wfi(igrid3,igrid2,iloc1-2,1) + wfi(igrid3,igrid2,iloc1+2,1)) * inv_ke1b + &
               (wfi(igrid3,igrid2-1,iloc1,1) + wfi(igrid3,igrid2+1,iloc1,1)) * inv_ke2a + &
               (wfi(igrid3,igrid2-2,iloc1,1) + wfi(igrid3,igrid2+2,iloc1,1)) * inv_ke2b + &
               (wfi(igrid3-1,igrid2,iloc1,1) + wfi(igrid3+1,igrid2,iloc1,1)) * inv_ke3a + &
               (wfi(igrid3-2,igrid2,iloc1,1) + wfi(igrid3+2,igrid2,iloc1,1)) * inv_ke3b + &
                wfi(igrid3,igrid2,iloc1,1) * diag_ke) + &
               V1 * wfi(igrid3,igrid2,iloc1,1) + &
               cR2 * (wfi(igrid3,igrid2,iloc1,2) + wfi(igrid3,igrid2,iloc1,3))
           wfr(igrid3,igrid2,iloc1,1) = wfr(igrid3,igrid2,iloc1,1) + halfdeltat * dval

           dval = -( &
               (wfi(igrid3,igrid2,iloc1-1,2) + wfi(igrid3,igrid2,iloc1+1,2)) * inv_ke1a + &
               (wfi(igrid3,igrid2,iloc1-2,2) + wfi(igrid3,igrid2,iloc1+2,2)) * inv_ke1b + &
               (wfi(igrid3,igrid2-1,iloc1,2) + wfi(igrid3,igrid2+1,iloc1,2)) * inv_ke2a + &
               (wfi(igrid3,igrid2-2,iloc1,2) + wfi(igrid3,igrid2+2,iloc1,2)) * inv_ke2b + &
               (wfi(igrid3-1,igrid2,iloc1,2) + wfi(igrid3+1,igrid2,iloc1,2)) * inv_ke3a + &
               (wfi(igrid3-2,igrid2,iloc1,2) + wfi(igrid3+2,igrid2,iloc1,2)) * inv_ke3b + &
                wfi(igrid3,igrid2,iloc1,2) * diag_ke) + &
               V2 * wfi(igrid3,igrid2,iloc1,2) + &
               cR2 * wfi(igrid3,igrid2,iloc1,1)
           wfr(igrid3,igrid2,iloc1,2) = wfr(igrid3,igrid2,iloc1,2) + halfdeltat * dval

           dval = -( &
               (wfi(igrid3,igrid2,iloc1-1,3) + wfi(igrid3,igrid2,iloc1+1,3)) * inv_ke1a + &
               (wfi(igrid3,igrid2,iloc1-2,3) + wfi(igrid3,igrid2,iloc1+2,3)) * inv_ke1b + &
               (wfi(igrid3,igrid2-1,iloc1,3) + wfi(igrid3,igrid2+1,iloc1,3)) * inv_ke2a + &
               (wfi(igrid3,igrid2-2,iloc1,3) + wfi(igrid3,igrid2+2,iloc1,3)) * inv_ke2b + &
               (wfi(igrid3-1,igrid2,iloc1,3) + wfi(igrid3+1,igrid2,iloc1,3)) * inv_ke3a + &
               (wfi(igrid3-2,igrid2,iloc1,3) + wfi(igrid3+2,igrid2,iloc1,3)) * inv_ke3b + &
                wfi(igrid3,igrid2,iloc1,3) * diag_ke) + &
               V3 * wfi(igrid3,igrid2,iloc1,3) + &
               cR2 * wfi(igrid3,igrid2,iloc1,1)
           wfr(igrid3,igrid2,iloc1,3) = wfr(igrid3,igrid2,iloc1,3) + halfdeltat * dval
         end do
       end do
     end do

     ! ----- Halo exchange for wfr (2-cell, for 4th-order stencil) -----
     call exchange_halo_wide(wfr)

     ! ----- Step 2: dwfidt from wfr, advance wfi by dt -----
     !$acc parallel loop collapse(3) present(wfr, wfi)
     do iloc1 = 1, local_n1
       do igrid2 = 3, ngrid2-2
         do igrid3 = 3, ngrid3-2
           igrid1 = global_start + iloc1 - 1
           R1  = x1min + (dble(igrid1)-0.5d0)*deltaR1
           R2  = x2min + (dble(igrid2)-0.5d0)*deltaR2
           V1  = -w1 * R1
           V2  =  w2 * R1 - delta
           V3  =  w2 * R1 - 2.0d0 * delta
           cR2 = c * R2

           dval = ( &
               (wfr(igrid3,igrid2,iloc1-1,1) + wfr(igrid3,igrid2,iloc1+1,1)) * inv_ke1a + &
               (wfr(igrid3,igrid2,iloc1-2,1) + wfr(igrid3,igrid2,iloc1+2,1)) * inv_ke1b + &
               (wfr(igrid3,igrid2-1,iloc1,1) + wfr(igrid3,igrid2+1,iloc1,1)) * inv_ke2a + &
               (wfr(igrid3,igrid2-2,iloc1,1) + wfr(igrid3,igrid2+2,iloc1,1)) * inv_ke2b + &
               (wfr(igrid3-1,igrid2,iloc1,1) + wfr(igrid3+1,igrid2,iloc1,1)) * inv_ke3a + &
               (wfr(igrid3-2,igrid2,iloc1,1) + wfr(igrid3+2,igrid2,iloc1,1)) * inv_ke3b + &
                wfr(igrid3,igrid2,iloc1,1) * diag_ke) - &
               V1 * wfr(igrid3,igrid2,iloc1,1) - &
               cR2 * (wfr(igrid3,igrid2,iloc1,2) + wfr(igrid3,igrid2,iloc1,3))
           wfi(igrid3,igrid2,iloc1,1) = wfi(igrid3,igrid2,iloc1,1) + deltat * dval

           dval = ( &
               (wfr(igrid3,igrid2,iloc1-1,2) + wfr(igrid3,igrid2,iloc1+1,2)) * inv_ke1a + &
               (wfr(igrid3,igrid2,iloc1-2,2) + wfr(igrid3,igrid2,iloc1+2,2)) * inv_ke1b + &
               (wfr(igrid3,igrid2-1,iloc1,2) + wfr(igrid3,igrid2+1,iloc1,2)) * inv_ke2a + &
               (wfr(igrid3,igrid2-2,iloc1,2) + wfr(igrid3,igrid2+2,iloc1,2)) * inv_ke2b + &
               (wfr(igrid3-1,igrid2,iloc1,2) + wfr(igrid3+1,igrid2,iloc1,2)) * inv_ke3a + &
               (wfr(igrid3-2,igrid2,iloc1,2) + wfr(igrid3+2,igrid2,iloc1,2)) * inv_ke3b + &
                wfr(igrid3,igrid2,iloc1,2) * diag_ke) - &
               V2 * wfr(igrid3,igrid2,iloc1,2) - &
               cR2 * wfr(igrid3,igrid2,iloc1,1)
           wfi(igrid3,igrid2,iloc1,2) = wfi(igrid3,igrid2,iloc1,2) + deltat * dval

           dval = ( &
               (wfr(igrid3,igrid2,iloc1-1,3) + wfr(igrid3,igrid2,iloc1+1,3)) * inv_ke1a + &
               (wfr(igrid3,igrid2,iloc1-2,3) + wfr(igrid3,igrid2,iloc1+2,3)) * inv_ke1b + &
               (wfr(igrid3,igrid2-1,iloc1,3) + wfr(igrid3,igrid2+1,iloc1,3)) * inv_ke2a + &
               (wfr(igrid3,igrid2-2,iloc1,3) + wfr(igrid3,igrid2+2,iloc1,3)) * inv_ke2b + &
               (wfr(igrid3-1,igrid2,iloc1,3) + wfr(igrid3+1,igrid2,iloc1,3)) * inv_ke3a + &
               (wfr(igrid3-2,igrid2,iloc1,3) + wfr(igrid3+2,igrid2,iloc1,3)) * inv_ke3b + &
                wfr(igrid3,igrid2,iloc1,3) * diag_ke) - &
               V3 * wfr(igrid3,igrid2,iloc1,3) - &
               cR2 * wfr(igrid3,igrid2,iloc1,1)
           wfi(igrid3,igrid2,iloc1,3) = wfi(igrid3,igrid2,iloc1,3) + deltat * dval
         end do
       end do
     end do

     ! ----- Halo exchange for wfi (updated by Step 2; 2-cell) -----
     call exchange_halo_wide(wfi)

     ! ----- Step 3: dwfrdt from wfi, advance wfr by dt/2 -----
     !$acc parallel loop collapse(3) present(wfr, wfi)
     do iloc1 = 1, local_n1
       do igrid2 = 3, ngrid2-2
         do igrid3 = 3, ngrid3-2
           igrid1 = global_start + iloc1 - 1
           R1  = x1min + (dble(igrid1)-0.5d0)*deltaR1
           R2  = x2min + (dble(igrid2)-0.5d0)*deltaR2
           V1  = -w1 * R1
           V2  =  w2 * R1 - delta
           V3  =  w2 * R1 - 2.0d0 * delta
           cR2 = c * R2

           dval = -( &
               (wfi(igrid3,igrid2,iloc1-1,1) + wfi(igrid3,igrid2,iloc1+1,1)) * inv_ke1a + &
               (wfi(igrid3,igrid2,iloc1-2,1) + wfi(igrid3,igrid2,iloc1+2,1)) * inv_ke1b + &
               (wfi(igrid3,igrid2-1,iloc1,1) + wfi(igrid3,igrid2+1,iloc1,1)) * inv_ke2a + &
               (wfi(igrid3,igrid2-2,iloc1,1) + wfi(igrid3,igrid2+2,iloc1,1)) * inv_ke2b + &
               (wfi(igrid3-1,igrid2,iloc1,1) + wfi(igrid3+1,igrid2,iloc1,1)) * inv_ke3a + &
               (wfi(igrid3-2,igrid2,iloc1,1) + wfi(igrid3+2,igrid2,iloc1,1)) * inv_ke3b + &
                wfi(igrid3,igrid2,iloc1,1) * diag_ke) + &
               V1 * wfi(igrid3,igrid2,iloc1,1) + &
               cR2 * (wfi(igrid3,igrid2,iloc1,2) + wfi(igrid3,igrid2,iloc1,3))
           wfr(igrid3,igrid2,iloc1,1) = wfr(igrid3,igrid2,iloc1,1) + halfdeltat * dval

           dval = -( &
               (wfi(igrid3,igrid2,iloc1-1,2) + wfi(igrid3,igrid2,iloc1+1,2)) * inv_ke1a + &
               (wfi(igrid3,igrid2,iloc1-2,2) + wfi(igrid3,igrid2,iloc1+2,2)) * inv_ke1b + &
               (wfi(igrid3,igrid2-1,iloc1,2) + wfi(igrid3,igrid2+1,iloc1,2)) * inv_ke2a + &
               (wfi(igrid3,igrid2-2,iloc1,2) + wfi(igrid3,igrid2+2,iloc1,2)) * inv_ke2b + &
               (wfi(igrid3-1,igrid2,iloc1,2) + wfi(igrid3+1,igrid2,iloc1,2)) * inv_ke3a + &
               (wfi(igrid3-2,igrid2,iloc1,2) + wfi(igrid3+2,igrid2,iloc1,2)) * inv_ke3b + &
                wfi(igrid3,igrid2,iloc1,2) * diag_ke) + &
               V2 * wfi(igrid3,igrid2,iloc1,2) + &
               cR2 * wfi(igrid3,igrid2,iloc1,1)
           wfr(igrid3,igrid2,iloc1,2) = wfr(igrid3,igrid2,iloc1,2) + halfdeltat * dval

           dval = -( &
               (wfi(igrid3,igrid2,iloc1-1,3) + wfi(igrid3,igrid2,iloc1+1,3)) * inv_ke1a + &
               (wfi(igrid3,igrid2,iloc1-2,3) + wfi(igrid3,igrid2,iloc1+2,3)) * inv_ke1b + &
               (wfi(igrid3,igrid2-1,iloc1,3) + wfi(igrid3,igrid2+1,iloc1,3)) * inv_ke2a + &
               (wfi(igrid3,igrid2-2,iloc1,3) + wfi(igrid3,igrid2+2,iloc1,3)) * inv_ke2b + &
               (wfi(igrid3-1,igrid2,iloc1,3) + wfi(igrid3+1,igrid2,iloc1,3)) * inv_ke3a + &
               (wfi(igrid3-2,igrid2,iloc1,3) + wfi(igrid3+2,igrid2,iloc1,3)) * inv_ke3b + &
                wfi(igrid3,igrid2,iloc1,3) * diag_ke) + &
               V3 * wfi(igrid3,igrid2,iloc1,3) + &
               cR2 * wfi(igrid3,igrid2,iloc1,1)
           wfr(igrid3,igrid2,iloc1,3) = wfr(igrid3,igrid2,iloc1,3) + halfdeltat * dval
         end do
       end do
     end do

     ! ----- Diagnostics every 100 steps -----
     if (mod(n,100)==0) then
       norm_local = 0.0d0
       pop1_l = 0.0d0;  pop2_l = 0.0d0;  pop3_l = 0.0d0
       !$acc parallel loop collapse(3) present(wfr, wfi) &
       !$acc reduction(+:norm_local, pop1_l, pop2_l, pop3_l)
       do iloc1 = 1, local_n1
         do igrid2 = 2, ngrid2-1
           do igrid3 = 2, ngrid3-1
             pop1_l = pop1_l + wfr(igrid3,igrid2,iloc1,1)**2 + wfi(igrid3,igrid2,iloc1,1)**2
             pop2_l = pop2_l + wfr(igrid3,igrid2,iloc1,2)**2 + wfi(igrid3,igrid2,iloc1,2)**2
             pop3_l = pop3_l + wfr(igrid3,igrid2,iloc1,3)**2 + wfi(igrid3,igrid2,iloc1,3)**2
           end do
         end do
       end do
       norm_local = pop1_l + pop2_l + pop3_l
       pop_local(1) = pop1_l;  pop_local(2) = pop2_l;  pop_local(3) = pop3_l

       call MPI_Allreduce(norm_local, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call MPI_Allreduce(pop_local,  pop,  3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       norm = norm * deltaR1 * deltaR2 * deltaR3
       pop  = pop  * deltaR1 * deltaR2 * deltaR3
       if (myrank == 0) then
         write(30,*) n, norm
         write(31,*) n, pop
       end if

       call compute_momentum(n)

       if (mod(n,1000)==0) then
         !$acc update self(wfr, wfi)
         call dump_density(n)
       end if
     end if

     n = n + 1_8
   end do

 ! =====================================================================
 ! CLEANUP
 ! =====================================================================
   !$acc exit data delete(wfr, wfi, sbuf_l, sbuf_r, rbuf_l, rbuf_r)

   if (myrank == 0) then
     close(30)
     close(31)
     close(32)
   end if

   finish = MPI_Wtime()
   if (myrank == 0) print*, "Program completes in ", finish-start, " seconds"

   deallocate(wfr, wfi, sbuf_l, sbuf_r, rbuf_l, rbuf_r)
   call MPI_Finalize(ierr)

 contains

 ! =====================================================================
 ! NARROW HALO EXCHANGE (single slab per side). Pack on GPU, MPI on host,
 ! unpack on GPU. Used by the propagator; only fills the inner halos at
 ! il1 = 0 and il1 = local_n1+1. The outer halos (-1, local_n1+2) are
 ! refreshed separately by exchange_halo_wide before compute_momentum.
 ! Slot 1 of the (3,2) buffers holds the innermost slab; slot 2 is unused
 ! here and left untouched.
 ! =====================================================================
   subroutine exchange_halo(arr)
     implicit none
     real*8, intent(inout) :: arr(ngrid3, ngrid2, -1:local_n1+2, 3)
     integer :: i2, i3, is, ierr_h
     integer :: stat_h(MPI_STATUS_SIZE)

     ! Pack innermost boundary slices into slot 1 of send buffers (on GPU)
     !$acc parallel loop collapse(3) present(arr, sbuf_l, sbuf_r)
     do is = 1, 3
       do i2 = 1, ngrid2
         do i3 = 1, ngrid3
           sbuf_l(i3, i2, is, 1) = arr(i3, i2, 1,        is)
           sbuf_r(i3, i2, is, 1) = arr(i3, i2, local_n1, is)
         end do
       end do
     end do

     ! GPU -> host (slot 1 only; one slab's worth = halo_size elements)
     !$acc update self(sbuf_l(:,:,:,1:1), sbuf_r(:,:,:,1:1))

     ! MPI exchange (MPI_PROC_NULL = no-op at domain boundaries -> halo stays 0).
     ! Count = halo_size sends exactly slot 1 thanks to Fortran column-major layout.
     call MPI_Sendrecv(sbuf_r, halo_size, MPI_DOUBLE_PRECISION, right_rank, 0, &
                       rbuf_l, halo_size, MPI_DOUBLE_PRECISION, left_rank,  0, &
                       MPI_COMM_WORLD, stat_h, ierr_h)
     call MPI_Sendrecv(sbuf_l, halo_size, MPI_DOUBLE_PRECISION, left_rank,  1, &
                       rbuf_r, halo_size, MPI_DOUBLE_PRECISION, right_rank, 1, &
                       MPI_COMM_WORLD, stat_h, ierr_h)

     ! host -> GPU (slot 1 only)
     !$acc update device(rbuf_l(:,:,:,1:1), rbuf_r(:,:,:,1:1))

     ! Unpack received halos into the inner halo slabs (on GPU)
     !$acc parallel loop collapse(3) present(arr, rbuf_l, rbuf_r)
     do is = 1, 3
       do i2 = 1, ngrid2
         do i3 = 1, ngrid3
           arr(i3, i2, 0,            is) = rbuf_l(i3, i2, is, 1)
           arr(i3, i2, local_n1 + 1, is) = rbuf_r(i3, i2, is, 1)
         end do
       end do
     end do
   end subroutine exchange_halo

 ! =====================================================================
 ! WIDE HALO EXCHANGE (two slabs per side). Fills all four R1 halo slabs:
 !   il1 = -1, 0, local_n1+1, local_n1+2.
 ! Slot 1 of the buffers = innermost slab, slot 2 = next-out slab. This
 ! ordering is consistent with exchange_halo (which uses only slot 1), so
 ! the inner-halo content agrees whether the narrow or wide call was the
 ! most recent. Only compute_momentum needs the outer halos.
 ! =====================================================================
   subroutine exchange_halo_wide(arr)
     implicit none
     real*8, intent(inout) :: arr(ngrid3, ngrid2, -1:local_n1+2, 3)
     integer :: i2, i3, is, ierr_h
     integer :: stat_h(MPI_STATUS_SIZE)

     ! Pack two slabs per side (on GPU)
     !$acc parallel loop collapse(3) present(arr, sbuf_l, sbuf_r)
     do is = 1, 3
       do i2 = 1, ngrid2
         do i3 = 1, ngrid3
           sbuf_l(i3, i2, is, 1) = arr(i3, i2, 1,            is)  ! innermost
           sbuf_l(i3, i2, is, 2) = arr(i3, i2, 2,            is)  ! next-out
           sbuf_r(i3, i2, is, 1) = arr(i3, i2, local_n1,     is)  ! innermost
           sbuf_r(i3, i2, is, 2) = arr(i3, i2, local_n1 - 1, is)  ! next-out
         end do
       end do
     end do

     ! GPU -> host (both slots)
     !$acc update self(sbuf_l, sbuf_r)

     call MPI_Sendrecv(sbuf_r, 2*halo_size, MPI_DOUBLE_PRECISION, right_rank, 0, &
                       rbuf_l, 2*halo_size, MPI_DOUBLE_PRECISION, left_rank,  0, &
                       MPI_COMM_WORLD, stat_h, ierr_h)
     call MPI_Sendrecv(sbuf_l, 2*halo_size, MPI_DOUBLE_PRECISION, left_rank,  1, &
                       rbuf_r, 2*halo_size, MPI_DOUBLE_PRECISION, right_rank, 1, &
                       MPI_COMM_WORLD, stat_h, ierr_h)

     ! host -> GPU (both slots)
     !$acc update device(rbuf_l, rbuf_r)

     ! Unpack into four halo slabs (on GPU)
     !$acc parallel loop collapse(3) present(arr, rbuf_l, rbuf_r)
     do is = 1, 3
       do i2 = 1, ngrid2
         do i3 = 1, ngrid3
           arr(i3, i2,  0,           is) = rbuf_l(i3, i2, is, 1)  ! neighbor's local_n1
           arr(i3, i2, -1,           is) = rbuf_l(i3, i2, is, 2)  ! neighbor's local_n1-1
           arr(i3, i2, local_n1 + 1, is) = rbuf_r(i3, i2, is, 1)  ! neighbor's 1
           arr(i3, i2, local_n1 + 2, is) = rbuf_r(i3, i2, is, 2)  ! neighbor's 2
         end do
       end do
     end do
   end subroutine exchange_halo_wide

 ! =====================================================================
 ! DUMP 2D MARGINAL DENSITIES (called with host-current wfr/wfi)
 ! =====================================================================
   subroutine dump_density(nstep)
     implicit none
     integer*8, intent(in) :: nstep

     real*8, allocatable :: rho12_l(:,:,:), rho13_l(:,:,:), rho23_l(:,:,:)
     real*8, allocatable :: rho12_g(:,:,:), rho13_g(:,:,:), rho23_g(:,:,:)
     integer :: ii, jj, kk, ss, ig1, ierr_d
     real*8  :: dd
     character(len=30) :: fname

     allocate(rho12_l(ngrid1, ngrid2, 3), rho12_g(ngrid1, ngrid2, 3))
     allocate(rho13_l(ngrid1, ngrid3, 3), rho13_g(ngrid1, ngrid3, 3))
     allocate(rho23_l(ngrid2, ngrid3, 3), rho23_g(ngrid2, ngrid3, 3))

     rho12_l = 0.0d0;  rho13_l = 0.0d0;  rho23_l = 0.0d0

     do ss = 1, 3
       do ii = 1, local_n1
         ig1 = global_start + ii - 1
         do jj = 1, ngrid2
           do kk = 1, ngrid3
             dd = wfr(kk, jj, ii, ss)**2 + wfi(kk, jj, ii, ss)**2
             rho12_l(ig1, jj, ss) = rho12_l(ig1, jj, ss) + dd
             rho13_l(ig1, kk, ss) = rho13_l(ig1, kk, ss) + dd
             rho23_l(jj,  kk, ss) = rho23_l(jj,  kk, ss) + dd
           end do
         end do
       end do
     end do

     call MPI_Reduce(rho12_l, rho12_g, ngrid1*ngrid2*3, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr_d)
     call MPI_Reduce(rho13_l, rho13_g, ngrid1*ngrid3*3, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr_d)
     call MPI_Reduce(rho23_l, rho23_g, ngrid2*ngrid3*3, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr_d)

     if (myrank == 0) then
       write(fname,'("rho_r1r2_t",I8.8,".bin")') nstep
       open(unit=99, file=trim(fname), status="replace", form="unformatted", access="stream")
       write(99) rho12_g
       close(99)

       write(fname,'("rho_r1r3_t",I8.8,".bin")') nstep
       open(unit=99, file=trim(fname), status="replace", form="unformatted", access="stream")
       write(99) rho13_g
       close(99)

       write(fname,'("rho_r2r3_t",I8.8,".bin")') nstep
       open(unit=99, file=trim(fname), status="replace", form="unformatted", access="stream")
       write(99) rho23_g
       close(99)
     end if

     deallocate(rho12_l, rho13_l, rho23_l)
     deallocate(rho12_g, rho13_g, rho23_g)
   end subroutine dump_density

 ! =====================================================================
 ! COMPUTE MOMENTUM EXPECTATION VALUES AND STD DEVIATIONS
 !
 ! 4th-order central finite differences (leading error O(deltaR^4)):
 !   d/dR    ~= [-f(+2) + 8 f(+1) - 8 f(-1) + f(-2)] / (12 deltaR)
 !   d^2/dR^2 ~= [-f(-2) + 16 f(-1) - 30 f + 16 f(+1) - f(+2)] / (12 deltaR^2)
 !
 !   <p_i>   = integral [Re * d Im/d R_i - Im * d Re/d R_i] dV
 !   <p_i^2> = -integral psi^* (d^2/d R_i^2) psi dV
 !           = -integral [Re * d^2 Re/d R_i^2 + Im * d^2 Im/d R_i^2] dV
 !
 ! Plane-wave eigenvalues: for psi = exp(i k x) on the grid,
 !   4th-order D1 gives k * [1 - (k h)^4 / 30 + ...]
 !   4th-order D2 gives -k^2 * [1 - (k h)^4 / 90 + ...]
 ! At k = 10 and deltaR = 0.012 this is ~7e-6 relative error for <p> and
 ! ~2e-4 for <p^2>, versus ~2.4e-3 and ~1.2e-3 from 2nd-order stencils.
 !
 ! Requires a 2-cell R1 halo (filled by exchange_halo_wide). In R2/R3 the
 ! loop is restricted to (3..ngrid{2,3}-2) so the 5-point stencil fits in
 ! the array; the two skipped outermost interior slabs contribute
 ! negligibly for any wavefunction that is well-localized away from the
 ! domain boundaries.
 ! =====================================================================
   subroutine compute_momentum(nstep)
     implicit none
     integer*8, intent(in) :: nstep

     real*8 :: p1_l, p2_l, p3_l, p1sq_l, p2sq_l, p3sq_l
     real*8 :: mom_loc(6), mom_glb(6)
     real*8 :: dwfr1, dwfi1, dwfr2, dwfi2, dwfr3, dwfi3
     real*8 :: d2wfr1, d2wfi1, d2wfr2, d2wfi2, d2wfr3, d2wfi3
     real*8 :: c12_r1, c12_r2, c12_r3
     real*8 :: c12_r1sq, c12_r2sq, c12_r3sq
     integer :: il1, ig2, ig3, is, ierr_m

     ! Two-cell halo exchange so the 4th-order stencils have neighbor data
     call exchange_halo_wide(wfr)
     call exchange_halo_wide(wfi)

     c12_r1   = 1.d0 / (12.d0 * deltaR1)
     c12_r2   = 1.d0 / (12.d0 * deltaR2)
     c12_r3   = 1.d0 / (12.d0 * deltaR3)
     c12_r1sq = 1.d0 / (12.d0 * deltaR1**2)
     c12_r2sq = 1.d0 / (12.d0 * deltaR2**2)
     c12_r3sq = 1.d0 / (12.d0 * deltaR3**2)

     p1_l = 0.d0;  p2_l = 0.d0;  p3_l = 0.d0
     p1sq_l = 0.d0;  p2sq_l = 0.d0;  p3sq_l = 0.d0

     !$acc parallel loop collapse(3) present(wfr, wfi) &
     !$acc reduction(+:p1_l,p2_l,p3_l,p1sq_l,p2sq_l,p3sq_l)
     do il1 = 1, local_n1
       do ig2 = 3, ngrid2-2
         do ig3 = 3, ngrid3-2
           !$acc loop seq
           do is = 1, 3
             ! ----- 4th-order first derivatives -----
             dwfr1 = (-wfr(ig3,ig2,il1+2,is) + 8.d0*wfr(ig3,ig2,il1+1,is) &
                      - 8.d0*wfr(ig3,ig2,il1-1,is) + wfr(ig3,ig2,il1-2,is)) * c12_r1
             dwfi1 = (-wfi(ig3,ig2,il1+2,is) + 8.d0*wfi(ig3,ig2,il1+1,is) &
                      - 8.d0*wfi(ig3,ig2,il1-1,is) + wfi(ig3,ig2,il1-2,is)) * c12_r1
             dwfr2 = (-wfr(ig3,ig2+2,il1,is) + 8.d0*wfr(ig3,ig2+1,il1,is) &
                      - 8.d0*wfr(ig3,ig2-1,il1,is) + wfr(ig3,ig2-2,il1,is)) * c12_r2
             dwfi2 = (-wfi(ig3,ig2+2,il1,is) + 8.d0*wfi(ig3,ig2+1,il1,is) &
                      - 8.d0*wfi(ig3,ig2-1,il1,is) + wfi(ig3,ig2-2,il1,is)) * c12_r2
             dwfr3 = (-wfr(ig3+2,ig2,il1,is) + 8.d0*wfr(ig3+1,ig2,il1,is) &
                      - 8.d0*wfr(ig3-1,ig2,il1,is) + wfr(ig3-2,ig2,il1,is)) * c12_r3
             dwfi3 = (-wfi(ig3+2,ig2,il1,is) + 8.d0*wfi(ig3+1,ig2,il1,is) &
                      - 8.d0*wfi(ig3-1,ig2,il1,is) + wfi(ig3-2,ig2,il1,is)) * c12_r3

             p1_l = p1_l + wfr(ig3,ig2,il1,is)*dwfi1 - wfi(ig3,ig2,il1,is)*dwfr1
             p2_l = p2_l + wfr(ig3,ig2,il1,is)*dwfi2 - wfi(ig3,ig2,il1,is)*dwfr2
             p3_l = p3_l + wfr(ig3,ig2,il1,is)*dwfi3 - wfi(ig3,ig2,il1,is)*dwfr3

             ! ----- 4th-order second derivatives -----
             d2wfr1 = (-wfr(ig3,ig2,il1-2,is) + 16.d0*wfr(ig3,ig2,il1-1,is) &
                      - 30.d0*wfr(ig3,ig2,il1,is)   + 16.d0*wfr(ig3,ig2,il1+1,is) &
                      - wfr(ig3,ig2,il1+2,is)) * c12_r1sq
             d2wfi1 = (-wfi(ig3,ig2,il1-2,is) + 16.d0*wfi(ig3,ig2,il1-1,is) &
                      - 30.d0*wfi(ig3,ig2,il1,is)   + 16.d0*wfi(ig3,ig2,il1+1,is) &
                      - wfi(ig3,ig2,il1+2,is)) * c12_r1sq
             d2wfr2 = (-wfr(ig3,ig2-2,il1,is) + 16.d0*wfr(ig3,ig2-1,il1,is) &
                      - 30.d0*wfr(ig3,ig2,il1,is)   + 16.d0*wfr(ig3,ig2+1,il1,is) &
                      - wfr(ig3,ig2+2,il1,is)) * c12_r2sq
             d2wfi2 = (-wfi(ig3,ig2-2,il1,is) + 16.d0*wfi(ig3,ig2-1,il1,is) &
                      - 30.d0*wfi(ig3,ig2,il1,is)   + 16.d0*wfi(ig3,ig2+1,il1,is) &
                      - wfi(ig3,ig2+2,il1,is)) * c12_r2sq
             d2wfr3 = (-wfr(ig3-2,ig2,il1,is) + 16.d0*wfr(ig3-1,ig2,il1,is) &
                      - 30.d0*wfr(ig3,ig2,il1,is)   + 16.d0*wfr(ig3+1,ig2,il1,is) &
                      - wfr(ig3+2,ig2,il1,is)) * c12_r3sq
             d2wfi3 = (-wfi(ig3-2,ig2,il1,is) + 16.d0*wfi(ig3-1,ig2,il1,is) &
                      - 30.d0*wfi(ig3,ig2,il1,is)   + 16.d0*wfi(ig3+1,ig2,il1,is) &
                      - wfi(ig3+2,ig2,il1,is)) * c12_r3sq

             ! <p_i^2> = -<psi | d^2/dR_i^2 | psi>
             p1sq_l = p1sq_l - wfr(ig3,ig2,il1,is)*d2wfr1 - wfi(ig3,ig2,il1,is)*d2wfi1
             p2sq_l = p2sq_l - wfr(ig3,ig2,il1,is)*d2wfr2 - wfi(ig3,ig2,il1,is)*d2wfi2
             p3sq_l = p3sq_l - wfr(ig3,ig2,il1,is)*d2wfr3 - wfi(ig3,ig2,il1,is)*d2wfi3
           end do
         end do
       end do
     end do

     mom_loc(1) = p1_l;    mom_loc(2) = p2_l;    mom_loc(3) = p3_l
     mom_loc(4) = p1sq_l;  mom_loc(5) = p2sq_l;  mom_loc(6) = p3sq_l
     call MPI_Allreduce(mom_loc, mom_glb, 6, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_WORLD, ierr_m)
     mom_glb = mom_glb * deltaR1 * deltaR2 * deltaR3

     ! Write: step  <p1> <p2> <p3>  sigma(p1) sigma(p2) sigma(p3)
     if (myrank == 0) then
       write(32,*) nstep, mom_glb(1), mom_glb(2), mom_glb(3), &
                   sqrt(max(0.d0, mom_glb(4) - mom_glb(1)**2)), &
                   sqrt(max(0.d0, mom_glb(5) - mom_glb(2)**2)), &
                   sqrt(max(0.d0, mom_glb(6) - mom_glb(3)**2))
     end if
   end subroutine compute_momentum

 end program exact
