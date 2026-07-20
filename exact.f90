program exact
   implicit none
 ! This program solves the TDSE of a nucleus on PESs of 3 coupled states
 
 !--------------- grid size and domain --------------------------------    
   integer, parameter :: ngrid1 = 1000, ngrid2 = 1000
   real*8, parameter :: x1min = -4.0d0, x1max =  8.0d0
   real*8, parameter :: x2min = -4.0d0, x2max =  8.0d0
 
 !--------------- model parameters ------------------------------------
   real*8, parameter :: w1    = 0.25d0
   real*8, parameter :: w2    = 0.025d0
   real*8, parameter :: delta = 0.01d0
   real*8, parameter :: epsil = 0.00d0
   real*8, parameter :: c     = 0.025d0
   real*8, parameter :: pmass = 1.845d3
 
 !--------------- integration parameters ------------------------------
   real*8, parameter :: deltat   = 0.01d0
   real*8, parameter :: tstepmax = 4.0d4
 
 !--------------- initial conditions ----------------------------------
   real*8, parameter :: R1bar = -1.0d0, R2bar = 0.0d0
   real*8, parameter :: P1bar = 10.0d0, P2bar = 10.0d0
 
 !--------------- other parameters ------------------------------------
   real*8, parameter :: pi     = 3.14159265359d0
   real*8, parameter :: alpha1 = 6.0d0, alpha2 = 6.0d0
 
 !--------------- wavefunction arrays ---------------------------------
   real*8 :: wfr(3,ngrid1,ngrid2)
   real*8 :: wfi(3,ngrid1,ngrid2)
   real*8 :: dwfrdt(3,ngrid1,ngrid2)
   real*8 :: dwfidt(3,ngrid1,ngrid2)
   real*8 :: V(3,ngrid1,ngrid2)
 
 !--------------- diagnostics arrays ----------------------------------
   real*8 :: pop(3)
 
 !--------------- workspace variables ---------------------------------
   real*8 :: deltaR1, deltaR2
   real*8 :: prefac1, prefac2
   real*8 :: gauss1, gauss2
   real*8 :: R1, R2, R1disp, R2disp
   real*8 :: PbarRdisp, cosPbarRdisp, sinPbarRdisp
   real*8 :: kedenom1, kedenom2, halfdeltat
   real*8 :: norm
   real*8 :: start, finish
   integer*8 :: igrid1, igrid2, istate, n
 
 !--------------- precompute constants --------------------------------
   deltaR1    = (x1max - x1min) / dble(ngrid1)
   deltaR2    = (x2max - x2min) / dble(ngrid2)
   prefac1    = sqrt(sqrt(2.0d0 * alpha1 / pi))
   prefac2    = sqrt(sqrt(2.0d0 * alpha2 / pi))
   kedenom1   = 2.0d0 * deltaR1**2 * pmass
   kedenom2   = 2.0d0 * deltaR2**2 * pmass
   halfdeltat = 0.5d0 * deltat
 
   call cpu_time(start)
 
 !--------------- initialize wavefunction and potentials --------------
   do igrid1 = 1, ngrid1
     do igrid2 = 1, ngrid2
 
       ! boundary conditions
       if (igrid1==1 .or. igrid1==ngrid1 .or. igrid2==1 .or. igrid2==ngrid2) then
         wfr(:,igrid1,igrid2) = 0.0d0
         wfi(:,igrid1,igrid2) = 0.0d0
       else
         R1 = x1min + (dble(igrid1)-0.5d0)*deltaR1
         R2 = x2min + (dble(igrid2)-0.5d0)*deltaR2
 
         R1disp      = R1 - R1bar
         R2disp      = R2 - R2bar
         gauss1      = exp(-alpha1 * R1disp**2)
         gauss2      = exp(-alpha2 * R2disp**2)
         PbarRdisp   = P1bar*R1disp + P2bar*R2disp
         cosPbarRdisp= cos(PbarRdisp)
         sinPbarRdisp= sin(PbarRdisp)
 
         wfr(1,igrid1,igrid2) = prefac1*prefac2*gauss1*gauss2*cosPbarRdisp
         wfi(1,igrid1,igrid2) = prefac1*prefac2*gauss1*gauss2*sinPbarRdisp
         wfr(2,igrid1,igrid2)= 0.0d0
         wfi(2,igrid1,igrid2)= 0.0d0
         wfr(3,igrid1,igrid2)= 0.0d0
         wfi(3,igrid1,igrid2)= 0.0d0
       end if
 
       V(1,igrid1,igrid2) = -w1 * R1
       do istate = 2, 3
         V(istate,igrid1,igrid2) = w2 * R1 - (istate-1) * delta
       end do
 
     end do
   end do
 
 !--------------- write initial wavefunction --------------------------
   open(unit=21, file="wf_init.dat", status="replace")
   do igrid1 = 1, ngrid1
     do igrid2 = 1, ngrid2
       R1 = x1min + (dble(igrid1)-0.5d0)*deltaR1
       R2 = x2min + (dble(igrid2)-0.5d0)*deltaR2
       write(21,*) R1, R2, wfr(1,igrid1,igrid2), wfi(1,igrid1,igrid2)
     end do
   end do
   close(21)
 
 !--------------- compute initial norm & populations ------------------
   norm = 0.0d0
   pop  = 0.0d0
   do igrid1 = 2, ngrid1-1
     do igrid2 = 2, ngrid2-1
       do istate = 1, 3
         norm = norm + wfr(istate,igrid1,igrid2)**2 + wfi(istate,igrid1,igrid2)**2
         pop(istate) = pop(istate) + wfr(istate,igrid1,igrid2)**2 + wfi(istate,igrid1,igrid2)**2
       end do
     end do
   end do
   norm = norm * deltaR1 * deltaR2
   pop  = pop  * deltaR1 * deltaR2
 
   open(unit=30, file="norm.dat", status="replace")
   open(unit=31, file="pop.dat",  status="replace")
   write(30,*) 0_8, norm
   write(31,*) 0_8, pop
 
   n = 1_8
 
 !--------------- main propagation loop -------------------------------
   do while (n < tstepmax)
 
     ! compute dwfrdt
     do igrid1 = 2, ngrid1-1
       do igrid2 = 2, ngrid2-1
         do istate = 1, 3
           dwfrdt(istate,igrid1,igrid2) = -(( &
                wfi(istate,igrid1-1,igrid2) - 2.d0*wfi(istate,igrid1,igrid2) + wfi(istate,igrid1+1,igrid2) ) / kedenom1 + &
                ( wfi(istate,igrid1,igrid2-1) - 2.d0*wfi(istate,igrid1,igrid2) + wfi(istate,igrid1,igrid2+1) ) / kedenom2 ) + &
                V(istate,igrid1,igrid2) * wfi(istate,igrid1,igrid2)
         end do
       end do
     end do
     do igrid1 = 2, ngrid1-1
       do igrid2 = 2, ngrid2-1
         R2 = x2min + (dble(igrid2)-0.5d0)*deltaR2
         do istate = 2, 3
           dwfrdt(1,igrid1,igrid2)    = dwfrdt(1,igrid1,igrid2)    + c * R2 * wfi(istate,igrid1,igrid2)
           dwfrdt(istate,igrid1,igrid2)= dwfrdt(istate,igrid1,igrid2) + c * R2 * wfi(1,igrid1,igrid2)
         end do
       end do
     end do
 
     ! propagate wfr by dt/2
     wfr = wfr + halfdeltat * dwfrdt
 
     ! compute dwfidt
     do igrid1 = 2, ngrid1-1
       do igrid2 = 2, ngrid2-1
         do istate = 1, 3
           dwfidt(istate,igrid1,igrid2) = (( &
                wfr(istate,igrid1-1,igrid2) - 2.d0*wfr(istate,igrid1,igrid2) + wfr(istate,igrid1+1,igrid2) ) / kedenom1 + &
                ( wfr(istate,igrid1,igrid2-1) - 2.d0*wfr(istate,igrid1,igrid2) + wfr(istate,igrid1,igrid2+1) ) / kedenom2 ) - &
                V(istate,igrid1,igrid2) * wfr(istate,igrid1,igrid2)
         end do
       end do
     end do
     do igrid1 = 2, ngrid1-1
       do igrid2 = 2, ngrid2-1
         R2 = x2min + (dble(igrid2)-0.5d0)*deltaR2
         do istate = 2, 3
           dwfidt(1,igrid1,igrid2)    = dwfidt(1,igrid1,igrid2)    - c * R2 * wfr(istate,igrid1,igrid2)
           dwfidt(istate,igrid1,igrid2)= dwfidt(istate,igrid1,igrid2) - c * R2 * wfr(1,igrid1,igrid2)
         end do
       end do
     end do
 
     ! propagate wfi by dt
     wfi = wfi + deltat * dwfidt
 
     ! compute dwfrdt again
     do igrid1 = 2, ngrid1-1
       do igrid2 = 2, ngrid2-1
         do istate = 1, 3
           dwfrdt(istate,igrid1,igrid2) = -(( &
                wfi(istate,igrid1-1,igrid2) - 2.d0*wfi(istate,igrid1,igrid2) + wfi(istate,igrid1+1,igrid2) ) / kedenom1 + &
                ( wfi(istate,igrid1,igrid2-1) - 2.d0*wfi(istate,igrid1,igrid2) + wfi(istate,igrid1,igrid2+1) ) / kedenom2 ) + &
                V(istate,igrid1,igrid2) * wfi(istate,igrid1,igrid2)
         end do
       end do
     end do
     do igrid1 = 2, ngrid1-1
       do igrid2 = 2, ngrid2-1
         R2 = x2min + (dble(igrid2)-0.5d0)*deltaR2
         do istate = 2, 3
           dwfrdt(1,igrid1,igrid2)    = dwfrdt(1,igrid1,igrid2)    + c * R2 * wfi(istate,igrid1,igrid2)
           dwfrdt(istate,igrid1,igrid2)= dwfrdt(istate,igrid1,igrid2) + c * R2 * wfi(1,igrid1,igrid2)
         end do
       end do
     end do
 
     ! propagate wfr by dt/2 again
     wfr = wfr + halfdeltat * dwfrdt
 
     ! diagnostics every 100 steps
     if (mod(n,100)==0) then
       norm = 0.0d0
       pop  = 0.0d0
       do igrid1 = 2, ngrid1-1
         do igrid2 = 2, ngrid2-1
           do istate = 1, 3
             norm = norm + wfr(istate,igrid1,igrid2)**2 + wfi(istate,igrid1,igrid2)**2
             pop(istate) = pop(istate) + wfr(istate,igrid1,igrid2)**2 + wfi(istate,igrid1,igrid2)**2
           end do
         end do
       end do
       norm = norm * deltaR1 * deltaR2
       pop  = pop  * deltaR1 * deltaR2
       write(30,*) n, norm
       write(31,*) n, pop
       if (mod(n,400)==0) call dump_density(n, wfr, wfi, ngrid1, ngrid2)
     end if
 
     n = n + 1_8
   end do
 
   call cpu_time(finish)
   print*, "Program completes in ", finish-start
 
 contains
 
   subroutine dump_density(nstep, wfr, wfi, ng1, ng2)
     implicit none
     integer*8, intent(in) :: nstep
     integer,   intent(in) :: ng1, ng2
     integer :: i, j, k
     real*8,   intent(in) :: wfr(3,ng1,ng2), wfi(3,ng1,ng2)
     real*8, allocatable  :: rho(:,:)
     character(len=20)    :: fname
 
     allocate(rho(ng1,ng2))
     rho = 0.0d0
     do k = 1, 3
       do i = 1, ng1
         do j = 1, ng2
           rho(i,j) = rho(i,j) + wfr(k,i,j)**2 + wfi(k,i,j)**2
         end do
       end do
     end do
 
     write(fname,'("rho_t",I8.8,".bin")') nstep
     open(unit=99, file=trim(fname), status="replace", form="unformatted", access="stream")
     write(99) rho
     close(99)
 
     deallocate(rho)
   end subroutine dump_density
 
 end program exact
 