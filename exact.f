      program exact
      implicit none
c This program sovles the TDSE of a nucleus on PESs of 9 coupled states

c=============== grid size and domain=================================
      integer ngrid1,ngrid2
      real*8 x1min, x2min, x1max, x2max
      parameter (ngrid1=1000,ngrid2=1000)
      parameter (x1min=-4.0d0,x2min=-4.0d0,x1max=8.0d0,x2max=8.0d0)

c=============== model parameters ====================================
      real*8 w1, w2, delta, c, pmass, epsil
             !slope of state1
             !slope of state2-9
             !spacing between state2-9
             !linear coupling constant
             !nuclear mass
      parameter (w1=0.25d0,w2=0.025d0,delta=0.0005d0,epsil=0.00)
      parameter (c=0.025d0,pmass=1.845d3)
c=============== integration parameters ==============================
      real*8 deltat,tstepmax
             !nuclear simulation time step
             !maximum number of nuclear time steps within a simulation     
      parameter (deltat=0.01d0,tstepmax=4.0d4)
c=============== initial conditions ==================================
      real*8 R1bar,R2bar,P1bar,P2bar
      parameter(R1bar=-1.0d0,R2bar=0.0d0,P1bar=10.0d0,P2bar=10.0d0)
c=============== other parameters ====================================
      real*8 pi,alpha1,alpha2 !width
      parameter (pi=3.14159265359d0,alpha1=6.0d0,alpha2=6.0d0)
      
      real*8 wfr(9,ngrid1,ngrid2) !real-part of wavefunction
      real*8 wfi(9,ngrid1,ngrid2) !imginary-part of wavefunction

      real*8 dwfrdt(9,ngrid1,ngrid2) !time-derivative of real-part of wavefunction
      real*8 dwfidt(9,ngrid1,ngrid2) !time-derivative of imginary-part of wavefunction

      real*8 V(9,ngrid1,ngrid2) !potential energy
      
      real*8 deltaR1,deltaR2 !spacing of grids on x1/x2 direction
      real*8 prefac1,prefac2 !prefactor of gaussian function
      real*8 gauss1, gauss2
      real*8 R1disp, R2disp !distance from central gaussian
      real*8 n !time step count
      real*8 halfdeltat !half nuclear time step
      real*8 R1,R2 !position of nulceus
      
      real*8 cosPbarRdisp,sinPbarRdisp
      real*8 PbarRdisp
      real*8 kedenom1,kedenom2 !denominator of kinetic energy operator      

      real*8 norm        !norm of wf

      real*8 pop(9)     !population of diabatic states

      integer*8 igrid1,igrid2,istate
      real*8 start,finish

      deltaR1=(x1max-x1min)/dble(ngrid1)
      deltaR2=(x2max-x2min)/dble(ngrid2)
      prefac1=sqrt(sqrt(2.0d0*alpha1/pi))
      prefac2=sqrt(sqrt(2.0d0*alpha2/pi))

      call cpu_time(start)            
      do igrid1=1,ngrid1
         do igrid2=1,ngrid2
           !boundary conditions
           if (igrid1==1 .or. igrid1==ngrid1 .or. igrid2==1 
     c        .or.igrid2==ngrid2) then
              wfr(1,igrid1,igrid2)=0.0d0
              wfi(1,igrid1,igrid2)=0.0d0
           else
              R1=x1min+(dble(igrid1)-0.5d0)*deltaR1
              R2=x2min+(dble(igrid2)-0.5d0)*deltaR2
              R1disp=R1-R1bar
              R2disp=R2-R2bar
              gauss1=exp(-1.0d0*alpha1*R1disp*R1disp)
              gauss2=exp(-1.0d0*alpha2*R2disp*R2disp)
              PbarRdisp=P1bar*R1disp+P2bar*R2disp
              cosPbarRdisp=cos(PbarRdisp)
              sinPbarRdisp=sin(PbarRdisp)
          
              wfr(1,igrid1,igrid2)=
     c        prefac1*prefac2*gauss1*gauss2*cosPbarRdisp
              wfi(1,igrid1,igrid2)=
     c        prefac1*prefac2*gauss1*gauss2*sinPbarRdisp
           end if
        
           do istate = 2,9
                wfr(istate,igrid1,igrid2)=0.0d0
                wfi(istate,igrid1,igrid2)=0.0d0
           enddo
                
           V(1,igrid1,igrid2)=-1.0d0*w1*R1
           do istate = 2,5
                V(istate,igrid1,igrid2)=w2*R1-(istate-1)*delta
           enddo
                
           do istate = 6,9
                V(istate,igrid1,igrid2)=w2*R1-(istate-1)*delta-epsil
           enddo
         enddo
      enddo

      open(21,file="wf_init.dat") !output of initial wavefunction
      do igrid1=1,ngrid1
         do igrid2=1,ngrid2
          R1=x1min+(dble(igrid1)-0.5d0)*deltaR1
          R2=x2min+(dble(igrid2)-0.5d0)*deltaR2
          write(21,9999) R1,R2,wfr(1,igrid1,igrid2),wfi(1,igrid1,igrid2)   
         enddo
      enddo
!========write initial norm into wf_init.dat=================
      norm=0.0
      do igrid1=2,ngrid1-1
         do igrid2=2,ngrid2-1
            do istate=1,9
               norm=norm+
     c           wfr(istate,igrid1,igrid2)*wfr(istate,igrid1,igrid2)
     c          +wfi(istate,igrid1,igrid2)*wfi(istate,igrid1,igrid2)
            enddo
         enddo
      enddo

      norm = norm*(deltaR1*deltaR2)

      write(21,*) "initial norm:", norm
!==========
      pop(1)=0.0
      do igrid1=2,ngrid1-1
         do igrid2=2,ngrid2-1
            pop(1)=pop(1)+wfr(1,igrid1,igrid2)*wfr(1,igrid1,igrid2)+
     c           wfi(1,igrid1,igrid2)*wfi(1,igrid1,igrid2)
         enddo
      enddo

      pop(1)=pop(1)*(deltaR1*deltaR2)

      write(21,*) "initial population on diabat 1:", pop(1)
      close(20)
      close(21)
!===========norm.dat===============================================
      open(30,file="norm.dat") !output of norm
      open(31,file="pop.dat") !output of population 
c=============== start propogation ================================
      n=1
      kedenom1=2.0*deltaR1*deltaR1*pmass
      kedenom2=2.0*deltaR2*deltaR2*pmass
      halfdeltat=0.5d0*deltat
      do while (n<tstepmax)
c---------- update wfr in 0.5dt -----------------------------------
!--- compute diabatic Hc ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              do istate=1,9          
                   dwfrdt(istate,igrid1,igrid2)=
     c             -1.0d0*(wfi(istate,igrid1-1,igrid2)
     c             -2.0*wfi(istate,igrid1,igrid2)
     c             +wfi(istate,igrid1+1,igrid2))/kedenom1
     c             -1.0d0*(wfi(istate,igrid1,igrid2-1)
     c             -2.0*wfi(istate,igrid1,igrid2)
     c             +wfi(istate,igrid1,igrid2+1))/kedenom2
     c             +V(istate,igrid1,igrid2)*wfi(istate,igrid1,igrid2)
              enddo
           enddo
        enddo    
!--- compute couplings in Hc ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              R2=x2min+(dble(igrid2)-0.5d0)*deltaR2
              do istate=2,9
                dwfrdt(1,igrid1,igrid2)=
     c               dwfrdt(1,igrid1,igrid2)
     c               +c*R2*wfi(istate,igrid1,igrid2)
                dwfrdt(istate,igrid1,igrid2)=
     c               dwfrdt(istate,igrid1,igrid2)
     c               +c*R2*wfi(1,igrid1,igrid2)
              enddo
           enddo
        enddo
        
!--- propagate wfr ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              do istate=1,9
                 wfr(istate,igrid1,igrid2)=
     c                wfr(istate,igrid1,igrid2)
     c                +halfdeltat*dwfrdt(istate,igrid1,igrid2)
              enddo
           enddo
        enddo
c---------- update wfi in dt --------------------------------------
!--- compute diabatic Hc ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              do istate=1,9
                   dwfidt(istate,igrid1,igrid2)=
     c            (wfr(istate,igrid1-1,igrid2)
     c            -2.0*wfr(istate,igrid1,igrid2)
     c            +wfr(istate,igrid1+1,igrid2))/kedenom1
     c            +(wfr(istate,igrid1,igrid2-1)
     c            -2.0*wfr(istate,igrid1,igrid2)
     c            +wfr(istate,igrid1,igrid2+1))/kedenom2
     c          -1.0d0*V(istate,igrid1,igrid2)*wfr(istate,igrid1,igrid2)
              enddo
           enddo
        enddo
!--- compute couplings in Hc ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              R2=x2min+(dble(igrid2)-0.5d0)*deltaR2
              do istate=2,9
                dwfidt(1,igrid1,igrid2)=
     c               dwfidt(1,igrid1,igrid2)
     c               -c*R2*wfr(istate,igrid1,igrid2)
                dwfidt(istate,igrid1,igrid2)=
     c               dwfidt(istate,igrid1,igrid2)
     c               -c*R2*wfr(1,igrid1,igrid2)
              enddo
           enddo
        enddo

!--- propagate wfi ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              do istate=1,9
                 wfi(istate,igrid1,igrid2)=
     c               wfi(istate,igrid1,igrid2)
     c               +deltat*dwfidt(istate,igrid1,igrid2)
              enddo
           enddo
        enddo
c---------- update wfr in the other dt/2 ---------------------------------------
!--- compute diabatic Hc ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              do istate=1,9
                   dwfrdt(istate,igrid1,igrid2)=
     c               -1.0d0*(wfi(istate,igrid1-1,igrid2)
     c               -2.0*wfi(istate,igrid1,igrid2)
     c               +wfi(istate,igrid1+1,igrid2))/kedenom1
     c               -1.0d0*(wfi(istate,igrid1,igrid2-1)
     c               -2.0*wfi(istate,igrid1,igrid2)
     c               +wfi(istate,igrid1,igrid2+1))/kedenom2
     c               +V(istate,igrid1,igrid2)*wfi(istate,igrid1,igrid2)
              enddo
           enddo
        enddo
!--- compute couplings in Hc ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              R2=x2min+(dble(igrid2)-0.5d0)*deltaR2
              do istate=1,9
                dwfrdt(1,igrid1,igrid2)=
     c             dwfrdt(1,igrid1,igrid2)
     c             +c*R2*wfi(istate,igrid1,igrid2)
                dwfrdt(istate,igrid1,igrid2)=
     c             dwfrdt(istate,igrid1,igrid2)
     c             +c*R2*wfi(1,igrid1,igrid2)
              enddo
           enddo
        enddo

!--- propagate wfr ---
        do igrid1=2,ngrid1-1
           do igrid2=2,ngrid2-1
              do istate=1,9
                 wfr(istate,igrid1,igrid2)=
     c               wfr(istate,igrid1,igrid2)
     c               +halfdeltat*dwfrdt(istate,igrid1,igrid2)
              enddo
           enddo
        enddo
        !write norm into norm.dat
        if (mod(n,1.0d2)==0) then
           norm=0.0
           do istate = 1,9
                pop(istate)=0
           enddo

           do igrid1=2,ngrid1-1
              do igrid2=2,ngrid2-1
                 do istate = 1,9
                   pop(istate)=pop(istate)+
     c             wfr(istate,igrid1,igrid2)*wfr(istate,igrid1,igrid2)
     c            +wfi(istate,igrid1,igrid2)*wfi(istate,igrid1,igrid2)
                 enddo
              enddo
           enddo
           
           do istate = 1,9
                pop(istate)=pop(istate)*deltaR1*deltaR2
           enddo
           norm = sum(pop)
           write(30,9999) n,norm
           write(31,9999) n,pop
        endif
        !End of a time step
        n = n + 1
      enddo
      close(30)
      close(31)
c=============== compute norm ======================================
      open(23,file="out.dat")
      norm=0.0
      do igrid1=2,ngrid1-1
         do igrid2=2,ngrid2-1
            do istate=1,9
               norm=norm+
     c           wfr(istate,igrid1,igrid2)*wfr(istate,igrid1,igrid2)
     c          +wfi(istate,igrid1,igrid2)*wfi(istate,igrid1,igrid2)
            enddo
         enddo
      enddo

      norm = norm*(deltaR1*deltaR2)
      
      write(23,*) "final norm:", norm
 
      call cpu_time(finish)
      
      print*, "Program completes in ",finish-start
 9999 format(100f16.7)  
      end
