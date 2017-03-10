! stefan problem analytical solution from Book
!"Mathematical modeling of Melting and freezing process"
! Code by Ayurzana.B
      program stefan

      parameter (m=100, t=200, dis=0.1 ) !grid size, time, slab (m)
      parameter (tm=0., th=25.0, tc=0.) !melt,hot,cold
      parameter (tdl=0.143e-6, tds=1.1819e-6) !thermal Dif (m2/s-1)
      parameter (La=333.4,  pi=3.1415926 ) !kJ/kg and pi
      parameter (chl=4.1868, chs=0.138) !specific heat (J/gCels)
      real tem(t,m), xs(t) ! temperature and melt front
      real dx,dt,time,lamda   ! spacing, time, parameter
      real St
      character*6 rr
	  
      open(1,file='interface.dat')
      open(2,file='temp2.dat')
      open(7,file='flowan.dat')
	  
      St=chl*(th-tm)/La
      !lamda=0.706*sqrt(St)*(1.-0.21*(0.5642*St)**(0.93-0.15*St))
      call liamda(St,lamda)
      t_t=dis**2/(4.*lamda**2*tdl)
      dt=t_t/t
      dx=dis/m
	  
      t_p1=0.01**2/(4.*lamda**2*tdl) !1cm
      t_p3=0.03**2/(4.*lamda**2*tdl) !3cm
      t_p5=0.05**2/(4.*lamda**2*tdl) !5cm
      t_p9=0.09**2/(4.*lamda**2*tdl) !9cm
	  
      print*, t_t/3600, 50*dt/3600, 100*dt/3600, 150*dt/3600, 200*dt/3600
      xs=0.
      tem=0.
	  
      do i=1,t
         time=i*dt
         xs(i)=2.*lamda*sqrt(tdl*time) ! free boundary
         do j=1,m
              if(j*dx.le.xs(i)) then
            tem(i,j)=th+(tm-th)*erf((j*dx)/(sqrt(4.*tdl*time)))
     &                         /erf(lamda)
              else
!         if(j*dx.ge.xs(i)) tem(i,j)=tm
            tem(i,j)=tm*erf((j*dx)/(sqrt(4.*tdl*time)))
     &                         /erf(lamda)
              end if
		 
         if(int(time).eq.int(t_p1)) then
      write(rr,300) 100000+int(t_p1) 
300   format(i6)
         open(3,file='1cm'//TRIM(rr)//'.dat')
         write(3,*) j*dx, tem(i,j)
         endif
		 
         if(int(time).eq.int(t_p3)) then
      write(rr,300) 100000+int(t_p3) 
         open(4,file='3cm'//TRIM(rr)//'.dat')
         write(4,*) j*dx, tem(i,j)
         endif
		 
         if(int(time).eq.int(t_p5)) then
      write(rr,300) 100000+int(t_p5) 
         open(5,file='5cm'//TRIM(rr)//'.dat')
         write(5,*) j*dx, tem(i,j)
         endif
		 
         if(int(time).eq.int(t_p9)) then
      write(rr,300) 100000+int(t_p9) 
         open(6,file='9cm'//TRIM(rr)//'.dat')
         write(6,*) j*dx, tem(i,j)
         endif
	 
         end do
		 
         write(1,*) i, xs(i)/dx !time/3600.
         write(7,*) (i*dt)/3600., tem(i,90)			 


      end do
	  
      do i=1,t
      write(2,200) (tem(i,j), j=1,m)
      enddo
  200 format(1000e11.3)
	  
      print*, St, lamda, t_p1, t_p3, t_p5, t_t
      
	  
      close (1)
      close (2)
      close (7)
      stop
      end program stefan
!------------------------------------------------------------!
      subroutine liamda (St,lamda)
      parameter (pi=3.1415926)
      parameter (iter=20)
      real lamda

      lamda=0.1 !initial guess
      do i=1,iter
      funct=St/sqrt(pi)-lamda*exp(lamda**2)*erf(lamda)
      derivative=-exp(lamda**2)*erf(lamda)-2.*lamda**2
     &           *exp(lamda**2)*erf(lamda)-lamda*2./sqrt(pi)
      lamda=lamda-funct/derivative
      end do

      return
      end
