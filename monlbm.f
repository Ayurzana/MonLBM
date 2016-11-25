! MonLBM
! computer code for flow around square obstacle
      program obstacle
      include 'param.h'
!      open(2,file='uvfield.dat')
      open(3,file='uvelx.dat')
      open(4,file='vvelx.dat')  ! main result
      open(8,file='timeu.dat')  ! time with u and v
      open(19,file='uvely.dat')
      open(20,file='vvely.dat')
! Initialazation 
      alpha=(1.0/ome-0.5)/3.0
      Re=uo*20/alpha         ! Lattice Reynolds
      print *, 'Re=', Re, 'alpha', alpha
      pause
      w(:)=(/4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36./)
!*************************************************************************
      cx(:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
      cy(:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
!*************************************************************************
      do j=0,m
      do i=0,n
      rho(i,j)=rhoo
      p(i,j)=rho(i,j)/3.
      do k=1,8
      f(k,i,j)=w(k)*rho(i,j)
      enddo
      u(i,j)=0.0
      v(i,j)=0.0
      tau(i,j)=1./ome         ! initial relaxation
      end do
      end do
      do j=1,m-1              ! changed to channel flow
      u(0,j)=uo
      v(0,j)=0.0
      end do
!  ************************************************************************
      write(8,*) Re
! main loop
      ii=0
      do kk=1,mstep
      ii=ii+1
      call collesion
      call streaming
      call sfbound
      call rhouv
!      call subgrid
      call vorti
      call force
	  
!      if (ii.eq.100) then
!      call resul
!      ii=0
!      endif

      print *, tau(130,m/2),v(0,m/2),rho(130,m/2),u(n,m/2),v(100,10)
      write(8,*) kk, cd, cl, cp, v(160,140)
      END DO
! end of the main loop
      call resul
      stop
      end
! end of the main program
! ************************************************************************!
      subroutine subgrid
      include 'param.h'
      cc=0.0002                ! smagorinsky constant 0.0009,0.0002,0.00015

      alpha=(1.0/ome-0.5)/3.0	  
      DO i=0,n
      DO j=0,m
      poten=0.0                ! nonequilibirium stress summation      
      t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
      DO k=0,8
      t2=u(i,j)*cx(k)+v(i,j)*cy(k)
      feq(k,i,j)=w(k)*rho(i,j)*(1.+3.0*t2+4.50*t2*t2-1.50*t1)
      poten=poten+cx(k)*cy(k)*(f(k,i,j)-feq(k,i,j))
      enddo
      smag=sqrt(alpha**2+18.*cc*sqrt(poten*poten))-alpha
      smag=abs(smag/(6.*cc))
      tau(i,j)=3.*(alpha+cc*smag)+0.5
!      if(tau(i,j).gt.1.0) tau(i,j)=0.52
      enddo
      enddo
      return
      end
! ************************************************************************!
      subroutine collesion
      include 'param.h'
      DO i=0,n
      DO j=0,m
      t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
      DO k=0,8
      omega=1/tau(i,j)
      t2=u(i,j)*cx(k)+v(i,j)*cy(k)
      feq(k,i,j)=w(k)*rho(i,j)*(1.+3.0*t2+4.50*t2*t2-1.50*t1)
      f(k,i,j)=omega*feq(k,i,j)+(1.-omega)*f(k,i,j)
      END DO
      END DO
      END DO
      return
      end
! ************************************************************************!
      subroutine streaming
      include 'param.h'
! streaming
      DO j=0,m
      DO i=n,1,-1           !right to left
      f(1,i,j)=f(1,i-1,j)
      END DO
      DO i=0,n-1            !left to right
      f(3,i,j)=f(3,i+1,j)
      END DO
      END DO
      DO j=m,1,-1           !bottom from top
      DO i=0,n
      f(2,i,j)=f(2,i,j-1)
      END DO
      DO i=n,1,-1
      f(5,i,j)=f(5,i-1,j-1)
      END DO
      DO i=0,n-1
      f(6,i,j)=f(6,i+1,j-1)
      END DO
      END DO
      DO j=0,m-1            !top from bottom
      DO i=0,n
      f(4,i,j)=f(4,i,j+1)
      END DO
      DO i=0,n-1
      f(7,i,j)=f(7,i+1,j+1)
      END DO
      DO i=n,1,-1
      f(8,i,j)=f(8,i-1,j+1)
      END DO
      END DO
      return
      end
! ************************************************************************!
      subroutine sfbound
      include 'param.h'
      do j=0,m
! inflow on west boundary, book-Eq.5.21-5.24
      rhow=f(0,0,j)+f(2,0,j)+f(4,0,j)+2.*(f(3,0,j)+f(6,0,j)+f(7,0,j))
      rhow=rhow/(1-uo)
      f(1,0,j)=f(3,0,j)+2.*rhow*uo/3.
      f(5,0,j)=f(7,0,j)+rhow*uo/6.
      f(8,0,j)=f(6,0,j)+rhow*uo/6.
      enddo
! bounce back on south boundary
      do i=0,n
      f(2,i,0)=f(4,i,0)
      f(5,i,0)=f(8,i,0)
      f(6,i,0)=f(7,i,0)
      end do
! bounce back, north boundary
      do i=0,n
      f(4,i,m)=f(2,i,m)
      f(8,i,m)=f(5,i,m)
      f(7,i,m)=f(6,i,m)
      end do
! account for open boundary condition at outlets
      do j=1,m
! Zou/He boundary at outlets
!      uxx=-1.+(f(0,n,j)+f(2,n,j)+f(4,n,j)+2.*(f(1,n,j)+f(5,n,j)
!     & +f(8,n,j)))/(rho(n,m/2))
!      f(3,i,j)=f(1,i,j)-4./6.*rho(n,m/2)*uxx
!      f(6,i,j)=f(8,i,j)+0.5*(f(2,i,j)-f(4,i,j))-1./6.*rho(n,m/2)*uxx
!      f(7,i,j)=f(5,i,j)-0.5*(f(2,i,j)-f(4,i,j))-1./6.*rho(n,m/2)*uxx
! isolated wall
      f(3,n,j)=f(3,n-1,j)
      f(7,n,j)=f(7,n-1,j)
      f(6,n,j)=f(6,n-1,j)
! open boundary condition
!      f(1,n,j)=2.*f(1,n-1,j)-f(1,n-2,j)
!      f(5,n,j)=2.*f(5,n-1,j)-f(5,n-2,j)
!      f(8,n,j)=2.*f(8,n-1,j)-f(8,n-2,j)
      enddo
! obstracle at the inlet x=90, y=130
!      do i=90,110
!      j=150
!      f(2,i,j)=f(4,i,j)
!      f(5,i,j)=f(7,i,j)
!      f(6,i,j)=f(8,i,j)
!      enddo
!      do j=130,150
!      i=110
!      f(1,i,j)=f(3,i,j)
!      f(5,i,j)=f(7,i,j)
!      f(8,i,j)=f(6,i,j)  !
!      enddo
!      do i=90,110
!      j=130
!      f(4,i,j)=f(2,i,j)
!      f(7,i,j)=f(5,i,j)
!      f(8,i,j)=f(6,i,j)
!      enddo
!      do j=130,150
!      i=90
!      f(3,i,j)=f(1,i,j)
!      f(6,i,j)=f(8,i,j)
!      f(7,i,j)=f(5,i,j)  !
!      enddo
! circle at 100x140
      do j=130,150
      x1=sqrt((100.)-(float(j)-140.)**2)+100.
      x2=-sqrt((100.)-(float(j)-140.)**2)+100.
      if(x2.lt.x1) then
          qx=abs(int(x2)-x2)
        if(qx.ge.0.5) then
!          do k=1,8
          i=int(x2)
!          op=opp(k)
!          f(op,i,j)=(1./2.*qx)*f(k,i,j)+(2.*qx-1.)/(2.*qx)*f(op,i,j)
          f(3,i,j)=(1./2.*qx)*f(1,i,j)+(2.*qx-1.)/(2.*qx)*f(3,i,j)
          f(6,i,j)=(1./2.*qx)*f(8,i,j)+(2.*qx-1.)/(2.*qx)*f(6,i,j)
          f(7,i,j)=(1./2.*qx)*f(5,i,j)+(2.*qx-1.)/(2.*qx)*f(7,i,j)
!          enddo
        else
!          do k=1,8
          i=int(x2)
!          op=opp(k)
!          f(op,i,j)=2.*qx*f(k,i,j)+(1.-2.*qx)*f(op,i-1,j)
          f(3,i,j)=2.*qx*f(1,i,j)+(1.-2.*qx)*f(3,i-1,j)
          f(6,i,j)=2.*qx*f(8,i,j)+(1.-2.*qx)*f(6,i-1,j+1)
          f(7,i,j)=2.*qx*f(5,i,j)+(1.-2.*qx)*f(7,i-1,j-1)
!          enddo
        endif
      else
          qx=abs(int(x1)-x1)
        if(qx.ge.0.5) then
!          do k=1,8
          i=int(x1)+1
!          op=opp(k)
!          f(op,i,j)=(1./2.*qx)*f(k,i,j)+(2.*qx-1.)/(2.*qx)*f(op,i,j)
          f(1,i,j)=(1./2.*qx)*f(3,i,j)+(2.*qx-1.)/(2.*qx)*f(1,i,j)
          f(5,i,j)=(1./2.*qx)*f(7,i,j)+(2.*qx-1.)/(2.*qx)*f(5,i,j)
          f(8,i,j)=(1./2.*qx)*f(6,i,j)+(2.*qx-1.)/(2.*qx)*f(8,i,j)
!          enddo
        else
!          do k=1,8
          i=int(x1)+1
!          op=opp(k)
!          f(op,i,j)=2.*qx*f(k,i,j)+(1.-2.*qx)*f(op,i+1,j)
          f(1,i,j)=2.*qx*f(3,i,j)+(1.-2.*qx)*f(1,i+1,j)
          f(5,i,j)=2.*qx*f(7,i,j)+(1.-2.*qx)*f(5,i+1,j+1)
          f(8,i,j)=2.*qx*f(6,i,j)+(1.-2.*qx)*f(8,i+1,j-1)
!          enddo
        endif
      endif
      enddo

      do i=90,110
      y1=sqrt((100.)-(float(i)-100.)**2)+140.
      y2=-sqrt((100.)-(float(i)-100.)**2)+140.
      if(y2.lt.y1) then
          qy=abs(int(y2)-y2)
        if(qy.ge.0.5) then
!          do k=1,8
          j=int(y2)
!          op=opp(k)
!          f(op,i,j)=(1./2.*qy)*f(k,i,j)+(2.*qy-1.)/(2.*qy)*f(op,i,j)
          f(4,i,j)=(1./2.*qy)*f(2,i,j)+(2.*qy-1.)/(2.*qy)*f(4,i,j)
          f(8,i,j)=(1./2.*qy)*f(6,i,j)+(2.*qy-1.)/(2.*qy)*f(8,i,j)
          f(7,i,j)=(1./2.*qy)*f(5,i,j)+(2.*qy-1.)/(2.*qy)*f(7,i,j)
!          enddo
        else
!          do k=1,8
          j=int(y2)
!          op=opp(k)
!          f(op,i,j)=2.*qy*f(k,i,j)+(1.-2.*qy)*f(op,i,j-1)
          f(4,i,j)=2.*qy*f(2,i,j)+(1.-2.*qy)*f(4,i,j-1)
          f(8,i,j)=2.*qy*f(6,i,j)+(1.-2.*qy)*f(8,i+1,j-1)
          f(7,i,j)=2.*qy*f(5,i,j)+(1.-2.*qy)*f(7,i-1,j-1)
!          enddo
        endif
      else
          qy=abs(int(y1)-y1)
        if(qy.ge.0.5) then
!          do k=1,8
          j=int(y1)+1
!          op=opp(k)
!          f(op,i,j)=(1./2.*qy)*f(k,i,j)+(2.*qy-1.)/(2.*qy)*f(op,i,j)
          f(2,i,j)=(1./2.*qy)*f(4,i,j)+(2.*qy-1.)/(2.*qy)*f(2,i,j)
          f(5,i,j)=(1./2.*qy)*f(7,i,j)+(2.*qy-1.)/(2.*qy)*f(5,i,j)
          f(6,i,j)=(1./2.*qy)*f(8,i,j)+(2.*qy-1.)/(2.*qy)*f(6,i,j)
!          enddo
        else
!          do k=1,8
          j=int(y1)+1
!          op=opp(k)
!          f(op,i,j)=2.*qy*f(k,i,j)+(1.-2.*qy)*f(op,i,j+1)
          f(2,i,j)=2.*qy*f(4,i,j)+(1.-2.*qy)*f(2,i,j+1)
          f(5,i,j)=2.*qy*f(7,i,j)+(1.-2.*qy)*f(5,i+1,j+1)
          f(6,i,j)=2.*qy*f(8,i,j)+(1.-2.*qy)*f(6,i-1,j+1)
!          enddo
        endif
      endif
      enddo	  

!      f(1,i,j)=f(1,i,j)
!      f(2,i,j)=f(4,i,j)
!      f(3,i,j)=f(1,i,j)
!      f(4,i,j)=f(2,i,j)
      f(5,110,150)=f(7,110,150)
      f(6,90,150)=f(8,90,150)
      f(7,90,130)=f(5,90,130)
      f(8,110,130)=f(6,110,130)
      return
      end
! *************************************************************************!
      real function opp(k)
      if (k.le.2) opp=k+2
      if (k.eq.3.or.k.eq.4) opp=k-2
      if (k.eq.5.or.k.eq.6) opp=k+2
      if (k.eq.7.or.k.eq.8) opp=k-2
      return
      end
! *************************************************************************!
      subroutine rhouv
      include 'param.h'
      do j=0,m
      do i=0,n
      ssum=0.0
      do k=0,8
      ssum=ssum+f(k,i,j)
      end do
      rho(i,j)=ssum 
      p(i,j)=rho(i,j)/3.
      end do
      end do
!      do i=1,n
!      rho(i,m)=2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
!      rho(i,m)=rho(i,m)+f(0,i,m)+f(1,i,m)+f(3,i,m)
!      end do
      DO i=1,n
      DO j=1,m-1
      usum=0
      vsum=0
      DO k=0,8
      usum=usum+f(k,i,j)*cx(k)
      vsum=vsum+f(k,i,j)*cy(k)
      END DO
      u(i,j)=usum/rho(i,j)
      v(i,j)=vsum/rho(i,j)
      END DO
      END DO
      do j=0,m
      v(n,j)=0.0
      enddo
      do j=130,150
      do i=90,110
      if ((real(i)-real(100))**2+(real(j)-real(140))**2.le.100.) then
!      rho(i,j)=0.0
      u(i,j)=0.0
      v(i,j)=0.0
      endif
      enddo
      enddo
      return
      end
! *************************************************************************
      subroutine vorti
      include 'param.h'
      real o(0:n,0:m)
! vorticity
      o=0.
      vor=0.
      do i=1,n-1
      do j=1,m-1
      o(i,j)=(v(i,j-1)-v(i-1,j-1))-(u(i-1,j)-u(i-1,j-1))
      enddo
      enddo
      do i=1,n-1
      do j=1,m-1
      vor(i,j)=(o(i,j)+o(i-1,j)+o(i,j-1)+o(i-1,j-1))/4
      enddo
      enddo
      return
      end
! *************************************************************************
      subroutine force
! force infront of bluff
      include 'param.h'
      sumr=0.
      sumv=0.
      do j=1,m-1
      i=10
      sumr=sumr+rho(i,j)
      sumv=sumv+u(i,j)
      enddo
      rh=sumr/(m-2)
      uu=sumv/(m-2)
      sumo=0.
      do j=130,150
      i=90
      sumi=0.
      do k=1,8
      if (k.eq.1) sumi=sumi+cx(k)*f(k,i,j)-cx(k)*f(3,i-1,j)
      if (k.eq.2) sumi=sumi+cx(k)*f(k,i,j)
      if (k.eq.3) sumi=sumi+cx(k)*f(k,i,j)
      if (k.eq.4) sumi=sumi+cx(k)*f(k,i,j)
      if (k.eq.5) sumi=sumi+cx(k)*f(k,i,j)-cx(k)*f(7,i-1,j-1)
      if (k.eq.6) sumi=sumi+cx(k)*f(k,i,j)
      if (k.eq.7) sumi=sumi+cx(k)*f(k,i,j)
      if (k.eq.8) sumi=sumi+cx(k)*f(k,i,j)-cx(k)*f(6,i-1,j+1)
      enddo
      sumo=sumo+sumi/7.
      enddo
      forx=sumo*2.
      cd=abs(forx)/(rh*uo**2*20.)
      sumo=0.
      do j=130,150
      i=90
      sumi=0.
      do k=1,8
      if (k.eq.1) sumi=sumi+cy(k)*f(k,i,j)-cy(k)*f(3,i-1,j)
      if (k.eq.2) sumi=sumi+cy(k)*f(k,i,j)
      if (k.eq.3) sumi=sumi+cy(k)*f(k,i,j)
      if (k.eq.4) sumi=sumi+cy(k)*f(k,i,j)
      if (k.eq.5) sumi=sumi+cy(k)*f(k,i,j)-cy(k)*f(7,i-1,j-1)
      if (k.eq.6) sumi=sumi+cy(k)*f(k,i,j)
      if (k.eq.7) sumi=sumi+cy(k)*f(k,i,j)
      if (k.eq.8) sumi=sumi+cy(k)*f(k,i,j)-cy(k)*f(6,i-1,j+1)
      enddo
      sumo=sumo+sumi/7.
      enddo
      fory=sumo*2.
      cl=fory/(rh*uo**2*20.)
      cp=2*(p(90,140)-rh/3.)/rh/uu**2
      print*, 'cd', cd, 'cl', cl, 'cp', cp
      do i=90,110
      cpp(i,130)=2*(p(i,130)-rh/3.)/rh/uu**2
      cpp(i,150)=2*(p(i,150)-rh/3.)/rh/uu**2
      enddo
      do j=130,150
      cpp(90,j)=2*(p(90,j)-rh/3.)/rh/uu**2
      cpp(110,j)=2*(p(110,j)-rh/3.)/rh/uu**2
      enddo
      return
      end
! *************************************************************************
      subroutine resul
      include 'param.h'
      real strf(0:n,0:m)
      character*6 rr
	  
! stream function calculations
      strf(0,0)=0.
      do i=0,n
      rhoav=0.5*(rho(i-1,0)+rho(i,0))
      if(i.ne.0) strf(i,0)=strf(i-1,0)-rhoav*0.5*(v(i-1,0)+v(i,0))
      do j=1,m
      rhom=0.5*(rho(i,j)+rho(i,j-1))
      strf(i,j)=strf(i,j-1)+rhom*0.5*(u(i,j-1)+u(i,j))
      end do
      end do
	  
      write(rr,300) 100000+kk/100 
300   format(i6)
      open(2,file='con'//TRIM(rr)//'.dat')
	  
      write(2,*)' VARIABLES =X, Y, U, V, T, O'
      write(2,*)'ZONE ','I=',n+1,'J=',m+1,',','F=BLOCK'
      do j=0,m
      write(2,*)(i,i=0,n)
      end do
      do j=0,m
      write(2,*)(j,i=0,n)
      end do
      do j=0,m
      write(2,*)(u(i,j),i=0,n)
      end do
      do j=0,m
      write(2,*)(v(i,j),i=0,n)
      end do
      do j=0,m
      write(2,*)(tau(i,j),i=0,n)
      end do
      do j=0,m
      write(2,*)(vor(i,j),i=0,n)
      end do
      do i=0,n
      write(3,*) i/float(20)-5,u(i,m/2)/uo
      end do
      do i=0,n
      write(4,*) i/float(20)-5,v(i,m/2)
      end do
! streamcross velocity at center of obstacle
      do j=0,m
      write(19,*) j/float(20),u(100,j)/uo,u(170,j)/uo,u(250,j)/uo
      end do
      do j=0,m
      write(20,*) j/float(20),v(100,j),v(170,j),v(250,j)
      end do
! printing cpp coefficient
      open(10,file='cpp.dat')
      do j=140,150
      write(10,*) j-140, cpp(90,j)
      enddo
      do i=90,110
      write(10,*) i-80, cpp(i,150)
      enddo
      do j=150,130,-1
      write(10,*) 180-j, cpp(110,j)
      enddo
      do i=110,90,-1
      write(10,*) 160-i, cpp(i,130)
      enddo
      do j=130,140
      write(10,*) j-60, cpp(90,j)
      enddo
      close (10)
      return
      end
!============end of the program
