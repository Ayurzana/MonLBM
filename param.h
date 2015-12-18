      parameter (n=840, m=280, mstep=10000)
      parameter (rhoo=2.0, uo=0.11 )
      parameter (ome=1.85)
      parameter (dx=1., dy=1., dt=1.     )

      real f(0:8,0:n,0:m), feq(0:8,0:n,0:m)
      real rho(0:n,0:m), vor(0:n,0:m)
      real w(0:8), cx(0:8), cy(0:8)
      real u(0:n,0:m), v(0:n,0:m), tau(0:n,0:m)
      real p(0:n,0:m), cd, cp, cl, cpp(90:110,130:150)
      integer kk
	  
      common /data/f,rho,vor,w,cx,cy,u,v,tau,p,cd,cp,cl,cpp
      common /data/ kk
	  

! Re=465 uo=0.12, ome=1.94
! Re=477   uo=0.1   , ome=1.951
! Re=514   uo=0.11   , ome=1.95
! Re=536   uo=0.11   , ome=1.952  
! Re=560   uo=0.1115   , ome=1.954 
