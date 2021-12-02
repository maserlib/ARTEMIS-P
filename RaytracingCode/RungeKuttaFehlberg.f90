MODULE RungeKuttaFehlberg
	use subs_raytracing
	implicit none
	
	contains
	subroutine RK4F(x,y,dt,dxdt,dydt,param1,param2,dtmin,dtmax,x_1,y_1)
		!-----------------------------------------------------------------
		! Runge–Kutta-Fehlberg fourth-order method
		! (Coleman 2008 radio science bulletin)
		! CALLING SEQUENCE : RK4F(x,y,dt,dxdt,dydt,param1,param2,x_1,y_1)
		! INPUTS:     x: 3*Nray dimension vector
		!             y: 3*Nray dimension vector
		!            dt: integration time step
		!          dxdt: the derivative of x
		!          dydt: the derivative of y
		!        param1: parameter --> f
		!	 param2: parameter --> mode ('X' ou 'O' ou 'V')
		! OUTPUTS:  x_1: 3*Nray dimension vector at n+1 step
		!           y_1: 3*Nray dimension vector at n+1 step
		!-----------------------------------------------------------------
	implicit none
	real(kind=8), dimension(:,:), intent(in)  :: x,y
	real(kind=8), dimension(:), intent(inout) :: dt
	real(kind=8), intent(in)                  :: param1, dtmin,dtmax
	character(len=1), intent(in),dimension(:) :: param2
	real(kind=8), dimension(:,:), intent(out) :: x_1,y_1
	 
	external                                  :: dxdt, dydt
	
	real(kind=8), dimension(3,size(dt))       :: dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4,dx5,dy5,dx6,dy6,Yp, dx,dy
	real(kind=8), dimension(3,size(dt))       :: dx1dt,dy1dt,dx2dt,dy2dt,dx3dt,dy3dt,dx4dt,dy4dt,dx5dt,dy5dt,dx6dt,dy6dt
	real(kind=8), dimension(size(dt))		  :: Xp
	integer									  :: i, again
	real(kind=8), dimension(size(dt))		  :: nmode, err1x,err2x,err3x, err1y, err2y, err3y
	real(kind=8)							  :: tol = 1.d-5, a1, a2, a3, a4, a5
	
	again = 1
	
	do while (again == 1)
	
	call dxdt(x,y,param1,param2,dx1dt,size(dt)) !drdt
	call varplasma(x,param1,Xp,Yp)
	forall (i=1:3) dx1(i,:)=dt*dx1dt(i,:)*sqrt((1.d0-Xp)/(dx1dt(1,:)**2+dx1dt(2,:)**2+dx1dt(3,:)**2))
	
	call dydt(x,y,param1,param2,dy1dt,size(dt)) !dkdt
	forall (i=1:3) dy1(i,:)=dt*dy1dt(i,:)*sqrt((1.d0-Xp)/(dx1dt(1,:)**2+dx1dt(2,:)**2+dx1dt(3,:)**2))
	
	call dxdt(x+dx1/4.d0,y+dy1/4.d0,param1,param2,dx2dt,size(dt)) !drdt
	call varplasma(x+dx1/4.,param1,Xp,Yp)
	forall (i=1:3) dx2(i,:)=dt*dx2dt(i,:)*sqrt((1.d0-Xp)/(dx2dt(1,:)**2+dx2dt(2,:)**2+dx2dt(3,:)**2))
	
	call dydt(x+dx1/4.d0,y+dy1/4.d0,param1,param2,dy2dt,size(dt)) !dkdt
	forall (i=1:3) dy2(i,:)=dt*dy2dt(i,:)*sqrt((1.d0-Xp)/(dx2dt(1,:)**2+dx2dt(2,:)**2+dx2dt(3,:)**2))
	
	call dxdt(x+3.d0*dx1/32.d0+9.*dx2/32.d0,y+3.d0*dy1/32.d0+9.d0*dy2/32.d0,param1,param2,dx3dt,size(dt)) !drdt
	call varplasma(x+3.d0*dx1/32.d0+9.d0*dx2/32.d0,param1,Xp,Yp)
	forall (i=1:3) dx3(i,:)=dt*dx3dt(i,:)*sqrt((1.d0-Xp)/(dx3dt(1,:)**2+dx3dt(2,:)**2+dx3dt(3,:)**2))
	
	call dydt(x+3.d0*dx1/32.d0+9.d0*dx2/32.d0,y+3.d0*dy1/32.d0+9.d0*dy2/32.d0,param1,param2,dy3dt,size(dt)) !drdt
	forall (i=1:3) dy3(i,:)=dt*dy3dt(i,:)*sqrt((1.d0-Xp)/(dx3dt(1,:)**2+dx3dt(2,:)**2+dx3dt(3,:)**2))
	
	call dxdt(x+1932.d0*dx1/2197.d0-7200.d0*dx2/2197.d0+7296.d0*dx3/2197.d0,&
		& y+1932.d0*dy1/2197.d0-7200.d0*dy2/2197.d0+7296.d0*dy3/2197.d0,&
		& param1,param2,dx4dt,size(dt)) !drdt
	call varplasma(x+1932.d0*dx1/2197.d0-7200.d0*dx2/2197.d0+7296.d0*dx3/2197.d0,param1,Xp,Yp)
	forall (i=1:3) dx4(i,:)=dt*dx4dt(i,:)*sqrt((1.d0-Xp)/(dx4dt(1,:)**2+dx4dt(2,:)**2+dx4dt(3,:)**2))
	
	call dydt(x+1932.d0*dx1/2197.d0-7200.d0*dx2/2197.d0+7296.d0*dx3/2197.d0,&
		& y+1932.d0*dy1/2197.d0-7200.d0*dy2/2197.d0+7296.d0*dy3/2197.d0,&
		& param1,param2,dy4dt,size(dt)) !drdt
	forall (i=1:3) dy4(i,:)=dt*dy4dt(i,:)*sqrt((1.d0-Xp)/(dx4dt(1,:)**2+dx4dt(2,:)**2+dx4dt(3,:)**2))
	
	call dxdt(x+439.d0*dx1/216.d0-8.d0*dx2+3680.d0*dx3/513.d0-845.d0*dx4/4104.d0,&
		& y+439.d0*dy1/216.d0-8.d0*dy2+3680.d0*dy3/513.d0-845.d0*dy4/4104.d0,&
		& param1,param2,dx5dt,size(dt)) !drdt
	call varplasma(x+439.d0*dx1/216.d0-8.d0*dx2+3680.d0*dx3/513.d0-845.d0*dx4/4104.d0,param1,Xp,Yp)
	forall (i=1:3) dx5(i,:)=dt*dx5dt(i,:)*sqrt((1.d0-Xp)/(dx5dt(1,:)**2+dx5dt(2,:)**2+dx5dt(3,:)**2))
	
	call dydt(x+439.d0*dx1/216.d0-8.d0*dx2+3680.d0*dx3/513.d0-845.d0*dx4/4104.d0,&
		& y+439.d0*dy1/216.d0-8.d0*dy2+3680.d0*dy3/513.d0-845.d0*dy4/4104.d0,&
		& param1,param2,dy5dt,size(dt)) !drdt
	forall (i=1:3) dy5(i,:)=dt*dy5dt(i,:)*sqrt((1.d0-Xp)/(dx5dt(1,:)**2+dx5dt(2,:)**2+dx5dt(3,:)**2))
	
	call dxdt(x-8.d0*dx1/27.d0+2.d0*dx2-3544.d0*dx3/2565.d0+1859.d0*dx4/4104.d0-11.d0*dx5/40.d0, &
		& y-8.d0*dy1/27.d0+2.d0*dy2-3544.d0*dy3/2565.d0+1859.d0*dy4/4104.d0-11.d0*dy5/40.d0,param1,param2,dx6dt,size(dt)) !drdt
	call varplasma(x-8.d0*dx1/27.d0+2.d0*dx2-3544.d0*dx3/2565.d0+1859.d0*dx4/4104.d0-11.d0*dx5/40.d0,param1,Xp,Yp)
	forall (i=1:3) dx6(i,:)=dt*dx6dt(i,:)*sqrt((1.d0-Xp)/(dx6dt(1,:)**2+dx6dt(2,:)**2+dx6dt(3,:)**2))
	
	call dydt(x-8.d0*dx1/27.d0+2.d0*dx2-3544.d0*dx3/2565.d0+1859.d0*dx4/4104.d0-11.d0*dx5/40.d0,&
		& y-8.d0*dy1/27.d0+2.d0*dy2-3544.d0*dy3/2565.d0+1859.d0*dy4/4104.d0-11.d0*dy5/40.d0,param1,param2,dy6dt,size(dt)) !drdt
	forall (i=1:3) dy6(i,:)=dt*dy6dt(i,:)*sqrt((1.d0-Xp)/(dx6dt(1,:)**2+dx6dt(2,:)**2+dx6dt(3,:)**2))
	
	!step adaptation
	a1 = (16./135-25./216.) ; a2 = (6656./12825.-1408./2565.) ; a3 = (28561./56430.-2197./4104.); a4 = (-9./50.+1./5.) ; a5 = 2./55.
	dx = a1*dx1+a2*dx3+a3*dx4+a4*dx5+a5*dx6
	dy = a1*dy1+a2*dy3+a3*dy4+a4*dy5+a5*dy6
	err1x = abs(dx(1,:)) ; err2x = abs(dx(2,:)) ; err3x = abs(dx(3,:))
	err1y = abs(dy(1,:)) ; err2y = abs(dy(2,:)) ; err3y = abs(dy(3,:))

	if (maxval([err1x,err2x,err3x,err1y,err2y,err3y]) > tol .and. minval(dt) > dtmin) then
		again = 1
		write(*,*)'------------------------------------'
		write(*,*) dt
		write(*,*)'-------'
		where ((err1x > tol) .and. dt > dtmin)
			dt = 0.84*dt*(tol*dt/maxval(err1x))**(1./4.)
		end where
		where ((err2x > tol) .and. dt > dtmin)
			dt = 0.84*dt*(tol*dt/maxval(err2x))**(1./4.)
		end where
		where ((err3x > tol) .and. dt > dtmin)
			dt = 0.84*dt*(tol*dt/maxval(err3x))**(1./4.)
		end where
		where ((err1y > tol) .and. dt > dtmin)
			dt = 0.84*dt*(tol*dt/maxval(err1y))**(1./4.)
		end where
		where ((err2y > tol) .and. dt > dtmin)
			dt = 0.84*dt*(tol*dt/maxval(err2y))**(1./4.)
		end where
		where ((err3y > tol) .and. dt > dtmin)
			dt = 0.84*dt*(tol*dt/maxval(err3y))**(1./4.)
		end where
		
		
		write(*,*)'------------------------------------'
		write(*,*) dt
		write(*,*)'------------------------------------'

	else
		again = 0
	endif		
	end do


	
	!change of sign of derivatives in Haselgrove_63 needed to get the right equations without B
	where (param2 .eq. 'X')
		nmode =-1.d0
	else where
		nmode = 1.d0
	end where
	
	forall (i=1:3) x_1(i,:) = x(i,:) + nmode(:)*(16.d0*dx1(i,:)/135.d0+6656.d0*dx3(i,:)/12825.d0+28561.d0*dx4(i,:)/56430.d0-&
	&9.d0*dx5(i,:)/50.d0+2.d0*dx6(i,:)/55.d0)
	forall (i=1:3) y_1(i,:) = y(i,:) + nmode(:)*(16.d0*dy1(i,:)/135.d0+6656.d0*dy3(i,:)/12825.d0+28561.d0*dy4(i,:)/56430.d0-&
	&9.d0*dy5(i,:)/50.d0+2.d0*dy6(i,:)/55.d0)
	end subroutine RK4F
END MODULE RungeKuttaFehlberg
