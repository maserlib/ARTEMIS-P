MODULE RungeKutta
	
	use subs_raytracing
	implicit none

	contains

	subroutine RK4(x,y,dt,dxdt,dydt,param1,param2,x_1,y_1)
		!-----------------------------------------------------------------
		! Methode de Runge Kutta d'ordre 4 coupl�e
		! CALLING SEQUENCE : RK4(x,y,dt,dxdt,dydt,param1,x_1,y_1)
		! INPUTS:     x: vecteur de dim 3*Nray
		!             y: vecteur de dim 3*Nray
		!            dt: pas d'integration
		!          dxdt: fonction donnant la derivee de x
		!          dydt: fonction donnant la derivee de y
		!        param1: param�tre --> f
		!		 param2: param�tre --> mode ('X' ou 'O' ou 'V')
		! OUTPUTS:  x_1: vecteur de dim 3*Nray au pas n+1
		!           y_1: vecteur de dim 3*Nray au pas n+1
		!-----------------------------------------------------------------
	implicit none
	real(kind=8), dimension(:,:), intent(in)  :: x,y
	real(kind=8), dimension(:), intent(in)    :: dt
	real(kind=8), intent(in)                  :: param1
	character(len=1), intent(in), dimension(:):: param2
	real(kind=8), dimension(:,:), intent(out) :: x_1,y_1
	 
	external                                  :: dxdt, dydt
	
	real(kind=8), dimension(3,size(dt))       :: dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4
	real(kind=8), dimension(3,size(dt))       :: dx1dt,dy1dt,dx2dt,dy2dt,dx3dt	
	real(kind=8), dimension(3,size(dt))       :: dy3dt,dx4dt,dy4dt, Yp
	real(kind=8), dimension(size(dt))		  :: Xp
	integer									  :: i
	real(kind=8), dimension(size(dt))		  :: nmode
	
	call dxdt(x,y,param1,param2,dx1dt,size(dt)) !drdt
	call varplasma(x,param1,Xp,Yp)	
	forall (i=1:3) dx1(i,:)=dt(:)*dx1dt(i,:)!*sqrt((1.d0-Xp)/(dx1dt(1,:)**2+dx1dt(2,:)**2+dx1dt(3,:)**2))
	
	call dydt(x,y,param1,param2,dy1dt,size(dt)) !dkdt
	forall (i=1:3) dy1(i,:)=dt(:)*dy1dt(i,:)!*sqrt((1.d0-Xp)/(dx1dt(1,:)**2+dx1dt(2,:)**2+dx1dt(3,:)**2))

	call dxdt(x+dx1/2.d0,y+1.d0/2.d0*dy1,param1,param2,dx2dt,size(dt)) !drdt
	call varplasma(x+dx1/2.d0,param1,Xp,Yp)
	forall (i=1:3) dx2(i,:)=dt(:)*dx2dt(i,:)!*sqrt((1.d0-Xp)/(dx2dt(1,:)**2+dx2dt(2,:)**2+dx2dt(3,:)**2))
		
	call dydt(x+dx1/2.d0,y+1.d0/2.d0*dy1,param1,param2,dy2dt,size(dt)) !dkdt
	forall (i=1:3) dy2(i,:)=dt(:)*dy2dt(i,:)!*sqrt((1.d0-Xp)/(dx2dt(1,:)**2+dx2dt(2,:)**2+dx2dt(3,:)**2))

	call dxdt(x+dx2/2.d0,y+1.d0/2.d0*dy2,param1,param2,dx3dt,size(dt)) !drdt
	call varplasma(x+dx2/2.d0,param1,Xp,Yp)
	forall (i=1:3) dx3(i,:)=dt(:)*dx3dt(i,:)!*sqrt((1.d0-Xp)/(dx3dt(1,:)**2+dx3dt(2,:)**2+dx3dt(3,:)**2))

	call dydt(x+dx2/2.d0,y+1.d0/2.d0*dy2,param1,param2,dy3dt,size(dt)) !drdt
	forall (i=1:3) dy3(i,:)=dt(:)*dy3dt(i,:)!*sqrt((1.d0-Xp)/(dx3dt(1,:)**2+dx3dt(2,:)**2+dx3dt(3,:)**2))

	call dxdt(x+dx3,y+dy3,param1,param2,dx4dt,size(dt)) !drdt
	call varplasma(x+dx3,param1,Xp,Yp)
	forall (i=1:3) dx4(i,:)=dt(:)*dx4dt(i,:)!*sqrt((1.d0-Xp)/(dx4dt(1,:)**2+dx4dt(2,:)**2+dx4dt(3,:)**2))

	call dydt(x+dx3,y+dy3,param1,param2,dy4dt,size(dt)) !drdt
	forall (i=1:3) dy4(i,:)=dt(:)*dy4dt(i,:)!*sqrt((1.d0-Xp)/(dx4dt(1,:)**2+dx4dt(2,:)**2+dx4dt(3,:)**2))

	!changement de signe des d�riv�es dans Haselgrove_63 n�cessaire pour avoir les bonnes �quations sans B
	where (param2 .eq. 'X')
		nmode =-1.d0
	else where
		nmode = 1.d0
	end where
	

	forall(i=1:3) x_1(i,:) = x(i,:) + nmode(:)*1.d0/6.d0*(dx1(i,:)+2.d0*dx2(i,:)+2.d0*dx3(i,:)+dx4(i,:))
	forall(i=1:3) y_1(i,:) = y(i,:) + nmode(:)*1.d0/6.d0*(dy1(i,:)+2.d0*dy2(i,:)+2.d0*dy3(i,:)+dy4(i,:))

	end subroutine RK4
END MODULE RungeKutta
