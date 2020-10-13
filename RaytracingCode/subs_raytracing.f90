MODULE subs_raytracing
	
	use constantes
	use environ
	implicit none
	
	contains

subroutine varplasma(V,f,X,Y)
	!------------------------------------
	! Calcul les variables X et Y du plasma
	! CALLING SEQUENCE: call varplasma(V,f)
	! INPUTS : V: vect position en coord. cart.
	!		   f: frequence en Hz
	! OUTPUTS: X: (fp/f)^2
	!	  	   Y: fc/f
	!------------------------------------
	real(kind=8), dimension(:,:), intent(in)  :: V 
	real(kind=8), intent(in)                  :: f 
	real(kind=8), dimension(:), intent(out)   :: X 
	real(kind=8), dimension(:,:), intent(out) :: Y 
	real(kind=8), dimension(3,size(x))        :: B 
	real(kind=8), dimension(size(X))          :: Ne, fp
		
	call density(V,Ne)
	fp = 1.d0/(2.d0*pi)*sqrt(e**2*Ne/(eps0*me))*1.d3 !Hz
	X = (fp/f)**2
	
	call magneticf(V,B)
	Y = e*B/(me*2.d0*pi*f)
end subroutine varplasma

subroutine index(X,Y,k,mode,mu)
	!------------------------------------
	! calcul l'indice de refraction mu dans un plasma magnetise
	! d'apres la formule d'Appleton-Hartree donnee dans
	!    [Haselgrove, 1960]
	! CALLING SEQUENCE: call index(X,Y,th,mu)
	! INPUTS: X: (fp/f)^2
	!		  Y: fc/f
	!		 th: angle entre k et B
	! OUTPUTS: mu
	!------------------------------------
	real(kind=8), dimension(:), intent(in)     	 :: X
	real(kind=8), dimension(:,:), intent(in)   	 :: Y,k
	character(len=1), intent(in), dimension(:) 	 :: mode
	real(kind=8), dimension(:), intent(out)    	 :: mu
	real(kind=8), dimension(size(mu))          	 :: normY,normk,vdotk,A,B,C,B2mAC
	real(kind=8), dimension(size(mu))		   	 :: nmode
	

	where (mode .eq. 'X')
		nmode =-1.d0
	else where
		nmode = 1.d0
	end where
	
	normY = sqrt(Y(1,:)**2+Y(2,:)**2+Y(3,:)**2)
	normk = sqrt(k(1,:)**2+k(2,:)**2+k(3,:)**2)
	where (normY .eq.0.d0)
		vdotk(:) = maxval(k(:,:))
	elsewhere
		vdotk(:) = (Y(1,:)*k(1,:)+Y(2,:)*k(2,:)+Y(3,:)*k(3,:))/normY(:)
	end where
	
	A = 1.d0-X-normY**2*(1.d0-X*vdotk**2/normk**2)
	B = 0.5d0*(-2.d0*(1.d0-X)*(1.d0-X-normY**2)+X*normY**2*(1.d0-vdotk**2/normk**2))
	C = (1.d0-X)*((1.d0-X)**2-normY**2)
	
	B2mAC = B**2-A*C
	
	!Qd champ magnŽtqiue -> 0 alors B2mAC doit -> 0 => ARTEFACTS NUMERIQUES => qd abs(B2mAC) < 1.d-15 => B2mAC=0
	where (abs(B2mAC) < 1.d-15)
		B2mAC = 0.d0
	end where
	
	mu = sqrt((-B+nmode*sqrt(B2mAC))/A)
	
	
	!	write(*,*) B2mAC, -B+sqrt(B2mAC), mu
	
end subroutine index

subroutine calcp(X,Y,k,mode,p)
	!------------------------------------
	! calcul la variable p solution de
	!	a*p^2+2b*p+c*=0
	! d'Haselgrove (1963) 
	! CALLING SEQUENCE: call calcp(X,Y,k,p)
	! INPUTS: X: (fp/f)^2
	!		  Y: fc/f
	!		  k: vecteur d'onde en coord. cart.
	! OUTPUTS: p [Haselgrove,63]
	!------------------------------------
	real(kind=8), dimension(:,:), intent(in)   :: Y,k 
	real(kind=8), dimension(:), intent(in)     :: X
	character(len=1), intent(in), dimension(:) :: mode
	real(kind=8), dimension(:), intent(out)    :: p
	real(kind=8), dimension(size(p))           :: normY, vdotk,alpha,beta,gamma,nmode
	
	!------------------------------------------
	!	mode O +/- mode X
	!p = (-beta +/- sqrt(beta**2-alpha*gamma))/alpha

	
	where (mode .eq. 'X')
		nmode =-1.d0
	else where
		nmode = 1.d0
	end where
	
	normY = sqrt(Y(1,:)**2+Y(2,:)**2+Y(3,:)**2)
	
	!si B=0 (Y=0) alors direction de Y = verticale
	where (normY .eq. 0.d0) 
		vdotk(:) = maxval(k(:,:))
	elsewhere
		vdotk(:) = (Y(1,:)*k(1,:)+Y(2,:)*k(2,:)+Y(3,:)*k(3,:))/normY(:)
	end where
	
	!determination de p: pY=(q+1) et q=(n^2-1)/X
	alpha = 1.d0-X-normY**2
	beta = normY*(1.+(vdotk)**2)/2.
	gamma = -(vdotk)**2
	
	p = (-beta+nmode*sqrt(beta**2-alpha*gamma))/alpha
	!p = gamma/(-beta -nmode*sqrt(beta**2-alpha*gamma))

end subroutine calcp

subroutine varpartiel(X,Y,k,mode,A,B,C,D)
	!------------------------------------
	! calcul les variables J*,K*,L* et M* 
	! d'Haselgrove (1963) decrivant les d_partielles
	! CALLING SEQUENCE: call varpartiel(X,Y,k,A,B,C,D)
	! INPUTS: X: (fp/f)^2
	!		  Y: fc/f
	!		  k: vecteur d'onde en coord. cart.
	! OUTPUTS: A == J* [Haselgrove,63]
	!          B == K* [Haselgrove,63]
	!          C == L* [Haselgrove,63]
	!          D == M* [Haselgrove,63]
	!------------------------------------
	real(kind=8), dimension(:,:), intent(in)  :: Y,k 
	real(kind=8), dimension(:), intent(in)    :: X 
	character(len=1), intent(in),dimension(:) :: mode
	real(kind=8), dimension(:), intent(out)   :: A,B,C,D
	real(kind=8),dimension(size(X))           :: normY, vdotk,p,nmode
	
	normY = sqrt(Y(1,:)**2+Y(2,:)**2+Y(3,:)**2)

	where (mode .eq. 'X')
		nmode =-1.d0
	else where
		nmode = 1.d0
	end where
		!si B=0 (Y=0) alors direction de Y = verticale (vdotk ne pas tre nul)
	where (normY .eq. 0.d0)  
		vdotk(:) = maxval(k(:,:))
	elsewhere
		vdotk(:) = (Y(1,:)*k(1,:)+Y(2,:)*k(2,:)+Y(3,:)*k(3,:))/normY(:)
	end where

	call calcp(X,Y,k,mode,p)

	A = nmode*4.d0*sqrt((normY*(1.d0+(vdotk)**2)/2.)**2-(1.d0-X-normY**2)*(-(vdotk)**2))
	B = -2.d0*X*vdotk*(p*normY-1.d0)
	C = (normY*(1.d0-normY**2)*p**2-2.d0*(1.-X-normY**2)*p-normY)
	D = 2.d0*X*p*(p*normY-1.d0)
	
	!write(*,*) '-----------------------------------------'
	!write(*,*)'X', X
	!write(*,*)'Y_x', Y(1,:)
	!write(*,*)'Y_y', Y(2,:)
	!write(*,*)'Y_z', Y(3,:)
	!write(*,*)'normY', normY
	!write(*,*) 'k', k
	!write(*,*) 'vdotk', vdotk
	!write(*,*) '-----------------------------------------'
	!write(*,*) 'p', p
	!write(*,*) '-----------------------------------------'
	!write(*,*) 'J',A
	!write(*,*) 'K',B
	!write(*,*) 'L',C
	!write(*,*) 'M',D
		
end subroutine varpartiel

subroutine dxdt(V,k,f,mode,dx_dt,nray)
	!------------------------------------
	! calcul de dxidt = J*ki-K*Yi
	! 	d'aprs Haselgrove (1960)
	! CALLING SEQUENCE: call dxdt(V,k,f,dx_dt)
	! INPUTS : V: vect position en coord. cart.
	!          k: vect d'onde en coord. cart.
	!		   f: frequence en Hz
	! OUTPUTS: dx_dt : tableau de 3*Nray dim.
	!------------------------------------
	integer, intent(in)							 :: nray
	real(kind=8), dimension(3,nray), intent(in)  :: V,k
	real(kind=8), intent(in)                     :: f
	character(len=1), intent(in), dimension(nray)   :: mode	
	real(kind=8), dimension(3,nray), intent(out) :: dx_dt
	real(kind=8), dimension(3,nray)              :: Y
	real(kind=8), dimension(nray)                :: X, A,B,C,D
	integer									     :: i
	
	call varplasma(V,f,X,Y)
	call varpartiel(X,Y,k,mode,A,B,C,D)

	forall (i=1:3) dx_dt(i,:)=A(:)*k(i,:)-B(:)*Y(i,:)
	!write(*,*) '---------------------------'
	!write(*,*) 'dxdt_x', dx_dt(1,:)
	!write(*,*) 'dxdt_y', dx_dt(2,:)
	!write(*,*) 'dxdt_z', dx_dt(3,:)

end subroutine dxdt

subroutine dkdt(V,k,f,mode,dk_dt,nray)
	!------------------------------------
	! calcul de dkidt = L*dX/dxi+sumj((K*kj+M*Yj)*dYj/dxi)
	! 	d'aprs Haselgrove (1960)
	! CALLING SEQUENCE: call dkdt(V,k,f,dk_dt)
	! INPUTS : V: vect position en coord. cart.
	!          k: vect d'onde en coord. cart.
	!		   f: frequence en Hz
	! OUTPUTS: dk_dt: tableau de 3 dim. 
	!------------------------------------
	integer, intent(in)							 :: nray
	real(kind=8), dimension(3,nray), intent(in)  :: V,k
	real(kind=8), intent(in)                     :: f
	character(len=1), intent(in),dimension(nray)    :: mode
	real(kind=8), dimension(3,nray), intent(out) :: dk_dt
	real(kind=8), dimension(3, nray)             :: Vav, Vap, Yav,Yap,Y,dYdx
	real(kind=8), dimension(nray)                :: h, Xav, Xap,A,B,C,D,X,dXdx
	integer                                      :: i,j

	call varplasma(V,f,X,Y)
	call varpartiel(X,Y,k,mode,A,B,C,D)	
		
	h = 1.d-6
	do i = 1,3
		Vav = V
		Vav(i,:) = Vav(i,:)-h
		call varplasma(Vav,f,Xav,Yav)
		Vap = V
		Vap(i,:) = Vap(i,:)+h
		call varplasma(Vap,f,Xap,Yap)
		dXdx = (Xap-Xav)/(2.d0*h)
		forall (j=1:3) dYdx(j,:) = (Yap(j,:)-Yav(j,:))/(2.d0*h)
		dk_dt(i,:) = C*dXdx+((B*k(1,:)+D*Y(1,:))*dYdx(1,:)+(B*k(2,:)+D*Y(2,:))*dYdx(2,:)+(B*k(3,:)+D*Y(3,:))*dYdx(3,:))
	end do
	
	!write(*,*) '---------------------------'
	!write(*,*) 'dudt_x', dk_dt(1,:)
	!write(*,*) 'dudt_y', dk_dt(2,:)
	!write(*,*) 'dudt_z', dk_dt(3,:)
end subroutine dkdt

subroutine error(k,X,Y,mode,err)
	!---------------------------------
	! calcul de l'erreur k^2-1-X*(pY-1)
	! 	d'aprs Haselgrove (1960)
	! CALLING SEQUENCE: call error(k,X,Y, err)
	! INPUTS : V: vect position en coord. cart.
	!          k: vect d'onde en coord. cart.
	!		   f: frequence en Hz
	! OUTPUTS: err: scalaire = 0
	!------------------------------------
	real(kind=8), dimension(:,:), intent(in)  :: Y,k
	real(kind=8), dimension(:), intent(in)    :: X
	character(len=1), intent(in),dimension(:) :: mode
	real(kind=8), dimension(:), intent(out)   :: err
	real(kind=8), dimension(size(X))          :: p, normk2,normY
	
	normk2=k(1,:)**2+k(2,:)**2+k(3,:)**2
	normY = sqrt(Y(1,:)**2+Y(2,:)**2+Y(3,:)**2)
	call calcp(X,Y,k,mode,p)
		
	err = abs((normk2-1.d0-X*(p*normY-1.d0))/(X*(p*normY-1.d0)))

end subroutine error


END MODULE subs_raytracing
