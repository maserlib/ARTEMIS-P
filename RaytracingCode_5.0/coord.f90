MODULE coord
	
	use constantes
	implicit none
	
	contains

subroutine cart_sph(x,r)
	!----------------------------------------------------------
	! Passage des coord. cartésiennes aux coord. sphériques 
	! CALLING SEQUENCE: call cart_sph(x,r)
	! INPUT:  x: vecteur position en coordonnées cart.
	! OUTPUT: r: vecteur position en coordonnées sph. 
	!				(rayon,longitude,latitude)			!
	!----------------------------------------------------------

	real(kind=8), dimension(:,:), intent(in)  :: x
	real(kind=8), dimension(:,:), intent(out) :: r

	r(1,:) = sqrt(x(1,:)**2+x(2,:)**2+x(3,:)**2)
	
	!longitude
	where (x(2,:) < 0.d0)
		r(2,:) = 2.d0*pi-acos(x(1,:)/sqrt(x(1,:)**2+x(2,:)**2)) 
	elsewhere  
		r(2,:) = acos(x(1,:)/sqrt(x(1,:)**2+x(2,:)**2)) 
	end where

	where (x(1,:) .eq. 0.d0 .and. x(2,:).eq. 0.d0)
		r(2,:) = pi/2.d0
	end where
	
	!latitude
	where (r(1,:) .ne. 0.d0)
		r(3,:) = pi/2.d0 - acos(x(3,:)/r(1,:))
	elsewhere
		r(3,:) = pi/2.d0
	end where

end subroutine cart_sph

subroutine cart_cyl(x,r)
	!----------------------------------------------------------
	! Passage des coord. cartésiennes aux coord. cylindriques
	! CALLING SEQUENCE: call cart_cyl(x,r)
	! INPUT:  x: vecteur position en coordonnées cart.
	! OUTPUT: r: vecteur position en coordonnées cyl. 
	!				(rayon,angle,altitude)			!
	!----------------------------------------------------------

	real(kind=8), dimension(:,:), intent(in)  :: x
	real(kind=8), dimension(:,:), intent(out) :: r
	
	! altitude
	r(3,:)=x(3,:)

	! rayon
	r(1,:) = sqrt(x(1,:)**2+x(2,:)**2)
	
	! angle polaire
	where (x(1,:) .ge. 0.d0 .and. x(2,:) .ne. 0.d0)
		r(2,:) = asin(x(2,:)/r(1,:)) 
	end where 
	where (x(1,:) .lt. 0.d0 .and. x(2,:) .ne. 0.d0)
		r(2,:) = pi - asin(x(2,:)/r(1,:)) 
	end where

	where (x(1,:) .eq. 0.d0 .and. x(2,:).eq. 0.d0)
		r(2,:) = 0.d0
	end where

end subroutine cart_cyl

subroutine sph_cart(r,x)
	!----------------------------------------------------------
	! Passage des coord. sphériques aux coord. cartésiennes 
	! CALLING SEQUENCE: call cart_sph(r,x)
	! INPUT:  r: vecteur position en coordonnées sph.
	!				(rayon,longitude,latitude)	
	! OUTPUT: x: vecteur position en coordonnées cart. 
	!----------------------------------------------------------

	real(kind=8), dimension(:,:), intent(in)  :: r
	real(kind=8), dimension(:,:), intent(out) :: x

	x(1,:) = r(1,:)*cos(r(3,:))*cos(r(2,:)) 
	x(2,:) = r(1,:)*cos(r(3,:))*sin(r(2,:)) 
	x(3,:) = r(1,:)*sin(r(3,:))

end subroutine sph_cart

subroutine cyl_cart(r,x)
	!----------------------------------------------------------
	! Passage des coord. sphériques aux coord. cartésiennes 
	! CALLING SEQUENCE: call cyl_cart(r,x)
	! INPUT:  r: vecteur position en coordonnées sph.
	!				(rayon,angle,altitude)	
	! OUTPUT: x: vecteur position en coordonnées cart. 
	!----------------------------------------------------------

	real(kind=8), dimension(:,:), intent(in)  :: r
	real(kind=8), dimension(:,:), intent(out) :: x

	x(1,:) = r(1,:)*cos(r(2,:))
	x(2,:) = r(1,:)*sin(r(2,:)) 
	x(3,:) = r(3,:)

end subroutine cyl_cart

subroutine sph_cart_vect(Vr,r,Vx)
	!-------------------------------------------------------------------------------------
	! Passage des coord. sphériques aux coord. cartésiennes pour un vecteur en un point r
	! CALLING SEQUENCE: call sph_cart_vect(Vr,r,Vx)
	! INPUT:  r : vecteur position en coord. sph. (rayon,longitude,latitude)	
	!		  Vr: vecteur en coordonnées sph. (Vrayon,Vlongitude,Vlatitude)	
	!				
	! OUTPUT: Vx: vecteur en coordonnées cart. (Vx,Vy,Vz)
	!-------------------------------------------------------------------------------------

	real(kind=8), dimension(:,:), intent(in)  :: r,Vr
	real(kind=8), dimension(:,:), intent(out) :: Vx
	
	Vx(1,:)=Vr(1,:)*cos(r(2,:))*cos(r(3,:))-Vr(2,:)*sin(r(2,:))-Vr(3,:)*sin(r(3,:))*cos(r(2,:))
	Vx(2,:)=Vr(1,:)*sin(r(2,:))*cos(r(3,:))+Vr(2,:)*cos(r(2,:))-Vr(3,:)*sin(r(2,:))*sin(r(3,:))
	Vx(3,:)=Vr(1,:)*sin(r(3,:))+Vr(3,:)*cos(r(3,:))

end subroutine sph_cart_vect
subroutine cyl_cart_vect(Vr,r,Vx)
	!-------------------------------------------------------------------------------------
	! Passage des coord. sphériques aux coord. cartésiennes pour un vecteur en un point r
	! CALLING SEQUENCE: call cyl_cart_vect(Vr,r,Vx)
	! INPUT:  r : vecteur position en coord. cyl. (rayon,angle,altitude)	
	!		  Vr: vecteur en coordonnées sph. (Vrayon,Vangle,Valtitude)	
	!				
	! OUTPUT: Vx: vecteur en coordonnées cart. (Vx,Vy,Vz)
	!-------------------------------------------------------------------------------------

	real(kind=8), dimension(:,:), intent(in)  :: r,Vr
	real(kind=8), dimension(:,:), intent(out) :: Vx
	
	Vx(1,:)=Vr(1,:)*cos(r(2,:))-Vr(2,:)*sin(r(2,:))
	Vx(2,:)=Vr(1,:)*sin(r(2,:))+Vr(2,:)*cos(r(2,:))
	Vx(3,:)=Vr(3,:)

end subroutine cyl_cart_vect
subroutine rot(Vin,theta,Vout)
	!------------------------------------------------------
	! Rotation du vecteur Vin d'un angle theta
	! CALLING SEQUENCE: call rot(Vin,theta,Vout)
	! INPUT:  Vin: vecteur initial en coord. cart.	
	!		  theta: angle en degre
	!				
	! OUTPUT: Vout: vecteur modifié en coordonnées cart.
	!------------------------------------------------------
	real(kind=8), dimension(:,:), intent(in)  :: Vin
	real(kind=8), dimension(:), intent(in)    :: theta
	real(kind=8), dimension(:,:), intent(out) :: Vout
	real(kind=8), dimension(3,3)			  :: matrot
	integer									  :: i
	
	write(*,*) theta
	do i=1,size(theta)
!		matrot(1,1) = cos(theta(i)*pi/180.) ; matrot(1,2) = -sin(theta(i)*pi/180.) ; matrot(1,3)=0.d0
!		matrot(2,1) = sin(theta(i)*pi/180.) ; matrot(2,2) = cos(theta(i)*pi/180.)  ; matrot(2,3)=0.d0
!		matrot(3,1) = 0.d0		 		    ; matrot(3,2) = 0.d0				   ; matrot(3,3)=1.d0

		matrot(1,1) = cos(theta(i)*pi/180.) ; matrot(1,2) =  0.d0		; matrot(1,3)=-sin(theta(i)*pi/180.)
		matrot(2,1) = 0.d0				    ; matrot(2,2) =  1.d0		; matrot(2,3)=0.d0
		matrot(3,1) = sin(theta(i)*pi/180.) ; matrot(3,2) = 0.d0		; matrot(3,3)=cos(theta(i)*pi/180.) 

		Vout(:,i) = matmul(matrot,Vin(:,i)) 

	enddo

end subroutine rot


END MODULE coord