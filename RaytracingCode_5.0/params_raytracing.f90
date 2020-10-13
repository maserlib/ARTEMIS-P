MODULE params_raytracing

	use constantes
	use environ
	use subs_raytracing
	
	implicit none
	
	contains
	
subroutine thetakB(k,V,thkB)
	!-----------------------------------------------------
	! Calcul de l'angle entre k et B
	! CALLING SEQUENCE: call thetakB(k,V,thkB)
	! INPUTS:     k: vecteur d'onde en coord. cart.
	!		      V: vecteur position en coord. cart.
	! OUTPUTS: thkB: angle entre k et B en radian 
	!-----------------------------------------------------
	real(kind=8), dimension(:,:),intent(in) :: k,V
	real(kind=8), dimension(:),intent(out)  :: thkB
	real(kind=8), dimension(3,size(thkb))   :: B
	real(kind=8),dimension(size(thkb))      :: normB,normk
		
	call magneticf(V,B)
	normB=sqrt(B(1,:)**2+B(2,:)**2+B(3,:)**2)
	normk=sqrt(k(1,:)**2+k(2,:)**2+k(3,:)**2)

	where (normB .ne. 0)
		thkB=acos((k(1,:)*B(1,:)+k(2,:)*B(2,:)+k(3,:)*B(3,:))/(normB*normk)) 
	elsewhere
		thkB=0.d0
	end where
	
end subroutine thetakB

subroutine thetarB(V,Vap,thrB)
	!-----------------------------------------------------
	! Calcul de l'angle entre le rayon et B
	! CALLING SEQUENCE: call thetarB(k,V,thrB)
	! INPUTS:     k: vecteur d'onde en coord. cart.
	!		      V: vecteur position en coord. cart.
	! OUTPUTS: thrB: angle entre le rayon et B en radian 
	!-----------------------------------------------------
	real(kind=8), dimension(:,:),intent(in) :: V,Vap
	real(kind=8), dimension(:),intent(out)  :: thrB
	real(kind=8), dimension(3,size(thrb))   :: B
	real(kind=8), dimension(size(thrb))     :: normB,normV
		
	call magneticf(V,B)
	normB=sqrt(B(1,:)**2+B(2,:)**2+B(3,:)**2)
	normV=sqrt((V(1,:)-Vap(1,:))**2+(V(2,:)-Vap(2,:))**2+(V(3,:)-Vap(3,:))**2)
	
	where (normB .ne. 0)
		thrB=acos(((Vap(1,:)-V(1,:))*B(1,:)+(Vap(2,:)-V(2,:))*B(2,:)+(Vap(3,:)-V(3,:))*B(3,:))/(normB*normV)) 
	elsewhere
		thrB=0.d0
	end where
	
end subroutine thetarB

subroutine thetakr(k,V,Y,Vap,thkr)
	!-----------------------------------------------------
	! Calcul de l'angle entre k et le rayon
	! CALLING SEQUENCE: call thetakr(k,V,Vap,thkr)
	! INPUTS:     k: vecteur d'onde en coord. cart.
	!             V: vecteur position en n en coord. cart.
	!             Y: vecteur "magnétique"
	!	    Vap: vecteur position en n+1 en coord. cart.
	! OUTPUTS: thkr: angle entre k et le rayon en radian 
	!-----------------------------------------------------
	real(kind=8), dimension(:,:),intent(in) :: k,V,Vap,Y
	real(kind=8), dimension(:),intent(out)  :: thkr
	real(kind=8), dimension(3,size(thkr))	:: d
	real(kind=8),dimension(size(thkr))      :: normV,normk,normY
	
	d=Vap-V
	normV=sqrt(d(1,:)**2+d(2,:)**2+d(3,:)**2)
	normk=sqrt(k(1,:)**2+k(2,:)**2+k(3,:)**2)
	normY=sqrt(Y(1,:)**2+Y(2,:)**2+Y(3,:)**2)

	where ((normV .ne. 0) .and. (normY /= 0.d0))
		thkr=acos((k(1,:)*d(1,:)+k(2,:)*d(2,:)+k(3,:)*d(3,:))/(normV*normk)) 
	elsewhere
		thkr=0.d0
	end where
	
end subroutine thetakr


subroutine vgroup(V,k,f,thkr,mode,vg)
	!-----------------------------------------------------
	! Calcul de la vitesse de groupe Budden p.140
	!      vg = 1/(n'*cos(alpha))
	! 			n' = d(f*n)df
	!			alpha = angle (k,ray)
	! CALLING SEQUENCE: call vgroup(V,k,f,thkr,mode,vg)
	! INPUTS:     V: vecteur position en coord.cart.
	!			  k: vecteur d'onde en coord. cart.
	!			  f: frŽquence de l'onde
	!			  mode: mode de porpagation
	!			  thkr: angle (k,rayon)
	! OUTPUTS: 	  vg: vitesse de groupe
	!-----------------------------------------------------
	real(kind=8), dimension(:,:), intent(in)	:: V,k
	real(kind=8), dimension(:), intent(in)		:: thkr
	real(kind=8), intent(in) 					:: f
	character(len=1), intent(in), dimension(:)  :: mode	
	real(kind=8), dimension(:), intent(out)		:: vg
	
	real(kind=8), dimension(3,size(V,2)) 		:: Y, Yav, Yap
	real(kind=8), dimension(size(V,2))			:: X,Xav,Xap,mu, muav, muap,muprime
	real(kind=8)								:: df
	
	df=1.d0 !Hz

	call varplasma(V,f-df,Xav,Yav)
	call varplasma(V,f+df,Xap,Yap)
	call varplasma(V,f,X,Y)
	
	call index(Xav,Yav,k,mode,muav)
	call index(Xap,Yap,k,mode,muap)
	call index(X,Y,k,mode,mu)
	
	muprime = mu +f*(muap-muav)/(2.d0*df)
	
	vg(:) = 1.d0/(muprime(:)*cos(thkr(:)))	
	
end subroutine vgroup


subroutine vphase(mu,vp)
	!-----------------------------------------------------
	! Calcul de la vitesse de phase
	!      vp = c/mu 
	! CALLING SEQUENCE: call vgroup(mu,vp)
	! INPUTS:     mu: indice
	! OUTPUTS: 	  vp: vitesse de phase
	!-----------------------------------------------------
	real(kind=8), dimension(:), intent(in) 	:: mu
	real(kind=8), dimension(:), intent(out) :: vp
	
	vp = 1.d0/mu

end subroutine vphase

subroutine distance(V,Vav,dist,numray)
	!-----------------------------------------------------
	! Calcul de la distance parcourue par le rayon depuis son pt de dŽpart
	!
	! CALLING SEQUENCE: call distance(d,V,dist)
	! INPUTS:     d: distance au pas d'avant
	!             V: vecteur position en coord. cart en position n
	!			  Vav: vecteur position en coord. cart en position n-1
	! OUTPUTS: 	  dist: distance
	!-----------------------------------------------------
	real(kind=8), dimension(:,:), intent(in) :: V, Vav
	integer(kind=8), dimension(:), intent(in)   :: numray
	real(kind=8), dimension(:), intent(inout):: dist
	
	dist(numray) = dist(numray) + sqrt((V(1,:)-Vav(1,:))**2+(V(2,:)-Vav(2,:))**2+(V(3,:)-Vav(3,:))**2)
	
end subroutine distance


subroutine fplasma(f,X,fp)
	!-----------------------------------------------------
	! Calcul de la frŽquence plasma fp = f*sqrt(X)
	!
	! CALLING SEQUENCE: call fplasma(f,X,fp)
	! INPUTS:     f: frŽquence de l'onde
	!             X: (fp/f)^2
	! OUTPUTS: 	  fp: frequence plasma 
	!-----------------------------------------------------
	real(kind=8), intent(in)  				:: f
	real(kind=8), dimension(:), intent(in)  :: X
	real(kind=8), dimension(:), intent(out) :: fp
	
	fp(:) = f*sqrt(X(:))

end subroutine fplasma

subroutine fcyclo(f,Y,fc)
	!-----------------------------------------------------
	! Calcul de la frŽquence cyclotron fc = f*norm(Y)
	!
	! CALLING SEQUENCE: call fcyclo(f,X,fp)
	! INPUTS:     f: frŽquence de l'onde
	!             Y: fc/f
	! OUTPUTS: 	  fc: frequence plasma 
	!-----------------------------------------------------
	real(kind=8), intent(in)  				:: f
	real(kind=8), dimension(:,:), intent(in):: Y
	real(kind=8), dimension(:), intent(out) :: fc
	
	fc(:) = f*sqrt(Y(1,:)**2+Y(2,:)**2+Y(3,:)**2)
	
end subroutine fcyclo

subroutine fuph(f,X,Y,fuh)
	!-----------------------------------------------------
	! Calcul de la frŽquence uper hybride fuh = sqrt(fp^2+fc^2)
	!
	! CALLING SEQUENCE: call fuh(f,X,Y,fuh)
	! INPUTS:     f: frŽquence de l'onde
	!             X: (fp/f)^2
	!             Y: fc/f
	! OUTPUTS: 	  fc: frequence plasma 
	!-----------------------------------------------------
	real(kind=8), intent(in)  				:: f
	real(kind=8), dimension(:), intent(in)  :: X
	real(kind=8), dimension(:,:), intent(in):: Y
	real(kind=8), dimension(:), intent(out) :: fuh
	real(kind=8), dimension(size(fuh)) 		:: fp,fc
	
	call fplasma(f,X,fp)
	call fcyclo(f,Y,fc)
	
	fuh = sqrt(fp**2+fc**2)
	
end subroutine fuph


subroutine polar_ratio(mode,f,k,V,LP, pratio,pratio_1)
	!-----------------------------------------------------
	! Calcul du ratio axial de polarisation: (Quemada)
	!	repre: (x',y',z') avec z' // k et B ds (x',z')
	!
	!	S = 1-X(1-Y^2) ; D = XY/(1-Y^2) ; P = 1-X
	!
	!	p = Ey'/Ex' = iPcos(th)(n^2-S)/D(n^2sin^2(th)-P)
	!	q = Ez'/Ex' = -(n^2-P)sin(th)/Pcos(th)(n^2-S)
	!	r = Ez'/Ey' = q/p
	!
	!	ATTENTION: BUDDEN Y=eB/mf avec e=-1.6e-19C donc \Theta_Budden = angle(k,Y) = pi - angle(k,B)
	!
	! CALLING SEQUENCE: call polar_ratio(mode,f,k,V,pratio)
	! INPUTS:     mode: mode de propagation ('X', 'O' ou 'V')
	!			  f: frŽquence de l'onde
	!             k: vecteur d'onde en coord. cart.
	!		      V: vecteur position en coord. cart.
	!		 pratio: rapport de polarisation ˆ l'etape n
	!			 LP: vecteur limitong polarization ('VP' -> variable polar, 'LP' -> limiting polar
	! OUTPUTS: 	  pratio_1: p,q,r (dim 3*Nray)
	!-----------------------------------------------------
	character(len=1), intent(in),dimension(:)		:: mode
	real(kind=8), intent(in)  				 		:: f
	real(kind=8), dimension(:,:), intent(in)  		:: V,k
	character(len=2),dimension(:), intent(in)		:: LP
	complex(kind=8), dimension(:,:), intent(in) 	:: pratio
	complex(kind=8), dimension(:,:), intent(out) 	:: pratio_1
	
	real(kind=8), dimension(size(V,2))   :: X, S,D,P, mu, th, normY, sinth,costh,racine,mu2, tanth
	real(kind=8), dimension(size(V,2))   :: sinth4, sinth2,costh2, nmode
	real(kind=8), dimension(3,size(V,2)) :: Y
	integer								 :: i
	
		call varplasma(V,f,X,Y)
		normY = sqrt(Y(1,:)**2+Y(2,:)**2+Y(3,:)**2)
		call index(X,Y,k,mode,mu)
		call thetakb(k,V,th)
	
		nmode(:)=1.d0
		where (mode .eq. 'X') 
			nmode=-1.d0
		end where
	
		S = 1-X/(1-normY**2)
		D = X*normY/(1-normY**2)
		P = 1-X
		
		where (LP .eq. 'VP')	
			where (normY /= 0.d0) ! clcul de polarisation
        			!p = Ey'/Ex'
        			pratio_1(1,:) = -ic*(D*(mu**2*sin(th)**2-P))/(P*cos(th)*(mu**2-S))
		
	        		!q = Ez'/Ex'
        			pratio_1(2,:) = -(mu**2-P)*sin(th)/(P*cos(th))
		
	        		!r = Ez'/Ey'
        			pratio_1(3,:) = -ic*(mu**2-S)*(mu**2-P)*sin(th)/(D*(mu**2*sin(th)**2-P))
			
			elsewhere !polarisation circulaire si B = 0
				pratio_1(1,:) = nmode*ic
				pratio_1(2,:) = 0.d0
				pratio_1(3,:) = 0.d0
			end where
			!DL autour de th=90¡ en mode O
			!-------------------------------

			!	where (mode .ne. 'X' .and. ((th .gt. pi/2.-0.0174533) .and. (th .lt. pi/2.+0.0174533)))
			!	where (((th .gt. pi/2.-0.0174533) .and. (th .lt. pi/2.+0.0174533)))
			!		sinth = 1.d0-(pi/2.d0-th)**2/2.d0
			!		costh = (pi/2.d0 -th)-(pi/2.d0-th)**3/6.d0
			!		costh2 = (pi/2.d0 -th)**2
			!		mu2 = 1-X/(1-X*costh2)
			!		pratio(1,:) = ic*(D/(P*(mu2-S)))*costh2
			!		pratio(2,:) = -((mu2-P)/P)*tan(th)
			!		pratio(3,:) = -ic*(mu2-S)*(mu2-P)*tan(th)**2/(D*sinth)
			!	end where
		else where
			pratio_1(1,:) = pratio(1,:)
			pratio_1(2,:) = pratio(2,:)
			pratio_1(3,:) = pratio(3,:)		
		end where
	
end subroutine polar_ratio

subroutine limiting_polar(f,V,k,V_1,k_1,dr,unit,LP)
	!-----------------------------------
	! Drapeau de limiting polarization
	! si abs(mux-muo) < c/2pif*dmux/dr ou c/2pif*dmux/dr alors la polarisation est figŽe
	!
	! ATTENTION AUX CHOIX DE DR ET AUX UNITES DE V (m,km,Rs,...)
	!-----------------------------------
	real(kind=8), intent(in)  				 		:: f !en Hz
	real(kind=8), dimension(:,:), intent(in)  		:: V,k,V_1,k_1
	real(kind=8), dimension(:),intent(in)			:: dr
	character(len=2), intent(in)                    :: unit
	character(len=2),dimension(:), intent(out)		:: LP
	
	character(len=1), dimension(size(V,2))			:: modeX,modeO
	real(kind=8), dimension(size(V,2))				:: X, mux,muo, mux_1,muo_1, dmux,dmuo,dmuxo
	real(kind=8),dimension(3,size(V,2))				:: Y, Vdr
	real(kind=8)									:: dr2, convunit
	integer 										:: i
	
	modeX(:)='X'
	modeO(:)='O'
	
	if ((unit .eq. 'KM') .or. (unit .eq. 'km')) then
		convunit = 1.d3
	else if ((unit .eq. 'RS') .or. (unit .eq. 'rs') .or. (unit .eq. 'Rs')) then
		convunit = RS*1.d3
	else if ((unit .eq. 'RSol') .or. (unit .eq. 'rsol') .or. (unit .eq. 'Rsol') .or. (unit .eq. 'RSOL')) then
		convunit=RSol*1.d3
	else if ((unit .eq. 'RT') .or. (unit .eq. 'rt') .or. (unit .eq. 'Rt')) then
		convunit=RT*1.d3
	else if ((unit .eq. 'RJ') .or. (unit .eq. 'rj') .or. (unit .eq. 'RJ')) then
		convunit=RJ*1.d3
	else
		write(*,*) 'PB UNITE'
		convunit = 1.d0
	endif
	
	!calcul de mux et muo en (r_n, k_n)
	call varplasma(V,f,X,Y)
	call index(X,Y,k,modeX,muX)
	call index(X,Y,k,modeO,muO)

	!calcul de mux et muo en (r_n+1,k_n+1)
	call varplasma(V_1,f,X,Y)
	call index(X,Y,k_1,modeX,muX_1)
	call index(X,Y,k_1,modeO,muO_1)

	!dmu = c/(2pif)*(mu1-mu)/dr !ATTENTION AUX UNITES
	dmux = c/(2.d0*pi*f)*abs(muX_1-muX)/(dr*convunit)
	dmuo = c/(2.d0*pi*f)*abs(muO_1-muO)/(dr*convunit)

	!calcul de muo_1-mux_1
	dmuxo = abs(muo_1-mux_1)	
		
	!Limiting polarization qd (muo-mux) < dmuo OU dmux < 1%
	where((abs(dmuxo) .le. dmux) .or. (abs(dmuxo) .le. dmuo))
		LP(:) ='LP' !limiting polarization
	else where
		LP(:)='VP'  !variable polarization
	end where

	do i=1,size(X)
		write(19,*) 'mux',mux_1(i),'muo',muo_1(i),'dmuxo',dmuxo(i), 'dmux', dmux(i), 'dmuo',dmuo(i), 'Limiting Polar ', LP(i)
		!write(19,*) '----------------------------------------------------------'
	enddo
		!write(19,*) '----------------------------------------------------------'


end subroutine limiting_polar

subroutine errordisp(k,V,X,Y,mu,errdisp)
	!---------------------------------
	! calcul de l'erreur k^2-1-X*(pY-1)
	! 	d'aprs Haselgrove (1960)
	! CALLING SEQUENCE: call error(k,X,Y, err)
	! INPUTS : V: vect position en coord. cart.
	!          k: vect d'onde en coord. cart.
	!		   f: frequence en Hz
	! OUTPUTS: err: scalaire = 0
	!------------------------------------
	real(kind=8), dimension(:,:), intent(in)  :: Y,k,V
	real(kind=8), dimension(:), intent(in)    :: X,mu
	real(kind=8), dimension(:), intent(out)   :: errdisp
	real(kind=8), dimension(size(X))          :: A,B,C, normY,p,s,d, th
	
	normY = sqrt(Y(1,:)**2+Y(2,:)**2+Y(3,:)**2)
	!normk2 = k(1,:)**2+k(2,:)**2+k(3,:)**2
	!call index(X,Y,k,mode,mu)
			!si B=0 (Y=0) alors direction de Y = verticale (vdotk ne pas tre nul)
	!where (normY .eq. 0.d0)  
!		vdotk(:) = maxval(k(:,:))
!	elsewhere
!		vdotk(:) = (Y(1,:)*k(1,:)+Y(2,:)*k(2,:)+Y(3,:)*k(3,:))/normY(:)
!	end where
!	
!	AH = 1-X-normY**2+X*normY**2*(vdotk)**2/normk2
!	!BH = 2B_60
!	BH = -2*(1-X)*(1-x-normY**2)+X*normY**2-X*normY**2*(vdotk)**2/normk2
!	CH = (1-X)*((1-X)**2-normY**2)
!	
!	errdisp = abs(AH*mu**4+BH*mu**2+CH)
!
p=1.d0-x
s=1.d0-x/(1.d0-normY**2)
d=x*normY/(1.d0-normY**2)

call thetakb(k,V,th)
A = p*cos(th)**2+s*sin(th)**2
B = (d**2-s**2)*sin(th)**2-p*s*(1+cos(th)**2)
C = p*(s**2-d**2)

errdisp = A*mu**4+B*mu**2+C

!write(*,*) errdisp

end subroutine errordisp

END MODULE params_raytracing
