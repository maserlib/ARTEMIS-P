PROGRAM raytracing
	
	use constantes
	use coord
	use subs_raytracing
	use environ
	use RungeKuttaFehlberg
	use RungeKutta
	use params_raytracing
	
	implicit none
	
	!DECLARATIONS
	integer(kind=8)                           	 :: Niter
	real(kind=8), dimension(:,:), allocatable 	 :: r,k,r_1,k_1,B,Y,r0
	complex(kind=8), dimension(:,:), allocatable :: pratio,pratio_1
	logical, dimension(:), allocatable	  	  	 :: mask
	real(kind=8), dimension(:), allocatable	  	 :: kiter,dt,dt_1,mu,mu_1,X,normk,err,errdisp,vp,thkb, thkr, thrb, dr
	real(kind=8), dimension(:), allocatable	  	 :: vg, dist, fp, fc, fuh, Nel, normB, mux, muo, dmux, dmuo, mux_1, muo_1
	integer(kind=8),dimension(:),allocatable  	 :: numray, numraytemp, flagmu, flagnan,flagerr, flagerrdisp
	real(kind=8)                              	 :: f,t1,t2, mumax, dtinit,dtmax,dtmin
	integer(kind=8)					          	 :: Nray,Nray2,Nray3 ,i ,j
	character(len=1), dimension(:), allocatable	 :: mode, mode_1, modeX, modeO
	character(len=2), dimension(:), allocatable	 :: LP
	character(len=4)				          	 :: integ, path, sysco
	character(len=2)                             :: unit
	character(len=1)							 :: erase
	character(len=24)				          	 :: filen
	logical									  	 :: is_existing

	call cpu_time(t1)



! INITIALISATION (fichier init_raytracing.txt')
	open(21, file='init_raytracing.txt')
	
	!Nom du dossier
	read(21,*) path
	write(*,*) 'MODELE:', path
	write(*,*) '----------------------------------'
	!FICHIER ECRITURE TEST
	!ecriture fichier
	inquire(directory=trim('../'//path), exist=is_existing)
	if (.not.is_existing) then
		write(*,*) is_existing, 'fichier existe pas'
		call system("mkdir "//'../'//path)
		open(19,file='../'//path//'/test_limitingpolar.txt')
	else
		write(*,'(a,$)') 'folder already exists. Erase?  (Y/N)'
		read(*,*) erase
		if (erase == 'Y') then 
			open(19,file='../'//path//'/test_limitingpolar.txt')
		else if (erase == 'N') then
			stop 'STOP folder already exists'
		else
			stop 'STOP incorrect answer'
		end if
	end if
		
	!Choix de l'integrateur
	read(21,*) integ
	write(*,*) 'INTGEGRATEUR ', integ	
	
	!initialisation frequence (kHz-->Hz)
	read(21,*) f
	f = f*1.d3
	write(*,*) 'frequence', f, ' Hz'
	
	!Nombre de rayons
	read(21,*) Nray
	write(*,*) 'Nombre de rayons:', Nray
	Nray2=Nray
	
	!Nombre d'iteration
	read(21,*) Niter
	write(*,*) 'Niter ',Niter
	write(*,*) '----------------------------------'
	
	!---------------------------------------------------------------------------
	!ALLOCATION MEMOIRE
	
	allocate(r(3,Nray),r0(3,Nray),r_1(3,Nray),k(3,Nray),k_1(3,Nray),B(3,Nray),Y(3,Nray),vg(Nray),pratio(3,Nray),pratio_1(3,Nray))
	allocate(X(Nray),thkb(Nray),mu(Nray),mu_1(Nray),normk(Nray),dt(Nray),dt_1(Nray), mode(Nray), mode_1(Nray))
	allocate(kiter(Nray),err(Nray), errdisp(Nray),vp(Nray), thkr(Nray), thrb(Nray),dist(Nray),mask(Nray),numray(Nray),numraytemp(Nray))
	allocate(flagmu(Nray),flagnan(Nray),flagerr(Nray),flagerrdisp(Nray))
	allocate(fp(Nray), fc(Nray), fuh(Nray), Nel(Nray), normB(Nray))
	allocate(mux(Nray),muo(Nray), dmux(Nray),dmuo(Nray),mux_1(Nray),muo_1(Nray),modeX(Nray),modeO(Nray), dr(Nray), LP(Nray))
	!---------------------------------------------------------------------------
	!vecteur indice de rayon
	forall(i=1:Nray) numray(i) = i
	
	!Système de coordonnées initial
	read(21,*) sysco
	write(*,*) 'SYSTEME DE COORD.:  ', sysco
	
	!Unite de longueur
	read(21,*) unit
	write(*,*) 'UNITE DE LONGUEUR: ', unit
	
	!initialisation position
	do i=1,Nray 
		read(21,*) r0(:,i)
	enddo
	
	if (sysco .eq. 'sphr') then
		r0(2,:) = r0(2,:)*pi/180.d0 ; r0(3,:) = r0(3,:)*pi/180.d0
		write(*,*) 'r ', r0(1,1), 'long. ', r0(2,1)*180.d0/pi, 'lat. ', r0(3,1)*180.d0/pi 
		call sph_cart(r0,r)
		write(*,*) r(:,1)
	else
		r = r0
		write(*,*) r(:,1)		
	endif
	
	!initialisation distance
	dist(:) = 0.d0
	
	!initialistaion vecteur d'onde
	do i=1,Nray 
		read(21,*) k(:,i)
	enddo
	write(*,*) 'k', k

	!Choix du mode
	do i=1,Nray
		read(21,*) mode(i)
	enddo
	write(*,*) 'MODE ',mode(1)


	!initialisation indice
	call varplasma(r,f,X,Y)
	write(*,*) 'X',X(1)
	write(*,*) 'Y',Y(:,1)

	call thetakB(k,r,thkb)
	write(*,*) 'thkB', thkB*180./pi
!	call index(X,Y,thkb,mode,mu)
	call index(X,Y,k,mode,mu)
	write(*,*) 'mu', mu(1)
	
	!lecture de l'indice max
	read(21,*) mumax
	write(*,*) 'mu max', mumax
	
	!normalisation de k par normk = mu
	normk = sqrt(k(1,:)**2+k(2,:)**2+k(3,:)**2)
	forall (i=1:3) k(i,:) = (k(i,:)/normk(:))*mu(:)
	write(*,*) 'k normalise', k

	
	!initialisation du pas de temps
	read(21,*) dtinit
	read(21,*) dtmax
	read(21,*) dtmin	
	write(*,*) 'dt init', dtinit, 'dt max', dtmax, 'dt min', dtmin
	dt(:)=dtinit
	
	!initialisation arbitraire de pratio en polar circulaire
	where (mode .eq. 'X')
		pratio(1,:) = -ic
	elsewhere
		pratio(1,:) = ic
	end where
	pratio(2,:) = 0.d0
	pratio(3,:) = 0.d0
	!initialisation arbitraire de LP='VP'
	LP(:)='VP'	
	
	!calcul de certains paramètres initaux
	call density(r, Nel)
	call vphase(mu,vp)
	call thetakb(k,r,thkb)
	call fplasma(f,X,fp)
	call fcyclo(f,Y,fc)
	call fuph(f,X,Y,fuh)
	call polar_ratio(mode, f,k,r,Lp,pratio,pratio_1)
	
	!ECRITURE FICHIER VISU
	open(18, file='../'//path//'/visu_param.txt')
		write(18,"(A4)") path
		write(18,"(I3.1)") Nray
		write(18,"(A4)") integ
		write(18,"(A2)") unit
		write(18,"(1f10.5)") f*1.d-3
		write(18,"(I7.5)") Niter
		do i=1,Nray
			write(18,"(A1)") mode(i)
		enddo
	close(18)	

!SAUVEGARDE DES PARAMETRES D'INITIALISATION
	
	if (integ .eq. 'RK4T') then
		filen='../'//path//'/dataRK4T_' 
	else 
		filen='../'//path//'/dataRK4F_'
	end if
		
	do i=1,nray
		write(filen(18:20),'(i2.2)') i
		filen(20:24) = '.dat'
		open(21+i,file=filen, form='formatted')
		write(21+i,"(A300)") 'x (km)    y(km)    z (km)    kx    ky    kz    indice    pas    vg (c)    vp (c)  &
			&  (k.B) (deg)    (r.B) (deg)    distance (km)    Bx (T)    By (T)    Bz (T)  &
			&  fp (Hz)    fc (Hz)    fUH (Hz)    Ne (cm-3)   &
			&  Re(Ey/Ex)    Im(Ey/Ex)    Re(Ez/Ex)    Im(Ez/Ex)    Re(Ez/Ey)    Im(Ez/Ey)     mode      LP      err      errdisp  '

	write(21+i,9991), r(:,i),k(:,i),mu(i),dt(i), -1.d0,vp(i), thkb(i)*180.d0/pi, -1.d0, &
			& 0.d0, B(:,i), fp(i), fc(i), fuh(i),Nel(i),&
			& real(pratio_1(1,i)), aimag(pratio_1(1,i)),&
			& real(pratio_1(2,i)), aimag(pratio_1(2,i)),&
			& real(pratio_1(3,i)), aimag(pratio_1(3,i)),&
			& mode(i), LP(i)
	enddo

	!initialisation kiter
	kiter(:) = 1.d0

	!BOUCLE

	do while (sum(kiter) .ne. Nray*Niter*1.d0)
		if (integ .eq. 'RK4T') then
			call RK4(r,k,dt,dxdt,dkdt,f,mode,r_1,k_1)
		else
			call RK4F(r,k,dt,dxdt,dkdt,f,mode,dtmin,dtmax,r_1,k_1)
		end if

		call varplasma(r_1,f,X,Y)
		call index(X,Y,k_1,mode,mu_1)
		!renormalisation de k a mu
		normk = sqrt(k_1(1,:)**2+k_1(2,:)**2+k_1(3,:)**2)
		forall (i=1:3) k_1(i,:) = (k_1(i,:)/normk(:))*mu_1(:)
		
		!calcul de dr=sqrt(r_1^2-r^2)
		dr(:)=sqrt(r_1(1,:)**2+r_1(2,:)**2+r_1(3,:)**2)-sqrt(r(1,:)**2+r(2,:)**2+r(3,:)**2)
		
		!Calcul des parametres
		call density(r_1,Nel)
		call vphase(mu_1,vp)
		call thetakb(k_1,r_1,thkb)
		call thetakr(k,r,Y,r_1,thkr)
		call vgroup(r_1,k_1,f, thkr,mode,vg)
		call thetarB(r,r_1,thrb)
		call distance(r,r_1,dist,numray)
		call magneticf(r_1,B)		
		call fplasma(f,X,fp)
		call fcyclo(f,Y,fc)
		call fuph(f,X,Y,fuh)
		!call limiting_polar(f,r,k,r_1,k_1,dr,unit,LP)
		!call polar_ratio(mode, f,k_1,r_1,LP,pratio,pratio_1)

		
		!adaptation du pas de temps		
		!---------------------------------
		where (mu_1 .ne. mu) 
			dt_1 = dtinit*abs(dr(:)/(mu_1(:)-mu(:)))
		elsewhere 
			dt_1 = dtmax
		end where

		where (dt_1 > dtmax) 
			dt_1 = dtmax
		end where
		
		where (dt_1 < dtmin) 
			dt_1 = dtmin
		end where
		!Conditions de sortie		
		!----------------------------------
		kiter(numray) = kiter(numray) + 1.d0
		
			! mu > mumax
			flagmu(:)=0
			where (mu .ge. mumax)  
				kiter(numray) = Niter*1.d0 
				flagmu = numray
			end where
			if(sum(flagmu) .ne. 0) then
				write(*,*) flagmu, 'STOP MU'
			endif
	
			!une coord. de r = NaN
			flagnan(:)=0
			where (((r(1,:)+1)/(r(1,:)+1) .ne. 1) .or. ((r(2,:)+1)/(r(2,:)+1) .ne. 1) .or. ((r(3,:)+1.d0)/(r(3,:)+1.d0) .ne. 1))
				kiter(numray) = Niter*1.d0
				flagnan = numray
			end where
			if(sum(flagnan) .ne. 0) then
				write(*,*) flagnan, 'STOP NaN'
			endif
		
			!erreur relative entre (n^2-1) et X(pY-1) > 1% => dt = dtmin
			flagerr(:)=0
			call error(k_1,X,Y,mode,err)
			where (err .gt. 1.d-2)  
				dt_1(numray) = dtmin
			end where
			!erreur dispersion > 1% => dt = dtmin
			flagerrdisp(:)=0
			call errordisp(k_1,r_1,X,Y,mu_1,errdisp)
			where (errdisp .gt. 1.d-2)  
				dt_1(numray) = dtmin
			end where
	
		mode_1=mode
		
			
		
		!Selection des paramètres des rayons qui continuent	
		mask = kiter(numray) .lt. Niter
		Nray3 = count(mask)

		deallocate (r,k,dt,mu,numraytemp, mode, mux, muo,pratio)
		allocate (r(3,Nray3), k(3,Nray3),dt(Nray3),mu(Nray3),numraytemp(Nray3), mode(Nray3), mux(Nray3), muo(Nray3),pratio(3,Nray3))

		j = 1
		do i=1,Nray2
			if (mask(i)) then
				!écriture dans le fichier
				write(21+numray(i),9991), r_1(:,i),k_1(:,i),mu_1(i),dt_1(i),&
					& vg(i), vp(i), thkb(i)*180.d0/pi, thrb(i)*180.d0/pi,&
					& dist(numray(i)), B(:,i), &
					& fp(i), fc(i), fuh(i),Nel(i),&
					& real(pratio_1(1,i)), aimag(pratio_1(1,i)),&
					& real(pratio_1(2,i)), aimag(pratio_1(2,i)),&
					& real(pratio_1(3,i)), aimag(pratio_1(3,i)),&
					& mode_1(i),&! abs(mux_1(i)-muo_1(i)), dmux(i), dmuo(i),&
					& LP(i), err(i), errdisp(i)								
								
				r(:,j) = r_1(:,i)
				k(:,j) = k_1(:,i)
				dt(j) = dt_1(i)
				mu(j) = mu_1(i)
				numraytemp(j) = numray(i)
				mode(j) = mode_1(i)
				mux(j) = mux_1(i)
				muo(j) = muo_1(i)
				pratio(:,j) = pratio_1(:,i)
				j = j+1
			end if
		end do

		9991 format(26(1e17.10,1x),A1,1x,A2,1x,1e17.10,1x,1e17.10)

		deallocate(numray)
		allocate(numray(Nray3))
		numray=numraytemp
		Nray2=Nray3
		deallocate(r_1,k_1,dt_1,mu_1,X,Y,normk,err, errdisp,vp,thkb,thkr,thrb,vg,B,mask, flagmu, flagnan,flagerr,flagerrdisp)
		deallocate(fp,fc,fuh,Nel,pratio_1)
		deallocate(mode_1, mux_1, muo_1, dmux, dmuo, modeX,modeO,dr,LP)
		allocate(r_1(3,Nray3),k_1(3,Nray3),Y(3,Nray3), B(3,Nray3),pratio_1(3,Nray3))
		allocate(dt_1(Nray3),mu_1(Nray3),X(Nray3),normk(Nray3),err(Nray3),errdisp(Nray3),vp(Nray3),vg(Nray3),LP(Nray3))
		allocate(thkb(Nray3),thkr(Nray3),thrb(Nray3),mask(Nray3), flagmu(Nray3), flagnan(Nray3),flagerr(Nray3),flagerrdisp(Nray3))
		allocate(fp(Nray3), fc(Nray3),fuh(Nray3), Nel(Nray3), mode_1(Nray3), mux_1(Nray3), muo_1(Nray3), dmux(Nray3),dmuo(Nray3))
		allocate(modeX(Nray3),modeO(Nray3), dr(Nray3))
	end do
	
	call cpu_time(t2)
	write(*,*) '----------------------------------'
	write(*,*) 'temps de calcul (s)', t2-t1

	close(21) ; close(20) ; close(19)
	deallocate(r_1,k_1,dt_1,mu_1,X,Y,normk,err,vp,thkb,thkr,thrb,vg,B,mask, flagmu, flagnan,flagerr,dr)
	deallocate(fp,fc,fuh, kiter,pratio,pratio_1,Nel, mode, mode_1, mux, muo, mux_1, muo_1,dmux,dmuo,modeX,modeO,LP)
END PROGRAM raytracing
