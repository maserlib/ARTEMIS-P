PROGRAM raytracing

  use constantes
  use coord
  use subs_raytracing
  use environ
  use RungeKuttaFehlberg
  use RungeKutta
  use params_raytracing
  use variables
  
  implicit none
  
  !DECLARATIONS
  integer(kind=8)                              :: Niter
  real(kind=8), dimension(:,:), allocatable    :: r,k,r_1,k_1,B,Y,r0
  complex(kind=8), dimension(:,:), allocatable :: pratio,pratio_1
  logical, dimension(:), allocatable           :: mask
  real(kind=8), dimension(:), allocatable       :: kiter,dt,dt_1,mu,mu_1,X,normk,err,errdisp,vp,thkb, thkr, thrb, dr
  real(kind=8), dimension(:), allocatable       :: vg, dist, fp, fc, fuh, Nel, normB, mux, muo, dmux, dmuo, mux_1, muo_1
  integer(kind=8),dimension(:),allocatable     :: numray, numraytemp, flagmu, flagnan,flagerr, flagerrdisp
  real(kind=8)                                 :: f,t1,t2, mumax, dtinit,dtmax,dtmin
  real                                         :: r_planet2
  integer                                      :: Nray2,Nray3 ,i ,j, icmd,ncmd
  character(len=1), dimension(:), allocatable  :: mode, mode_1, modeX, modeO
  character(len=2), dimension(:), allocatable  :: LP
  character(len=4)                             :: integ, path, sysco
  character(len=2)                             :: unit
  character(len=1)                             :: erase
  character(len=19)                            :: filen
  character(len=53):: init_file
  character(len=50) :: arg
  character(len=500) :: cmd
  logical                                       :: is_existing

  call cpu_time(t1)


 ! call read_environ()
  call read_environ_netcdf()
  call read_magnetic_field_netcdf()
  r_planet2 = r_planet*r_planet
! INITIALISATION (file init_raytracing.txt')
 
  call get_command(cmd)
  ncmd = command_argument_count()

  icmd=1
  do while (icmd<=ncmd)
    call get_command_argument(icmd,arg)
    write(*,*)arg
    select case(arg)
    case("-f","--file")
      icmd=icmd+1
      call get_command_argument(icmd,arg)
      init_file=arg
    case default
      write(*,*)"init_raytracing file not specified"
      stop
    end select
    icmd=icmd+1
  end do

  write(*,*)'init_ray file :',init_file

  open(21, file=trim(init_file)) !'init_raytracing.txt')
  
  ! results directory name
  read(21,*) path
!  write(*,*) 'MODELE:', path
!  write(*,*) '----------------------------------'
  !FICHIER ECRITURE TEST
  !ecriture fichier
  path=""
!    inquire(file=trim('./'//path), exist=is_existing) !SYNTAXE POUR GFORTRAN
!    !inquire(directory=trim('../'//path), exist=is_existing) !SYNTAXE POUR IFORT
!  if (.not.is_existing) then
!    write(*,*) is_existing, 'file don t exist'
!    call system("mkdir "//'./'//path)
!    open(19,file='./'//path//'/test_limitingpolar.txt')
!  else
!    write(*,'(a,$)') 'folder already exists. Erase?  (Y/N)'
!    read(*,*) erase
!    if (erase == 'Y') then 
!      open(19,file='./'//path//'/test_limitingpolar.txt')
!    else if (erase == 'N') then
!      stop 'STOP folder already exists'
!    else
!      stop 'STOP incorrect answer'
!    end if
!  end if
  ! integrator choice
  read(21,*) integ
  write(*,*) 'INTGEGRATOR ', integ  
  
  ! frequency initialisation (kHz-->Hz)
  read(21,*) f
  f = f*1.d3
  write(*,*) 'frequency', f, ' Hz'
  
  !Number of rays
  read(21,*) Nray
  write(*,*) 'Number of rays:', Nray
  Nray2=Nray
  
  !Numbre of iterations
  read(21,*) Niter
  write(*,*) 'Niter ',Niter
  write(*,*) '----------------------------------'
  
  !---------------------------------------------------------------------------
  ! MEMORY ALLOCATION
  
  allocate(r(3,Nray),r0(3,Nray),r_1(3,Nray),k(3,Nray),k_1(3,Nray),B(3,Nray),Y(3,Nray),vg(Nray),pratio(3,Nray),pratio_1(3,Nray))
  r(:,:)=0.; r0(:,:)=0.; r_1(:,:)=0.; k(:,:)=0.; k_1(:,:)=0.; B(:,:)=0; Y(:,:)=0.; vg(:)=0.; pratio(:,:)=0.; pratio_1(:,:)=0. 
  allocate(X(Nray),thkb(Nray),mu(Nray),mu_1(Nray),normk(Nray),dt(Nray),dt_1(Nray), mode(Nray), mode_1(Nray))
  X(:)=0.; thkb(:)=0.; mu(:)=0.; mu_1(:)=0.; normk(:)=0.; dt(:)=0.; dt_1(:)=0.; mode(:)=''; mode_1(:)=''
  allocate(kiter(Nray),err(Nray), errdisp(Nray),vp(Nray), thkr(Nray), thrb(Nray),dist(Nray),mask(Nray),numray(Nray),&
          numraytemp(Nray))
  kiter(:)=0.; err(:)=0.; errdisp(:)=0.; vp(:)=0.; thkr(:)=0.; thrb(:)=0.; dist(:)=0.; mask(:)=0; numray(:)=0; numraytemp(:)=0
  allocate(flagmu(Nray),flagnan(Nray),flagerr(Nray),flagerrdisp(Nray))
  flagmu(:)=0; flagnan(:)=0; flagerr(:)=0; flagerrdisp(:)=0
  allocate(fp(Nray), fc(Nray), fuh(Nray), Nel(Nray), normB(Nray))
  fp(:)=0.; fc(:)=0.; fuh(:)=0.; Nel(:)=0.; normB(:)=0.
  allocate(mux(Nray),muo(Nray), dmux(Nray),dmuo(Nray),mux_1(Nray),muo_1(Nray),modeX(Nray),modeO(Nray), dr(Nray), LP(Nray))
  mux(:)=0.;muo(:)=0.; dmux(:)=0.; dmuo(:)=0.; mux_1(:)=0.; muo_1(:)=0.; modeX(:)=''; modeO(:)=''; dr(:)=0.; LP(:)=''
  !---------------------------------------------------------------------------
  ! index vector for rays
  forall(i=1:Nray) numray(i) = i
  
  ! initial coordinate system
  read(21,*) sysco
  write(*,*) 'COORDINATE SYSTEM : ', sysco
  
  ! length unit
  read(21,*) unit
  write(*,*) ' LENGTH UNIT : ', unit
  
  ! position initialisation
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
  
  ! distance initialisation 
  dist(:) = 0.d0
  
  ! wave vector initialistaion
  do i=1,Nray 
    read(21,*) k(:,i)
  enddo
  write(*,*) 'k', k

  ! propagation mode
  do i=1,Nray
    read(21,*) mode(i)
  enddo
  write(*,*) 'MODE ',mode(1)


  ! refractive index initialisation
  call varplasma(r,f,X,Y)
  write(*,*) 'X',X(1)
  write(*,*) 'Y',Y(:,1)

  call thetakB(k,r,thkb)
  write(*,*) 'thkB', thkB*180./pi
!  call index(X,Y,thkb,mode,mu)
  call index(X,Y,k,mode,mu)
  write(*,*) 'mu', mu(1)
  
  !lecture de l'indice max
  read(21,*) mumax
  write(*,*) 'mu max', mumax
  
  !normalisation of k by normk = mu
  normk = sqrt(k(1,:)**2+k(2,:)**2+k(3,:)**2)
  forall (i=1:3) k(i,:) = (k(i,:)/normk(:))*mu(:)
  write(*,*) 'k normalise', k

  
  ! time step initialisation
  read(21,*) dtinit
  read(21,*) dtmax
  read(21,*) dtmin  
  write(*,*) 'dt init', dtinit, 'dt max', dtmax, 'dt min', dtmin
  dt(:)=dtinit
  
  ! arbitrary initialisation of pratio in circular polar
  where (mode .eq. 'X')
    pratio(1,:) = -ic
  elsewhere
    pratio(1,:) = ic
  end where
  pratio(2,:) = 0.d0
  pratio(3,:) = 0.d0
  ! arbitrary initialisation LP='VP'
  LP(:)='VP'  
  
  ! compute some initial parameters
  call density(r, Nel)
  call vphase(mu,vp)
  call thetakb(k,r,thkb)
  call fplasma(f,X,fp)
  call fcyclo(f,Y,fc)
  call fuph(f,X,Y,fuh)
  call polar_ratio(mode, f,k,r,LP,pratio,pratio_1)
  
  !WRITING IN VISU FILE
!  open(18, file='./'//path//'/visu_param.txt')
  open(18, file='./visu_param.txt')
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

  !SAVE INITIAL PARAMETERS
  
  9991 format(26(1e17.10,1x),a1,1x,a2,1x,1e17.10,1x,1e17.10)
  if (integ .eq. 'RK4T') then
!    filen='./'//path//'/dataRK4T_' 
    filen='./dataRK4T_' 
  else 
!    filen='./'//path//'/dataRK4F_'
    filen='./dataRK4F_'
  end if
  filen(15:19) = '.dat'
   
  do i=1,nray
    write(filen(12:14),'(i3.3)') i
    !open(21+i,file=filen, form='formatted')
    open(21+i,file=filen)!, form='formatted')
    write(21+i,'(30(a,1x))') 'x (km)','y (km)','z (km)','kx','ky','kz','indice','pas','vg (c)','vp (c)', &
    !write(21+i,*) 'x (km)','y (km)','z (km)','kx','ky','kz','indice','pas','vg (c)','vp (c)', &
      &  '(k.B) (deg)','(r.B) (deg)','distance (km)','Bx (T)','By (T)','Bz (T)', &
      &  'fp (Hz)','fc (Hz)','fUH (Hz)','Ne (cm-3)',   &
      &  'Re(Ey/Ex)','Im(Ey/Ex)','Re(Ez/Ex)','(Ez/Ex)','Re(Ez/Ey)','Im(Ez/Ey)','mode','LP','err','errdisp'

  write(21+i,9991) r(:,i),k(:,i),mu(i),dt(i), -1.d0,vp(i), thkb(i)*180.d0/pi, -1.d0, &
  !write(21+i,*) r(:,i),k(:,i),mu(i),dt(i), -1.d0,vp(i), thkb(i)*180.d0/pi, -1.d0, &
      & 0.d0, B(:,i), fp(i), fc(i), fuh(i),Nel(i),&
      & real(pratio_1(1,i)), aimag(pratio_1(1,i)),&
      & real(pratio_1(2,i)), aimag(pratio_1(2,i)),&
      & real(pratio_1(3,i)), aimag(pratio_1(3,i)),&
      & mode(i), LP(i),0.0,0.0
  close(21+i)
  enddo
 
  !initialisation kiter
  kiter(:) = 1.d0

  !LOOP

  do while (sum(kiter) .ne. Nray*Niter*1.d0)
    if (integ .eq. 'RK4') then
      call RK4(r,k,dt,dxdt,dkdt,f,mode,r_1,k_1)
    else
      call RK4F(r,k,dt,dxdt,dkdt,f,mode,dtmin,dtmax,r_1,k_1)
    end if

    call varplasma(r_1,f,X,Y)
    call index(X,Y,k_1,mode,mu_1)
    !renormalisation of k by mu
    normk = sqrt(k_1(1,:)**2+k_1(2,:)**2+k_1(3,:)**2)
    forall (i=1:3) k_1(i,:) = (k_1(i,:)/normk(:))*mu_1(:)
    
    ! dr=sqrt(r_1^2-r^2)
    dr(:)=sqrt(r_1(1,:)**2+r_1(2,:)**2+r_1(3,:)**2)-sqrt(r(1,:)**2+r(2,:)**2+r(3,:)**2)
    
    !Compute parameters
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

    
    ! time step adaptation
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
    ! exit conditions
    !----------------------------------
    kiter(numray) = kiter(numray) + 1.d0
    flagmu(:)=0
    ! out of bounds
    where(((r(1,:)+dr(:)) .ge. x_axis(size_x)) .or. ((r(1,:)-dr(:)) .le. x_axis(1)) .or. &
         ((r(2,:)+dr(:)) .ge. y_axis(size_y)) .or. ((r(2,:)-dr(:)) .le. y_axis(1)) .or. &
         ((r(3,:)+dr(:)) .ge. z_axis(size_z)) .or. ((r(3,:)-dr(:)) .le. z_axis(1)))
      kiter(numray)=Niter*1.d0
!      write(*,*)'STOP ray out of bound'
      flagmu = numray
    end where 
    if(sum(flagmu) .ne. 0) then
      write(*,*)flagmu,'STOP ray out of bound (out of box)'
!      write(*,*) flagmu, 'STOP MU'
    endif
    
    flagmu(:)=0
    where((r(1,:)*r(1,:) + r(2,:)*r(2,:) + r(3,:)*r(3,:))<= r_planet2)
      kiter(numray)=Niter*1.d0
      flagmu = numray
    end where
    if(sum(flagmu) .ne. 0) then
      write(*,*)flagmu,'STOP ray out of bound (in planet)'
      write(*,*)"x ",r(1,:)
      write(*,*)"y ",r(2,:)
      write(*,*)"z ",r(3,:)
      write(*,*)"radius2 ",(r(1,:)*r(1,:) + r(2,:)*r(2,:) + r(3,:)*r(3,:))
      write(*,*)"r_planet2 ",r_planet2,r_planet

    endif

      


    ! mu > mumax
    flagmu(:)=0
    where (mu .ge. mumax)  
      kiter(numray) = Niter*1.d0 
      flagmu = numray
    end where
    if(sum(flagmu) .ne. 0) then
      write(*,*) flagmu, 'STOP MU'
    endif
  
    !coord. of r = NaN
    flagnan(:)=0
    where (((r(1,:)+1)/(r(1,:)+1) .ne. 1) .or. ((r(2,:)+1)/(r(2,:)+1) .ne. 1) .or. ((r(3,:)+1.d0)/(r(3,:)+1.d0) .ne. 1))
      kiter(numray) = Niter*1.d0
      flagnan = numray
    end where
    if(sum(flagnan) .ne. 0) then
      write(*,*) flagnan, 'STOP NaN'
    endif
    
    ! relative error between (n^2-1) et X(pY-1) > 1% => dt = dtmin
    flagerr(:)=0
    call error(k_1,X,Y,mode,err)
    where (err .gt. 1.d-2)  
      dt_1(numray) = dtmin
    end where
    ! dispersion error > 1% => dt = dtmin
    flagerrdisp(:)=0
    call errordisp(k_1,r_1,X,Y,mu_1,errdisp)
    where (errdisp .gt. 1.d-2)  
      dt_1(numray) = dtmin
    end where
  
    mode_1=mode
    
      
    
    !Selection of the parameters of the continuing rays  
    mask = kiter(numray) .lt. Niter
    Nray3 = count(mask)

    deallocate (r,k,dt,mu,numraytemp, mode, mux, muo,pratio)
    allocate (r(3,Nray3), k(3,Nray3),dt(Nray3),mu(Nray3),numraytemp(Nray3), mode(Nray3), mux(Nray3), muo(Nray3),pratio(3,Nray3))
    r=0 ; k=0 ; dt=0 ; mu=0 ; numraytemp=0 ; mode='' ; mux=0 ; muo=0 ; pratio=0

    j = 1
    do i=1,Nray2
      if (mask(i)) then
        ! writing in file
        write(filen(12:14),'(i3.3)') i
        open(21+i,file=filen,position='append')
        write(21+numray(i),9991) r_1(:,i),k_1(:,i),mu_1(i),dt_1(i),&
          & vg(i), vp(i), thkb(i)*180.d0/pi, thrb(i)*180.d0/pi,&
          & dist(numray(i)), B(:,i), &
          & fp(i), fc(i), fuh(i),Nel(i),&
          & real(pratio_1(1,i)), aimag(pratio_1(1,i)),&
          & real(pratio_1(2,i)), aimag(pratio_1(2,i)),&
          & real(pratio_1(3,i)), aimag(pratio_1(3,i)),&
          & mode_1(i),&! abs(mux_1(i)-muo_1(i)), dmux(i), dmuo(i),&
          & LP(i), err(i), errdisp(i)                
        close(21+i)
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
    r_1=0; k_1=0; Y=0; B=0; pratio_1=0
    dt_1=0; mu_1=0; X=0; normk=0; err=0; errdisp=0; vp=0; vg=0; LP(:)='VP'
    thkb=0; thkr=0; thrb=0; mask=.FALSE.; flagmu=0; flagnan=0; flagerr=0; flagerrdisp=0
    fp=0; fc=0; fuh=0; Nel=0; mode_1=''; mux_1=0; muo_1=0; dmux=0; dmuo=0
    modeX=''; modeO=''; dr=0
  end do ! END LOOP
  
  call cpu_time(t2)
  write(*,*) '----------------------------------'
  write(*,*) 'temps de calcul (s)', t2-t1

!  do i=1,nray
!    close(21+i)
!  enddo

  close(21) ; close(20) ; close(19)
  deallocate(r_1,k_1,dt_1,mu_1,X,Y,normk,err,vp,thkb,thkr,thrb,vg,B,mask, flagmu, flagnan,flagerr,dr)
  deallocate(fp,fc,fuh, kiter,pratio,pratio_1,Nel, mode, mode_1, mux, muo, mux_1, muo_1,dmux,dmuo,modeX,modeO,LP)
  call dealloc_environ()
END PROGRAM raytracing
