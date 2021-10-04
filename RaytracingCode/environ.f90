MODULE environ
	
	use constantes
	implicit none
	
	contains

subroutine magneticf(V,B)
	!----------------------------------------------------------
	! Calcul du champ magnetique
	! CALLING SEQUENCE: call magneticf(V,B)
	! INPUTS:  V: vecteur position en coordonnées cart.
	! OUTPUTS: B: vecteur champ magnetique en coord cart. en T
	!----------------------------------------------------------
	real(kind=8), dimension(:,:),intent(in)  :: V
	real(kind=8), dimension(:,:),intent(out) :: B

	B(1,:)=0. 
	B(2,:)=0.
	B(3,:)=B0*(V(3,:)/z0B)**(-3)

end subroutine magneticf


subroutine read_environ()
        !----------------------------------------------------
        ! Lit le fichier init_environ.txt et range les donnees
        !  dans un fichier
        !----------------------------------------------------
        use constantes

        open(40,file='init_environ.txt')
        read(40,*)n0
        read(40,*)z0
        read(40,*)B0
        read(40,*)z0B
        close(40)
        write (*,*)"n0 ",n0," z0 ",z0," B0 ",B0," z0B ",z0B


end subroutine read_environ

subroutine density(V,Ne)
	!-----------------------------------------------------
	! Calcul de la densite d'electron
	! CALLING SEQUENCE: call density(V,Ne)
	! INPUTS:   V: vecteur position en coordonnées cart.
	! OUTPUTS: Ne: densite d'electron en cm-3 
	!-----------------------------------------------------
        use constantes
        real(kind=8),dimension(:,:),intent(in)    :: V
	real(kind=8), dimension(:), intent(out)   :: Ne



!	real(kind=8)                              :: n0,z0

!	open(40, file='init_environ.txt')
!	read(40,*) n0
!	read(40,*) z0

	Ne(:)=n0*8.0*(V(3,:)/z0)**(-2)
	
!	close(40)



end subroutine density

END MODULE environ
