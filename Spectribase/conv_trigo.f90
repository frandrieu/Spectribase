!**************************************************************
!**************************************************************
! $Id: conv_trigo.f90,v 1.1 2004/04/07 14:05:35 madec Exp doute $
!
!
! History:
!
! $Log: conv_trigo.f90,v $
! Revision 1.1  2004/04/07 14:05:35  madec
! Initial revision
!
!
!**************************************************************
!**************************************************************

module conv_trigo

implicit none 

contains

  ! conversion deg->rad
  real(kind=4) function rad(angl)

    implicit none

    real(kind=4), parameter:: pi=3.141592653589793
    real(kind=4), intent(in):: angl

    rad=angl/180.*pi 

  end function rad

  ! conversion rad->deg
  real(kind=4) function deg(angl)

    implicit none

    real(kind=4), parameter:: pi=3.141592653589793
    real(kind=4), intent(in):: angl

    deg=angl*180./pi 

  end function deg

  ! produit vectoriel de deux vecteurs	

  subroutine prodvect(vect1,vect2,vect3)

    real(kind=4), dimension(3):: vect1,vect2,vect3

    vect3(1)=vect1(2)*vect2(3)-vect2(2)*vect1(3)		
    vect3(2)=vect1(3)*vect2(1)-vect2(3)*vect1(1)
    vect3(3)=vect1(1)*vect2(2)-vect2(1)*vect1(2)

    return

  end subroutine prodvect

end module conv_trigo
