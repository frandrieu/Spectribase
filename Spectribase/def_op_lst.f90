!**************************************************************
!**************************************************************
! $Id: def_op_lst.f90,v 1.1 2004/04/07 14:05:35 madec Exp doute $
!
!
! History:
!
! $Log: def_op_lst.f90,v $
! Revision 1.1  2004/04/07 14:05:35  madec
! Initial revision
!
!
!**************************************************************
!**************************************************************
module def_op_lst

  implicit none 

  integer(kind=2), parameter:: nbordm=20

  type wnode
     real(kind=8):: wave
     type(wnode), pointer:: next
  end type wnode

  type indnode
     real(kind=8):: indreal,indimag
     type(indnode), pointer:: next
  end type indnode

  type geo1node
     integer(kind=2):: u,v
     real(kind=4):: muo_a,mu_a,muo_b,mu_b, phi_a
     real(kind=4):: thetauo,thetau,phi,alpha,depth
     integer(kind=2):: p,q
     type(geo1node), pointer:: next
  end type geo1node

  type geo2node
     real(kind=4):: muo,mu,phi,alpha,depth
     integer(kind=2):: p,q
     type(geo2node), pointer:: next
  end type geo2node

  type wpt
    type(wnode), pointer:: head
    type(wnode), pointer:: tail
  end type wpt

  type indpt
     type(indnode), pointer:: head
     type(indnode), pointer:: tail
  end type indpt

  type geo1pt
    type(geo1node), pointer:: head
    type(geo1node), pointer:: tail
  end type geo1pt

  type geo2pt
    type(geo2node), pointer:: head
    type(geo2node), pointer:: tail
  end type geo2pt

contains

  subroutine init_w(head,tail)

    ! initialise une liste vide

    type(wnode), pointer:: head, tail

    nullify(head,tail)

  end subroutine init_w

  subroutine add_w(new,head,tail)
    ! Ajoute un nouveau noeud à la fin de la liste

    type(wnode), pointer:: new, head, tail

    ! Vérification pour savoir si la liste est vide
    if(associated(head)) then
       ! La liste n'est pas vide
       tail%next => new
       nullify(new%next)
       tail => new
    else 
       ! La liste est vide
       head => new
       tail => new
       nullify(tail%next)
    end if
  end subroutine add_w

  subroutine delete_w(head,tail,first)
    ! Retourne un pointeur pointant sur le premier noeud de la liste
    ! chaînée et l'ôte de la liste

    type(wnode), pointer::head,tail,first

    ! Vérification pour savoir si la liste est vide

    if(associated(head)) then
       ! La liste n'est pas vide
       ! Vérification nombre d'éléments dans la liste
       if(associated(head%next)) then
          ! il y a plus d'un élément dans la liste
          first=>head
          head=>head%next
       else 
          ! il y a un seul élément dans la liste
          first=>head
          nullify(head,tail) ! liste vide maintenant
       end if
    else 
       ! La liste est vide
       nullify(first) ! aucun élément retourné
    end if
  end subroutine delete_w

  subroutine init_ind(head,tail)

    ! initialise une liste vide

    type(indnode), pointer:: head, tail

    nullify(head,tail)

  end subroutine init_ind

  subroutine add_ind(new,head,tail)
    ! Ajoute un nouveau noeud à la fin de la liste

    type(indnode), pointer:: new, head, tail

    ! Vérification pour savoir si la liste est vide
    if(associated(head)) then
       ! La liste n'est pas vide
       tail%next => new
       nullify(new%next)
       tail => new
    else 
       ! La liste est vide
       head => new
       tail => new
       nullify(tail%next)
    end if
  end subroutine add_ind

  subroutine delete_ind(head,tail,first)
    ! Retourne un pointeur pointant sur le premier noeud de la liste
    ! chaînée et l'ôte de la liste

    type(indnode), pointer::head,tail,first

    ! Vérification pour savoir si la liste est vide

    if(associated(head)) then
       ! La liste n'est pas vide
       ! Vérification nombre d'éléments dans la liste
       if(associated(head%next)) then
          ! il y a plus d'un élément dans la liste
          first=>head
          head=>head%next
       else 
          ! il y a un seul élément dans la liste
          first=>head
          nullify(head,tail) ! liste vide maintenant
       end if
    else 
       ! La liste est vide
       nullify(first) ! aucun élément retourné
    end if
  end subroutine delete_ind

  subroutine init_geo1(head,tail)

    ! initialise une liste vide

    type(geo1node), pointer:: head, tail

    nullify(head,tail)

  end subroutine init_geo1

  subroutine add_geo1(new,head,tail)
    ! Ajoute un nouveau noeud à la fin de la liste

    type(geo1node), pointer:: new, head, tail

    ! Vérification pour savoir si la liste est vide
    if(associated(head)) then
       ! La liste n'est pas vide
       tail%next => new
       nullify(new%next)
       tail => new
    else 
       ! La liste est vide
       head => new
       tail => new
       nullify(tail%next)
    end if
  end subroutine add_geo1

  subroutine delete_geo1(head,tail,first)
    ! Retourne un pointeur pointant sur le premier noeud de la liste
    ! chaînée et l'ôte de la liste

    type(geo1node), pointer::head,tail,first

    ! Vérification pour savoir si la liste est vide

    if(associated(head)) then
       ! La liste n'est pas vide
       ! Vérification nombre d'éléments dans la liste
       if(associated(head%next)) then
          ! il y a plus d'un élément dans la liste
          first=>head
          head=>head%next
       else 
          ! il y a un seul élément dans la liste
          first=>head
          nullify(head,tail) ! liste vide maintenant
       end if
    else 
       ! La liste est vide
       nullify(first) ! aucun élément retourné
    end if
  end subroutine delete_geo1

  subroutine init_geo2(head,tail)

    ! initialise une liste vide

    type(geo2node), pointer:: head, tail

    nullify(head,tail)

  end subroutine init_geo2

  subroutine add_geo2(new,head,tail)
    ! Ajoute un nouveau noeud à la fin de la liste

    type(geo2node), pointer:: new, head, tail

    ! Vérification pour savoir si la liste est vide
    if(associated(head)) then
       ! La liste n'est pas vide
       tail%next => new
       nullify(new%next)
       tail => new
    else 
       ! La liste est vide
       head => new
       tail => new
       nullify(tail%next)
    end if
  end subroutine add_geo2

  subroutine delete_geo2(head,tail,first)
    ! Retourne un pointeur pointant sur le premier noeud de la liste
    ! chaînée et l'ôte de la liste

    type(geo2node), pointer::head,tail,first

    ! Vérification pour savoir si la liste est vide

    if(associated(head)) then
       ! La liste n'est pas vide
       ! Vérification nombre d'éléments dans la liste
       if(associated(head%next)) then
          ! il y a plus d'un élément dans la liste
          first=>head
          head=>head%next
       else 
          ! il y a un seul élément dans la liste
          first=>head
          nullify(head,tail) ! liste vide maintenant
       end if
    else 
       ! La liste est vide
       nullify(first) ! aucun élément retourné
    end if
  end subroutine delete_geo2

end module def_op_lst
