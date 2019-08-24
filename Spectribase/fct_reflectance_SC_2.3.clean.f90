!MODIFICATION HISTORY : 12/12/2013
!Modification of est_fact : 
!	tabrap : cos(v)*sin(v)/|...|
!	/cos(dzeta*pi/180.)**2 suppressed in est_fact value
!
!Modification of fct_relectance_SC_spec  : 
!	replaced proba_spec=2./(pi*sin(dzeta)**2)*(1+tgv2)*exp(-tgv2/(pi*tan(dzeta)**2)) by
!	 	 proba_spec=2./(pi*tan(dzeta)**2)*sqrt(tgv2)*sqrt(1+tgv2)*exp(-tgv2/(pi*tan(dzeta)**2)) 
!	removed /domega in reflectance_SC_spec value
!	modified the calculation on sin(v)dpsi by only dpsi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function fct_reflectance_SC_diff(R0,T0,R0p,T0p,w)

  ! |------------------------------------------------------------|
  ! |  Parameters :                                              |
  ! |  - in   R0 coefficient for diffuse reflection in the case  |
  ! |          of an isotropic illumination                      |
  ! |         T0 coefficient for diffuse transmission in the case|
  ! |          of an isotropic illumination                      |
  ! |         R0p coefficient for diffuse reflection in the case |
  ! |          of an collimated illumination                     |
  ! |         T0p coefficient for diffuse transmission in the    |
  ! |          case of an collimated illumination                |
  ! |         w single scattering albedo                         |
  ! |                                                            |
  ! |  - out  fct_reflectance_SC_diff                            |
  ! |         isotropic contribution of the reflectance of a slab|
  ! |         of ice superimposed on a granular substrat         |
  ! |                                                            |
  ! |  Other quantities :                                        |
  ! |         rs isotropic reflectance of the substrat           |
  ! |                                                            |
  ! |------------------------------------------------------------|

  implicit none

  ! interface
  real:: R0,T0,R0p,T0p,w
  real:: fct_reflectance_SC_diff

  ! local variables
  real:: rs, gama

  gama=sqrt(1.-w)

  rs=(1.-gama)/(1.+gama)
  !print *, R0, T0, rs, R0p, T0p
  fct_reflectance_SC_diff=R0p+T0p*T0*rs/(1.-R0*rs)

  return

end function fct_reflectance_SC_diff

  ! **************************************************************
  ! **************************************************************

subroutine slab_refl_trans(long, n, k, singalbc, tauc, rdepth, fact, sep, R0p, T0p)

  ! |------------------------------------------------------------|
  ! |  Parameters :                                              |
  ! |  - in   long wavelength (micron)                           |
  ! |         n real index                                       |
  ! |         k imaginary index                                  |
  ! |         singalbc single scattering of the slab             |
  ! |         tauc optical depth of the slab                     |
  ! |         rdepth real metric thickness of the current        |
  ! |           layer (mm)                                       |
  ! |         fact mean "airmass" of the slab                    |
  ! |         sep external reflection coefficient for collimated |
  ! |           radiation                                        |
  ! |                                                            |
  ! |  - out  R0p coefficient for diffuse reflection in the case |
  ! |          of an collimated illumination                     |
  ! |         T0p coefficient for diffuse transmission in the    |
  ! |          case of an collimated illumination                |
  ! |                                                            |
  ! |  Other quantities :                                        |
  ! |         si internal reflection coefficient                 |
  ! |         se external reflection coefficient for isotropic   |
  ! |           radiation                                        |
  ! |         teta mean transmission factor for isotropic        |
  ! |           radiation                                        |
  ! |         tetap mean transmission factor for collimated      |
  ! |           radiation                                        |
  ! |         iota absorption factor                             |
  ! |                                                            |
  ! |------------------------------------------------------------|

  implicit none

  !parameters
  real, parameter:: pi=3.141592653589793

  ! interface
  real:: long, n, k, singalbc, tauc, rdepth, R0p, T0p

  !local variables
  real:: iota, x, y, ri, fact, teta, tetap, sep
  real:: est_fact, coefreflgp, factexp1, factexp2  
  real(kind=8)::  se, si
  !real,parameter:: eps=1.e-9
  !       *** calculation of the reflection and transmission functions
  !           for a rough specularly reflecting matrix
  !        case 2 : collimated illumination

  ri=(1-sqrt(1-singalbc))/(1+sqrt(1-singalbc))
  
  y=sqrt(1-singalbc)*tauc

  teta=(ri+exp(-2.*y))/(1+ri*exp(-2.*y))

  factexp1=exp(-fact*y)
  factexp2=exp(-2.*y)

  !tetap=(ri+exp(-fact*y))/(1+ri*exp(-fact*y))

  call coefreflcp(real(n,8), real(k,8) , se, si) !Attention c'était coeflreflg
!  if( (1.-teta*si).lt.eps )then
!    R0p=(1.-sep)*si/(1.+si)!+sep
!    T0p=(1.-sep)/(1.+si)
!  else
!    T0p=(1-sep)*tetap*(1-si)/(1-(teta*si)**2)
!    R0p=si*tetap*teta*(1-sep)*(1-si)/(1-(teta*si)**2)!+sep
!  endif

!if(abs(R0p+T0p+sep-1).gt.1.e-4) then
!  print *, R0p,T0p,sep, singalbc, tauc, se, si 
!endif

  R0p=(ri*(1-si*ri)+(si-ri)*factexp1*factexp2)/ &
  &  ((1-si*ri)**2-(si-ri)**2*factexp2**2)* &
  &  (1-sep)*(1-si)

  T0p=((1-si*ri)*factexp1+ri*(si-ri)*factexp1)/ &
  &  ((1-si*ri)**2-(si-ri)**2*factexp2**2)* &
  &  (1-sep)*(1-si)
!  T0p=1-R0p-sep
  return

end subroutine slab_refl_trans

  ! **************************************************************
  ! **************************************************************

function est_fact(n,theta0,dzeta)

  ! |------------------------------------------------------------|
  ! |  Parameters :                                              |
  ! |  - in   n real index                                       |
  ! |         theta0 incidence angle (in degree)                 |
  ! |         dzeta mean roughness slope (in degree)             |
  ! |                                                            |
  ! |  - out  est_fact mean "airmass" of the slab                |
  ! |                                                            |  
  ! |  Other quantities :                                        |
  ! |         v array of facet tilt angles                       |
  ! |         xhi array of facet azimuth angles                  |
  ! |                                                            |
  ! |------------------------------------------------------------|

  implicit none

  ! interface
  real:: n,theta0,dzeta,est_fact 

  ! parameters
  real, parameter:: pi=3.141592653589793
  integer, parameter:: nb_max=1000 ! number of samples for estimation

  ! local variables :
  REAL, DIMENSION(nb_max) :: u1, u2, xhi, v, tabcosi, tabrap
  real:: muo
  if(dzeta.gt.5.) dzeta=5.
  muo=cos(theta0*pi/180.)

  ! generate nb_max random numbers uniformally distributed in (0,1(
  CALL RANDOM_SEED
  CALL RANDOM_NUMBER(u1)

  ! generate nb_max random numbers uniformally distributed in (0,1(
  CALL RANDOM_SEED
  CALL RANDOM_NUMBER(u2)

  ! generate nb_max random facet orientations distributed according to an
  ! importance sampling distribution funtion

  where(u1/=0.)
     v=atan(sqrt(-pi*tan(dzeta*pi/180.)**2*log(u1)))
  elsewhere
     v=pi/2.
  end where

  xhi=2*pi*u2

  ! calculation of the corresponding table of cosines of factor incidence angle
  tabcosi=sqrt(1-muo**2)*sin(v)*cos(xhi)+cos(v)*muo

  tabrap=cos(v)*sin(v)/abs(-1./n*muo+cos(v)*(1./n*tabcosi-sqrt(1-1./n**2*(1-tabcosi**2))))

  ! integration
  est_fact=1./nb_max*sum(tabrap)!/cos(dzeta*pi/180.)**2

  return

end function est_fact

  ! **************************************************************
  ! **************************************************************

!	**************************************************************
!	*                                                            *
!	*  Calculation of the external and internal reflection       *
!	*  coefficients for an irregular equant particle             *
!	*                                                            *
!	**************************************************************
!
function coefreflgp(n, k, theta0, dzeta)

!n, k, theta0, dzeta

  ! |------------------------------------------------------------|
  ! |  Parameters :                                              |
  ! |  - in   n real index                                       |
  ! |         k imaginary index                                  |
  ! |         theta0 incidence angle (in degree)                 |
  ! |         dzeta mean roughness slope (in degree)             |
  ! |                                                            |
  ! |  - out  coefreflgp external reflection coefficient for     |
  ! |         collimated radiation                               |
  ! |                                                            |  
  ! |  Other quantities :                                        |
  ! |         theta emission angle (in degree)                   |
  ! |         alpha phase angle (in degree)                      |
  ! |         nphong exponent of the phong specular model        |
  ! |         beta array of zenith angles for the emergence      |
  ! |           direction in the specular coordinate system      |
  ! |         fi azimuth of the current emergence direction      |
  ! |           in the specular coordinate system                |
  ! |                                                            |
  ! |------------------------------------------------------------|

  implicit none

  ! interface
  real::  n, k, theta0, dzeta

  ! parameters
  integer, parameter:: nb_max=50, nfi=20
  real, parameter:: pi=3.141592653589793
  real, parameter:: eps=5.e-4
  ! local variables :
  INTEGER :: i, j, nb_fi, nerr
  real, DIMENSION(nb_max) :: u1, beta
  real, DIMENSION(:), ALLOCATABLE :: u2, fi
  real:: muo, mu, x1, y1, z1, theta, alpha, c_phi, s_phi
  real:: prodt, filim
  real:: est, coefreflgp, phi, detjac, reflectance_SC_spec
  real:: dbeta, dfi, dx1beta, dx1fi, dy1beta, dy1fi, dg1beta, dg1fi
  real:: dz1beta, dz1fi, dg2beta, dg2fi, detjacbetafi, tgv2, proba_spec
  real:: mu_a, muo_a, mu_b, muo_b, rf, rpe, rpa, fct_S
  complex:: m1, m2, beta1, beta2
  real, parameter::sample=0.8
nerr=0
  ! cosine of the incidence angle
  muo=cos(theta0*pi/180.)
  !if(dzeta.gt.5.) dzeta=5.

  ! exponent of the specular Phong model
  if(dzeta > 45.) dzeta=45. ! protection
!  nphong=log(0.5)/log(cos(2.*pi/180.*dzeta))

  u1=(/(i,i=1,nb_max,1)/)

  if(dzeta > 5.) then
    u1=real(u1)/nb_max 
    ! zenith angle with respect to the specular direction 
    !UNIFORM SAMPLING OF BETA
    beta=(pi/2.+theta0*pi/180.-2.*eps)*u1+eps
  else
    !PREFERENTIAL SAMPLING OF BETA ('sample' percent in a 10*dzeta interval)
    do i=1, int(sample*nb_max) 
      beta(i)=eps+real(i)*(10.*dzeta*pi/180.)/(sample*real(nb_max))
    enddo
    do i=int(sample*nb_max)+1,nb_max
      beta(i)=eps+10.*dzeta*pi/180.+real(i-int(sample*nb_max))*&
          &(theta0*pi/180.+pi/2.-2*eps-10.*dzeta*pi/180.)/((1.-sample)*nb_max)
    enddo
  endif
  ! generate nb_max random emergence directions around the specular direction 
  ! eliminating those that are below the local horizon




  ! resetting the estimator
  est=0.

  do i=1, nb_max
   
     nb_fi=min(int(abs(sin(beta(i))*nfi/2.-1))+10,int(nfi/2.))

     if((theta0==0.).or.(beta(i)==0.).or.(beta(i)==pi/2.)) then
        prodt=2.
     else
        prodt=1./(tan(theta0*pi/180.)*tan(beta(i)))
     endif

     filim=0.
     if( prodt .le. 1.) then
     filim=acos(prodt)
     nb_fi=int((pi-filim)/pi*nb_fi)
     endif
  
     allocate(u2(2*nb_fi-1))
     allocate(fi(2*nb_fi-1))
  
     u2=(/(j,j=1,2*nb_fi-1,1)/)
     u2=real(u2)/(2.*nb_fi-1)


     ! emergence always above the horizon
     if( prodt > 1) then
        fi=(2.*pi-eps)*u2+pi ! azimuth angle in the specular coordinate system
        fi=mod(fi, (2.*pi) )
        dfi=(2.*pi-eps)/(2.*nb_fi-1.)
        !UNIFORM SAMPLING OF FI
     else
        ! emergence sometimes below the horizon
        fi=(2.*(pi-filim)-2.*eps)*u2+filim+eps ! azimuth angle in the specular coordinate system
        dfi=(2.*pi-2.*filim-2.*eps)/(2.*nb_fi-1.)
        !UNIFORM SAMPLING OF FI
     endif


     if(dzeta > 5.) then
       !UNIFORM SAMPLING OF BETA	
       dbeta=(pi/2.+theta0*pi/180.-2.*eps)/real(nb_max) !rad
     else
       !PREFERENTIAL SAMPLING OF BETA ('sample' percent in a 10*dzeta interval)
       if(i.le.int(sample*nb_max)) then
         dbeta=10.*dzeta*pi/180./(sample*real(nb_max))
       else
         dbeta=(theta0*pi/180.+pi/2.-10.*dzeta*pi/180.-2*eps)/((1.-sample)*nb_max)
       endif
     endif

    do  j=1, 2*nb_fi-1
      ! direction of emergence in the local coordinate system
      x1=-muo*sin(beta(i))*cos(fi(j))-sqrt(1-muo**2)*cos(beta(i))
      y1=-sin(beta(i))*sin(fi(j))
      z1=-sqrt(1-muo**2)*sin(beta(i))*cos(fi(j))+muo*cos(beta(i))   !z1=cos(theta_emergence) in the local coordinate (non spec) system
      theta=acos(z1)*180./pi ! in deg : emergence angle

     ! calculation of the phase angle in degree
      mu=z1
      c_phi=x1/sqrt(x1**2+y1**2)
      s_phi=y1/sqrt(x1**2+y1**2)
      if(c_phi.gt.1.) c_phi=1.
      if(c_phi.lt.-1.) c_phi=-1.
      alpha=acos(muo*mu+sqrt(1.-muo**2)*sqrt(1.-mu**2)*c_phi)*180./pi
      phi=acos(c_phi)
      if(s_phi.lt.0.) phi=2.*pi-phi
 
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Partial derivatives of x1, xy, z1
      dx1beta=-muo*cos(beta(i))*cos(fi(j))+sqrt(1-muo**2)*sin(beta(i))
      dx1fi=muo*sin(beta(i))*sin(fi(j))
      dy1beta=-cos(beta(i))*sin(fi(j))
      dy1fi=-sin(beta(i))*cos(fi(j))
      dz1beta=-sqrt(1.-muo**2)*cos(beta(i))*cos(fi(j))-muo*sin(beta(i))
      dz1fi=sqrt(1.-muo**2)*sin(beta(i))*sin(fi(j))

      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Terms of the jacobian for the change beta, fi --> v, ksi
      ! using g1(beta, fi)=v, g2(beta, fi)=ksi
      !dg'N''X' is the partial derivative of gN by the var X
      dg1beta=1./((x1+sqrt(1.-muo**2))**2+y1**2+(z1+muo)**2)*((z1+muo)/sqrt((x1+sqrt(1.-muo**2))**2&
         &+y1**2)*((x1+sqrt(1.-muo**2))*dx1beta+y1*dy1beta)-sqrt((x1+sqrt(1.-muo**2))**2+y1**2)*dz1beta)
      dg1fi=1./((x1+sqrt(1.-muo**2))**2+y1**2+(z1+muo)**2)*((z1+muo)/sqrt((x1+sqrt(1.-muo**2))**2+y1**2)*&
         &((x1+sqrt(1.-muo**2))*dx1fi+y1*dy1fi)-sqrt((x1+sqrt(1.-muo**2))**2+y1**2)*dz1fi)
      dg2beta=((x1+sqrt(1.-muo**2))*dy1beta-y1*dx1beta)/((x1+sqrt(1.-muo**2))**2+y1**2)
      dg2fi=((x1+sqrt(1.-muo**2))*dy1fi-y1*dx1fi)/((x1+sqrt(1.-muo**2))**2+y1**2)
      detjacbetafi=abs(dg1beta*dg2fi-dg1fi*dg2beta)
      !!!!!!!!!!!!!!!!!!!!!!!!!

      tgv2=((x1+sqrt(1.-muo**2))**2+y1**2)/(z1+muo)**2
      proba_spec=2./(pi*tan(dzeta*pi/180.)**2)*sqrt(tgv2)*sqrt(1+tgv2)*exp(-tgv2/(pi*tan(dzeta*pi/180.)**2)) 

      m1=cmplx(1.,0.)
      m2=cmplx(n,k)
      beta1=m1*cos(alpha*pi/360.)
      beta2=csqrt(m2**2-m1**2*(1-cos(alpha*pi/360.)**2))
      if (aimag(beta2).lt.0) beta2=-beta2	      
      rpe=cabs((beta1-beta2)/(beta1+beta2))
      rpa=cabs((beta1/m1**2-beta2/m2**2)/&
          & (beta1/m1**2+beta2/m2**2))

      rf=rpe**2+rpa**2 

      reflectance_SC_spec=cos(alpha*pi/360.)/(2.*muo)&
         & *rf*proba_spec*dbeta*dfi*detjacbetafi
     
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Correction of shadowing effect
      if (dzeta/=0.) then
        call correct(theta0, theta, phi, dzeta,& 
           &muo_a, mu_a, muo_b, mu_b)
        reflectance_SC_spec=reflectance_SC_spec*&
           &fct_S(theta0, theta, muo_a, mu_a, muo_b, mu_b, phi, dzeta)
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!


     if((theta0.lt.0.).or.(theta0.ge.90.)) then  
        reflectance_SC_spec=0.
        nerr=nerr+1
     endif

     if((theta.lt.0.).or.(theta.ge.90.)) then
        reflectance_SC_spec=0.
        nerr=nerr+1
     endif

     if((alpha.lt.0.).or.(alpha.gt.180.)) then
        reflectance_SC_spec=0.
        nerr=nerr+1
     endif

    !!!!!!!!!!!!!!!!!!!!!!!!!
    ! Integration
    if ((.not. isnan(reflectance_SC_spec)).and.(reflectance_SC_spec.le.huge(reflectance_SC_spec))) est=est+reflectance_SC_spec
    if (( isnan(reflectance_SC_spec)).or.(reflectance_SC_spec.gt.huge(reflectance_SC_spec))) nerr=nerr+1
    !!!!!!!!!!!!!!!!!!!!!!!!!

    enddo
    deallocate(u2)
    deallocate(fi)
  enddo

  coefreflgp=est/pi
  if( nerr.ne.0 ) print *, nerr, "error(s) occured in estimation of sep for n, k, theta0=", n, k, theta0
  !print *, coefreflgp
  !if( coefreflgp.gt.1.) print *, 'Sep error : Sep = ', coefreflgp , 'Theta0= ', theta0
  return  


end function coefreflgp

  ! **************************************************************
  ! **************************************************************


function fct_reflectance_SC_spec(n, k, dzeta, theta0, theta, phi, alpha, err_flag)

  ! |------------------------------------------------------------|
  ! |  Parameters :                                              |
  ! |  - in   n real index                                       |
  ! |         k imaginary index                                  |
  ! |         dzeta mean roughness slope (in degree)             |
  ! |         theta0 incidence angle (in degree)                 |
  ! |         theta emission angle (in degree)                   |
  ! |         alpha phase angle (in degree)                      |
  ! |                                                            |
  ! |  - out  fct_reflectance_SC_spec                            |
  ! |         specular contribution of the reflectance of a slab |
  ! |         of ice superimposed on a granular substrat         |
  ! |                                                            |
  ! |  Other quantities :                                        |
  ! |         dtheta0 inf. quantity linked with the sun apparent |
  ! |           angle                                            |
  ! |         dtheta inf. quantity linked with the instrument    |
  ! |           solid angle of acceptance                        |
  ! |         dphi inf. quantity linked with the sun apparent    |
  ! |           angle and with the instrument solid angle of     |
  ! |           acceptance                                       |
  ! |         mu cosinus of theta                                |
  ! |         muo cosinus of theta0                              |
  ! |         muo_b factor                                       |
  ! |         mu_b factor                                        |
  ! |         c_phi  cosinus of the azimuth angle                |
  ! |         phi azimuth angle (in radian)                      |
  ! |         v tilt of a facet to be in specular conditions     |
  ! |         psi azimuth of a facet to be in specular conditions|
  ! |         proba_spec probability for a facet to be in        |
  ! |         orientation (v, psi)                               |
  ! |         rf Fresnel reflection factor for the given geometry|
  ! |                                                            |
  ! |------------------------------------------------------------|


  implicit none

  ! interface
  real:: n,k,theta0, theta, alpha, dzeta, phi
  integer:: err_flag
  real:: fct_reflectance_SC_spec
  integer:: nb_max
  ! parameters
  ! apparent diameter of the sun from Mars (0°21')
  real, parameter:: dtheta0=0.0061086523819802 ! (rad)
  ! apparent diameter of the source from IPAG spectro-goniometer (0,2°)
  !real, parameter:: dtheta0=0.00349066 ! (rad)
  ! apparent diameter of the source for simulation purposes (0,001°)
  !real, parameter:: dtheta0=0.!0000175 ! (rad)
  ! field of view of an OMEGA pixel (0°4.15')
  real, parameter:: dtheta=0.0012071860659627 ! (rad)
  ! field of view of a CRISM pixel (0.00174°)
  !real, parameter:: dtheta=3.03e-5 ! (rad)
  ! field of view of the IPAG spectro-goniometer (4.1° diameter) 
  !real, parameter:: dtheta=0.0715585 ! (rad)
  ! field of view for simulation purposes (1° diameter) 
  !real, parameter:: dtheta=0.0174533 ! (rad)
  ! field of view for simulation purposes (180° diameter) 
  !real, parameter:: dtheta=3.141592653589793
  ! field of view for simulation purposes (0.25° diameter) 
  !real, parameter:: dtheta=0.00436332 ! (rad)
  ! solid angle sustained by an OMEGA pixel 
  !real, parameter:: domega=1.4572981978546e-06 ! (rad**2)
  ! solid angle sustained by an CRISM pixel 
  !real, parameter:: domega=2.9e-09 ! (sr)
  !real, parameter:: dphi=0.0073158360659627 ! (sum of the 2 previous values)
  !real, parameter:: dphi=0.00613895 ! (sum of the 2 previous values)
  real, parameter:: pi=3.141592653589793
  real, parameter:: eps=1.e-9
  ! local variables :
  INTEGER :: i, j, nerr, nb_h, nb_wd, integ
  real, DIMENSION(:), ALLOCATABLE  :: u1,u2,u3,u4
  real, DIMENSION(:), ALLOCATABLE  :: x1,x2,x3,x4,x5
  real:: c_phi, s_phi, phil, corjac, ener_rep
  real:: prodt, dphi, phi2, theta2,x, dphii
  real:: est, coefreflgp, detjac, reflectance_SC_spec
  real:: dx1e, dx1psi, dg1e, dg1psi
  real:: dg2e,dg2psi,detjacepsi,tgv2,tgv22,tgv23,tgv24,tgv25
  real:: mu,muo,mu_a,muo_a,mu_b,muo_b, rf, rpe, rpa, fct_S
  real:: a1,a2,a3,c,b1,b2,b3, elli, elli1, elli2, elli3, elli4
  real:: proba_spec, c_alpha, omegac, omegaci, ellimax, ellimin
  real:: theta_C, theta0_c, dzeta_c, dtheta0p, dthetap
  complex:: m1, m2, beta1, beta2

  if(dzeta.gt.5.) dzeta=5.
 
  dthetap=dtheta+dtheta0
  ! cosine of the incidence angle
  muo=cos(theta0*pi/180.)
  if(abs(muo-1.).lt.eps) muo=1.-(theta0*pi/180.)**2/2.
  ! cosine of the emergence angle
  mu=cos(theta*pi/180.)
  if(abs(mu-1.).lt.eps) mu=1.-(theta*pi/180.)**2/2.

  ! if far from sepecular, avoid calculations
  x=(sin(theta*pi/180.))**2+1.-muo**2+2.*( sin(theta*pi/180.)*sqrt(1.-muo**2)*cos(phi) )
  if(x/((mu+muo)**2).gt.7.*dzeta*pi/180.+dthetap) then
     fct_reflectance_SC_spec=0.
     return
  endif 
  ener_rep=1.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! For spectro-imaging simulations
  integ=0
  !Captor's solid angle
  omegac=2.*pi*(1.-cos(dtheta/2.))
  !omegac=2.*pi*(1.-cos(dthetaip/2.))
  !omegac=dtheta**2!*pi/4.
  !calculating the field of view in emergence and azimuth
  if(theta*pi/180..lt.dtheta*0.5) then
    dphi=2.*pi
  else
    dphii=atan(tan(dtheta/2.)/sin(theta*pi/180.))*2.
    dphi=atan(tan((dtheta+dtheta0)/2.)/sin(theta*pi/180.))*2.
    if ( isnan(dphii)) dphii=pi
    if ( isnan(dphi)) dphi=pi
  endif
  omegaci=pi*dtheta*dphii/4.
  !if( omegaci.gt.2.*pi*dtheta) omegaci=2.*pi*dtheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! WARNING : ONLY FOR ENERGY 
!!CONSERVATION TESTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  integ=1
!  dphi=dtheta
!  omegac=sin(theta*pi/180.)*dtheta**2
!  if( sin(theta*pi/180.).eq.0.) omegac=pi*theta*dtheta**2/4.
!  if( omegac.eq.0.) omegac=eps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtheta0p=dtheta0*2./3.
  
  !if( dzeta*pi/180..le.dtheta0/2.) dtheta0p=2.*dzeta*pi/180. !to avoid "missing" the peak of a(v)
  ! **************************************************************
  !determination of the number of integration points inside the pixel
  nb_max=(int(2.*dthetap/(dzeta*pi/180.))+5)
  nb_h=int(real(nb_max)/2.)
  nb_max=2*nb_h+1 !making nb_max uneven to make sure to have the central point
  if(nb_max.lt.3) nb_max=3
  
  
  !Allocating space for the differents variables arrays
  allocate(u1(nb_max**2))
  allocate(u2(nb_max**2))
  allocate(u3(nb_max**2))
  allocate(u4(nb_max**2))
  allocate(x1(nb_max**2))
  allocate(x2(nb_max**2))
  allocate(x3(nb_max**2))
  allocate(x4(nb_max**2))
  allocate(x5(nb_max**2))

  !Creation of the gridding inside the pixel
  do i=1, nb_max
    do j=1,nb_max
      theta2=abs( theta*pi/180.+real(2*i-nb_max-1)*dtheta/(2.*nb_max) )
      if(theta2.gt.pi/2.) then
        u1((i-1)*nb_max+j)=pi/2.
      else
        u1((i-1)*nb_max+j)=theta2
      endif
      phi2=abs(phi+real(2*j-nb_max-1)*dphii/(2.*nb_max))
      if(phi2.gt.pi) then
        u2((i-1)*nb_max+j)=2.*pi-phi2
      else
        u2((i-1)*nb_max+j)=phi2
      endif 
      if( (u2((i-1)*nb_max+j)+dtheta0p/(2.*nb_max)).gt.pi ) then
        u3((i-1)*nb_max+j)=abs(2.*pi-(u2((i-1)*nb_max+j)+dtheta0p/(2.*nb_max)))
      else
        u3((i-1)*nb_max+j)=abs(u2((i-1)*nb_max+j)+dtheta0p/(2.*nb_max))
      endif
      u4((i-1)*nb_max+j)=abs(u2((i-1)*nb_max+j)-dtheta0p/(2.*nb_max))
    enddo
  enddo
  
  ! resetting the estimator
  est=0.

  if((theta0.lt.0.).or.(theta0.ge.90.)) then 
     fct_reflectance_SC_spec=0.
     err_flag=1
     return
  endif

  if((theta.lt.0.).or.(theta.ge.90.)) then
     fct_reflectance_SC_spec=0.
     err_flag=1
     return
  endif

  if((alpha.lt.0.).or.(alpha.gt.180.)) then
     fct_reflectance_SC_spec=0.
     err_flag=1
     return
  endif


  !terms that will be useful for future calculations
  x1=(sin(u1))**2+1.-muo**2+2.*( sin(u1)*sqrt(1.-muo**2)*cos(u2) )
  if( dtheta0p.gt.0.) then
    x2=(sin(u1))**2+(sin(theta0*pi/180.+dtheta0p/2.))**2+2.*( sin(u1)*sin(theta0*pi/180.+dtheta0p/2.)*cos(u2) )
    x3=(sin(u1))**2+1.-muo**2+2.*( sin(u1)*sqrt(1.-muo**2)*cos(u3) )
    x4=(sin(u1))**2+(sin(theta0*pi/180.-dtheta0p/2.))**2+2.*( sin(u1)*sin(theta0*pi/180.-dtheta0p/2.)*cos(u2) )
    x5=(sin(u1))**2+1.-muo**2+2.*( sin(u1)*sqrt(1.-muo**2)*cos(u4) )
  endif

  nerr=0
  if(dzeta > 45.) dzeta=45. ! protection
  fct_reflectance_SC_spec=0.
  nb_wd=0
!! Integration
  do i=1,nb_max**2    
    elli1=( 2.*(theta*pi/180.-u1(i)+dtheta/2./nb_max)/dtheta )**2+( 2.*(phi-u2(i)+dphii/2./nb_max)/dphii )**2
    elli2=( 2.*(theta*pi/180.-u1(i)+dtheta/2./nb_max)/dtheta )**2+( 2.*(phi-u2(i)-dphii/2./nb_max)/dphii )**2
    elli3=( 2.*(theta*pi/180.-u1(i)-dtheta/2./nb_max)/dtheta )**2+( 2.*(phi-u2(i)+dphii/2./nb_max)/dphii )**2
    elli4=( 2.*(theta*pi/180.-u1(i)-dtheta/2./nb_max)/dtheta )**2+( 2.*(phi-u2(i)-dphii/2./nb_max)/dphii )**2
    ellimax=max(elli1, elli2, elli3, elli4) ! if all 4 corners inside ellipse ok
    ellimin=min(elli1, elli2, elli3, elli4) ! if 1 corner inside ellipse ok
    elli=( 2.*(theta*pi/180.-u1(i))/dtheta )**2+( 2.*(phi-u2(i))/dphii )**2 !if center inside ellipse ok
    if(elli.ge.1..and.integ.eq.0) then
      nb_wd=nb_wd+1
    else
      if(ellimax.gt.1) then
          ener_rep=0.5
      endif
      mu=cos(u1(i))
      c_phi=cos(u2(i))
      s_phi=sin(u2(i))
      if(c_phi.gt.1.) c_phi=1.
      if(c_phi.lt.-1.) c_phi=-1.
      if(c_phi.eq.-1.) c_phi=-1.+(pi-u2(i))**2/2.
      if(mu.gt.1.) mu=1.
      if(abs(mu-1.).lt.eps) mu=1.-(u1(i))**2/2.

     
      c_alpha=muo*mu+sqrt(1.-muo**2)*sqrt(1.-mu**2)*c_phi
      if(c_alpha.gt.1.) c_alpha=1.
      alpha=acos(c_alpha)*180./pi
      if( x1(i).lt.0.) x1(i)=0.
      if( dtheta0p.gt.0.) then
        if( x2(i).lt.0.) x2(i)=0.
        if( x3(i).lt.0.) x3(i)=0.
        if( x4(i).lt.0.) x4(i)=0.
        if( x5(i).lt.0.) x5(i)=0.
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Partial derivatives of x1
      ! +
      ! Terms of the jacobian for the change e, psi --> v, ksi
      ! using g1(e, psi)=v, g2(e, psi)=ksi
      !dg'N''X' is the partial derivative of gN by the var X
      !!!!!!!!!!!!!!!!!!!!!!!!!

      if(abs(muo-1.).lt.eps) then 
        detjacepsi=0.5
      else
        if(abs(mu-1.).lt.eps) then
          detjacepsi=0.
        else
          if(muo.eq.mu) then
            if(c_phi.eq.-1.) then
              detjacepsi=0.5
            else
              if(c_phi.eq.1.) then
              detjacepsi=0.5
              else
                dx1e=2.*mu*( sqrt(1.-mu**2)+sqrt(1.-muo**2)*c_phi )
                dx1psi=-2.*sqrt(1.-mu**2)*sqrt(1.-muo**2)*s_phi
                dg1e=sqrt(x1(i))/( x1(i)+(mu+muo)**2 )*( (mu+muo)*dx1e/(2.*x1(i)) + sqrt(1.-mu**2) )
                dg1psi=(mu+muo)/( x1(i)+(mu+muo)**2 )*dx1psi/(2.*sqrt(x1(i)))
                dg2e=-1./sqrt(( 1.-c_phi**2 )*(1.-mu**2))*( mu*c_phi-&
                        & ( sqrt(1.-muo**2)+sqrt(1.-mu**2)*c_phi )*dx1e/(2.*x1(i)) )
                dg2psi=1./sqrt(( 1.-c_phi**2 )*(1.-mu**2))*( sqrt(1.-mu**2)*s_phi+&
                        & ( sqrt(1.-muo**2)+sqrt(1.-mu**2)*c_phi )*dx1psi/(2.*x1(i)) )
                detjacepsi=abs(dg1e*dg2psi-dg1psi*dg2e)
              endif  
            endif
          else       
            if(c_phi.eq.-1.) then
              dg1e=( cos(2.*u1(i))+cos(u1(i)-theta0*pi/180.) )/( 2.*( 1.+cos(u1(i)+theta0*pi/180.) ) )
              if(muo.ge.mu) dg1e=-dg1e
              dg1psi=0.
              dg2e=0.
              dg2psi=1.-sqrt(1.-muo**2)/(sqrt(1.-muo**2)-sqrt(1.-mu**2))        
              detjacepsi=abs(dg1e*dg2psi-dg1psi*dg2e)
            else
              if(c_phi.eq.1.) then
                detjacepsi=0.5
              else
                dx1e=2.*mu*( sqrt(1.-mu**2)+sqrt(1.-muo**2)*c_phi )
                dx1psi=-2.*sqrt(1.-mu**2)*sqrt(1.-muo**2)*s_phi
                dg1e=sqrt(x1(i))/( x1(i)+(mu+muo)**2 )*( (mu+muo)*dx1e/(2.*x1(i)) + sqrt(1.-mu**2) )
                dg1psi=(mu+muo)/( x1(i)+(mu+muo)**2 )*dx1psi/(2.*sqrt(x1(i)))
                dg2e=-1./sqrt(( 1.-c_phi**2 )*(1.-mu**2))*( mu*c_phi-&
                         & ( sqrt(1.-muo**2)+sqrt(1.-mu**2)*c_phi )*dx1e/(2.*x1(i)) )
                dg2psi=1./sqrt(( 1.-c_phi**2 )*(1.-mu**2))*( sqrt(1.-mu**2)*s_phi+&
                         & ( sqrt(1.-muo**2)+sqrt(1.-mu**2)*c_phi )*dx1psi/(2.*x1(i)) )
                detjacepsi=abs(dg1e*dg2psi-dg1psi*dg2e)
              endif
            endif
          endif
        endif
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      if( dtheta0p.eq.0.) then 
        tgv2=x1(i)/(mu+muo)**2
        proba_spec=1./(pi**2*tan(dzeta*pi/180.)**2)*sqrt(tgv2)*sqrt(1+tgv2)*exp(-tgv2/(pi*tan(dzeta*pi/180.)**2))
      else
        !Calculation of (tan(v)**2 where v is the slope satisfying the specular conditions
        tgv2=x1(i)/(mu+muo)**2
        tgv22=x2(i)/(mu+cos(theta0*pi/180.+dtheta0p/2.))**2
        tgv23=x3(i)/(mu+muo)**2
        tgv24=x4(i)/(mu+cos(theta0*pi/180.-dtheta0p/2.))**2
        tgv25=x5(i)/(mu+muo)**2

        !Value of the probability density funcion of the distribution o slopes (Hapke) a(v)     
        proba_spec=1./(pi**2*tan(dzeta*pi/180.)**2)*( sqrt(tgv2)*sqrt(1+tgv2)*exp(-tgv2/(pi*tan(dzeta*pi/180.)**2))+&
             & sqrt(tgv22)*sqrt(1+tgv22)*exp(-tgv22/(pi*tan(dzeta*pi/180.)**2))+&
             & sqrt(tgv23)*sqrt(1+tgv23)*exp(-tgv23/(pi*tan(dzeta*pi/180.)**2))+&
             & sqrt(tgv24)*sqrt(1+tgv24)*exp(-tgv24/(pi*tan(dzeta*pi/180.)**2))+&
             & sqrt(tgv25)*sqrt(1+tgv25)*exp(-tgv25/(pi*tan(dzeta*pi/180.)**2)) )/5.

      endif
      ! Reflexion coefficients
      m1=cmplx(1.,0.)
      m2=cmplx(n,k)
      beta1=m1*cos(alpha*pi/360.)
      beta2=csqrt(m2**2-m1**2*(1-cos(alpha*pi/360.)**2))
      if (aimag(beta2).lt.0) beta2=-beta2	      
      rpe=cabs((beta1-beta2)/(beta1+beta2))
      rpa=cabs((beta1/m1**2-beta2/m2**2)/&
          & (beta1/m1**2+beta2/m2**2))

      rf=rpe**2+rpa**2 

      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Integration
    

      reflectance_SC_spec=cos(alpha*pi/360.)*rf*proba_spec*detjacepsi*ener_rep

      
      if ((.not. isnan(reflectance_SC_spec)).and.(reflectance_SC_spec.le.huge(reflectance_SC_spec))) est=est+reflectance_SC_spec
      if (( isnan(reflectance_SC_spec)).or.(reflectance_SC_spec.gt.huge(reflectance_SC_spec))) nerr=nerr+1
      !!!!!!!!!!!!!!!!!!!!!!!!!
    endif
  enddo
  if ( isnan(est)) est=0.
  fct_reflectance_SC_spec=pi/(omegac*muo*mu)*est*omegaci/(nb_max**2-nb_wd)!omegac=captor's solid angle
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! Correction of shadowing effect
  if (dzeta/=0.) then
    call correct(theta0, theta, phi, dzeta,& 
       &muo_a, mu_a, muo_b, mu_b)
   fct_reflectance_SC_spec=fct_reflectance_SC_spec*&
       &fct_S(theta0, theta, muo_a, mu_a, muo_b, mu_b, phi, dzeta)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!


  return  

end function fct_reflectance_SC_spec

  ! **************************************************************
  ! **************************************************************



