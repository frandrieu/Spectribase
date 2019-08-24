C       **************************************************************
C	$Id: calc_texprop.f,v 1.1 2004/04/07 14:05:35 madec Exp doute $
C	**************************************************************
C	**************************************************************
C	**                                                          **
C	**  Calculation of the single scattering albedo and the     **
C	**  optical depth of a granular texture layer : intimate    **
C	**  mixture of different grain types                        **
C	**                                                          **
C	**************************************************************
C	**************************************************************

C	**************************************************************
C	History:
C
C	$Log: calc_texprop.f,v $
C	Revision 1.1  2004/04/07 14:05:35  madec
C	Revised F.Andrieu 2014: Adding routines for compact slab structures
C
C
C	**************************************************************

C	
	subroutine proptexg(nb,n,k,prop,diam,s,compa,long,rdepth,
     &  singalbg,taug)
C
C	|------------------------------------------------------------|
C	|  Parameters :                                              |
C	|  - in   nb number of grain types                           |
C	|         n vector of real indexes                           |
C	|         k vector of imaginary indexes                      |
C	|         prop vector of proportions (according to density)  |
C	|         diam vector of grain sizes (micron)                |
C	|         s vector of the internal diffusion parameters of   |
C	|           the grains (micron-1)                            |
C       |         compa compactness of the current layer             |
C	|         long wavelength (micron)                           |
C       |         rdepth real metric thickness of the current        |
C       |           layer (mm)                                       |
C	|                                                            |
C	|  - out  singalbg single scattering albedo                  |
C       |         taug optical depth of the current layer            |
C	|                                                            |
C	|  Other quantities :                                        |
C       |         diamp grain size parameter of Hapke model (micron) |
C       |         sfact grain size factor                            |
C	|         sects normalized mean scattering cross section     |
C	|         secte normalized mean geometrical cross section    |
C	|         qs scattering efficiency for an equant grain       |
C	|         si internal reflection coefficient                 |
C	|         se external reflection coefficient                 |
C	|         teta mean transmission factor                      |
C	|         iota absorption factor                            |
C	|                                                            |
C	|------------------------------------------------------------|
C	
	implicit none
	integer*2 nb,i
	double precision se,si,iota,ri,teta,qs,sects
	double precision long,singalbg,x,y,pi,rdepth
	double precision compa,taug,dens,diamp,sfact
	double precision secte,prop(nb),diam(nb),n(nb),k(nb),s(nb)
	parameter(pi=3.14159,sfact=0.5)

C	**************************************************************

	secte=0.
	sects=0.
	dens=0.
	
	do 10 i=1,nb
	
C	*** calculation of the scattering efficiency for the type i	
	   iota=4.*pi/long*k(i)
	   if ((iota.eq.0.).and.(s(i).eq.0.)) then
	      x=0.
	   else
	      x=dsqrt(iota/(iota+s(i)))
	   endif
	   y=dsqrt(iota*(iota+s(i)))
	   ri=(1-x)/(1+x)
	   diamp=diam(i)*sfact
	   teta=(ri+dexp(-y*diamp))/
     &     (1+ri*dexp(-y*diamp))

	   call coefreflg(n(i),k(i),se,si)
	   qs=se+((1-se)*(1-si)*teta)/(1-si*teta)
	   
C	*** calculation of the mean scattering and geometrical
C	cross section	   
	   sects=sects+prop(i)*pi*diam(i)**2*qs/4.
	   secte=secte+prop(i)*pi*diam(i)**2/4.
	   
C	*** calculation of the total density	   
	   dens=dens+diam(i)**3*prop(i)
	   
10	continue

	dens=6.*compa/pi/dens	   
         
        taug=-dens*dlog(1-compa)/compa*secte*rdepth*1.d03
         
	singalbg=sects/secte
	!singalbg=1.
	return
	end
C
C	**************************************************************

C	**************************************************************
C	*                                                            *
C	*  Calculation of the external and internal reflection       *
C	*  coefficients for an irregular equant particle             *
C	*                                                            *
C	**************************************************************
C	
	subroutine coefreflg(n,k,se,si)
	
C	|------------------------------------------------------------|
C	|                                                            |
C	|  Parameters :                                              |
C	|  - in   (n,k) complex index of the material                |
C	|                                                            |
C	|  - out  se external reflection coefficient                 |
C	|         si internal reflection coefficient                 |
C	|                                                            |
C	|  Other quantities :                                        |
C	|         rpe perpendicular Fresnel coefficient              |
C	|         rpa parallel Fresnel coefficient                   |
C	|         m1 complex index of the first medium               |
C	|         m2 complex index of the second medium              |
C	|         zlegen vector of the Legendre Polynomial roots     |
C	|         a vector of the associated Christoffel numbers     |
C	|                                                            |
C	|------------------------------------------------------------|
	
	implicit none
	integer*2 i
	double precision se,si,rpe,rpa,zlegen(32),a(32),x,y,pi,mul
        double precision n, k, mu
	complex m1,m2,beta1,beta2
	parameter(pi=3.14159)
	common/legendre16/zlegen,a
C
C   	**************************************************************
C	
	se=0.
	si=0.
	
	do 10 i=1,32
C	*** calculation of the Fresnel coefficients for a
C	void-material interface with the incidence arcos(zlegen(i))	
C       Integration using Gauss quadrature approximation mu=zlegen(i)
	mu=0.5*zlegen(i)+0.5
	   m1=(1.,0.)
	   m2=cmplx(n,k)
	   beta1=m1*mu
	   beta2=csqrt(m2**2-m1**2*(1-sngl(mu)**2))
	   if (aimag(beta2).lt.0) beta2=-beta2	      
	   rpe=cabs((beta1-beta2)/(beta1+beta2))
	   rpa=cabs((beta1/m1**2-beta2/m2**2)/
     &     (beta1/m1**2+beta2/m2**2))
     
	   se=se+a(i)*mu*(rpe**2+rpa**2)    !Integration using Gauss quadrature approximation
	   
C	*** calculation of the Fresnel coefficients for a
C	material-void interface with the incidence arcos(mu) 
C       Integration using Gauss quadrature approximation mu=zlegen(i)	
C	(Hapke 2012 sec 4.3.4 eq 4.37-4.38 pp55-58)
	   m1=cmplx(sngl(n),sngl(k))
	   m2=(1.,0.)
	   mul=cos(asin(real(m2)/real(m1)))
	   if (mu.gt.mul) then
	      x=real(m1)/real(m2)*sqrt(1-mu**2)         
	      y=sqrt(1-x**2)
	      beta1=csqrt(m1**2-m2**2*real(x)**2)
	      if (aimag(beta1).lt.0) beta1=-beta1
	      beta2=m2*y
	      rpe=cabs((beta1-beta2)/(beta1+beta2))
	      rpa=cabs((beta1/m1**2-beta2/m2**2)/
     &        (beta1/m1**2+beta2/m2**2))
	   else
C	* total reflection	   
	      rpe=1.
	      rpa=1.
	   endif
	   
	   si=si+a(i)*mu*(rpe**2+rpa**2)*0.5   !Integration using Gauss quadrature approximation
	      
10	continue

	return
	end
	
	subroutine coefreflcp(n,k,se,si)
	
C	|------------------------------------------------------------|
C	|  Creation :             06/2014, F. Andrieu                |
C	|         calculation of external and internal reflexion     |
C	|         coefficients se and si for a flat interface under  |
C	|         an isotropic radiation.                            |
C	|  Parameters :                                              |
C	|  - in   (n,k) complex index of the material                |
C	|                                                            |
C	|  - out  se external reflection coefficient                 |
C	|         si internal reflection coefficient                 |
C	|                                                            |
C	|  Other quantities :                                        |
C	|         rpe perpendicular Fresnel coefficient              |
C	|         rpa parallel Fresnel coefficient                   |
C	|         m1 complex index of the first medium               |
C	|         m2 complex index of the second medium              |
C	|         zlegen vector of the Legendre Polynomial roots     |
C	|         a vector of the associated Christoffel numbers     |
C	|                                                            |
C	|------------------------------------------------------------|
	
	implicit none
	integer*2 i
	double precision se,si,rpe,rpa,zlegen(32),a(32),x,y,pi,mul
        double precision n, k, mu
	complex m1,m2,beta1,beta2
	parameter(pi=3.14159)
	common/legendre16/zlegen,a
C
C   	**************************************************************
C	
	se=0.
	si=0.
	
	do 10 i=1,32
C	*** calculation of the Fresnel coefficients for a
C	void-material interface with the incidence arcos(zlegen(i))	
C       Integration using Gauss quadrature approximation mu=0.5*zlegen(i)+0.5
	!mu=0.5*zlegen(i)+0.5
        mu=cos(pi/4.*zlegen(i)+pi/4.)
	   m1=(1.,0.)
	   m2=cmplx(n,k)
	   beta1=m1*mu
	   beta2=csqrt(m2**2-m1**2*(1-sngl(mu)**2))
	   if (aimag(beta2).lt.0) beta2=-beta2	      
	   rpe=cabs((beta1-beta2)/(beta1+beta2))
	   rpa=cabs((beta1/m1**2-beta2/m2**2)/
     &     (beta1/m1**2+beta2/m2**2))
     
	   se=se+a(i)*mu*(rpe**2+rpa**2)*pi/4.!/dsqrt(1.-(mu)**2)    !Integration using Gauss quadrature approximation
	   
C	*** calculation of the Fresnel coefficients for a
C	material-void interface with the incidence arcos(zlegen(i)) 
C       Integration using Gauss quadrature approximation mu=zlegen(i)	
C	(Hapke 2012 sec 4.3.4 eq 4.37-4.38 pp55-58)
	   m1=cmplx(sngl(n),sngl(k))
	   m2=(1.,0.)
	   mul=cos(asin(real(m2)/real(m1)))
	   if (mu.gt.mul) then
	      x=real(m1)/real(m2)*sqrt(1-mu**2)         
	      y=sqrt(1-x**2)
	      beta1=csqrt(m1**2-m2**2*real(x)**2)
	      if (aimag(beta1).lt.0) beta1=-beta1
	      beta2=m2*y
	      rpe=cabs((beta1-beta2)/(beta1+beta2))
	      rpa=cabs((beta1/m1**2-beta2/m2**2)/
     &        (beta1/m1**2+beta2/m2**2))
	   else
C	* total reflection	   
	      rpe=1.
	      rpa=1.
	   endif
	   
	   si=si+a(i)*mu*(rpe**2+rpa**2)*pi/4.!/dsqrt(1.-(mu)**2)   !Integration using Gauss quadrature approximation
	      
10	continue

	return
	end
	
C	**************************************************************

C	**************************************************************
C	**************************************************************
C	**                                                          **
C	**  Calculation of the single scattering albedo and the     **
C	**  optical depth of a compact texture layer with an        **
C	**  intimate mixture of different inclusion types           **
C	**                                                          **
C	**************************************************************
C	**************************************************************
C	
	subroutine proptexc(nb,n,k,prop,diam,s,compa,long,rdepth,
     &  singalbc,tauc, R0, T0)
C
C	|------------------------------------------------------------|
C	|  Parameters :                                              |
C	|  - in   nb number of inclusion types                       |
C	|         n vector of real indexes                           |
C	|         k vector of imaginary indexes                      |
C	|         prop vector of proportions (according to density)  |
C	|         diam vector of inclusion sizes (micron)            |
C	|         s vector of the internal diffusion parameters of   |
C	|           the inclusions (micron-1)                        |
C	|         compa compactness of the ice matrix                |
C	|         long wavelength (micron)                           |
C       |         rdepth real metric thickness of the current        |
C       |           layer (mm)                                       |
C	|                                                            |
C	|  - out  singalbc single scattering albedo                  |
C       |         tauc optical depth of the current layer            |
C       |         R0 diffuse isotropic reflectance of the slab       |
C       |         T0 diffuse isotropic transmittance of the slab     |
C	|                                                            |
C	|  Other quantities :                                        |
C       |         diamp inclusion size parameter of the modified     |
C       |           Hapke model (micron)                             |
C       |         sfact inclusion size factor                        |
C	|         sects mean scattering cross section                |
C	|         secte mean geometrical cross section               |
C	|         dens total density of the inclusions               |
C	|         qs scattering efficiency for an equant inclusion   |
C	|         si internal reflection coefficient                 |
C	|         se external reflection coefficient                 |
C	|         teta mean transmission factor                      |
C	|         iota absorption factor                             |
C       |         stot scattering coefficient of the layer           |
C       |         etot attenuation coefficient of the layer          |
C	|                                                            |
C	|  Modification hystory : 06/2014, F. Andrieu                |
C	|         calling new subroutine coefreflcp instead of       |
C	|         coefreflg for the calculation of se and si for the |
C	|         slab                                               |
C	|------------------------------------------------------------|
C	
	implicit none
	integer*2 nb,i
	double precision se,si,iota,teta,qs,sects
	double precision long,singalbc,x,y,pi,coefreflc
	double precision secte,prop,diam,n,k,s,compa,dens
	double precision tauc,rdepth,diamp,sfact
	double precision R0, T0, ri, factexp1, factexp2, sep, stot, etot
	dimension prop(nb),diam(nb),n(0:nb),k(0:nb),s(nb)
	parameter(pi=3.14159,sfact=0.5)

C	**************************************************************

C	*** initialisation
	secte=0.
	sects=0.
	dens=0.
	
	do 10 i=1,nb
	
C	*** calculation of the scattering efficiency for the type i	
	   iota=4.*pi/long*k(i)
	   if ((iota.eq.0.).and.(s(i).eq.0.)) then
	      x=0.
	   else
	      x=dsqrt(iota/(iota+s(i)))
	   endif
	   y=dsqrt(iota*(iota+s(i)))
	   ri=(1.-x)/(1.+x)
	   diamp=diam(i)*sfact
	   teta=(ri+dexp(-y*diamp))/
     &     (1.+ri*dexp(-y*diamp))

           if((n(0).eq.n(i)).and.(k(0).eq.k(i))) then
C       matrix and inclusions are identical
              se=0.
              si=0.
           else
C       matrix and inclusions are different
    	      se=coefreflc(1,n(0),k(0),n(i),k(i),diam(i),long)
	      si=coefreflc(2,n(i),k(i),n(0),k(0),diam(i),long)
           endif
           
C	Equivalent slab approximation (Hapke eq 5.52a pp. 95-97)
           qs=se+((1-se)*(1-si)*teta)/(1-si*teta)
	   	   
C	*** calculation of the mean scattering and geometrical
C	cross section	   
	   sects=sects+prop(i)*pi*diam(i)**2*qs/4.
	   secte=secte+prop(i)*pi*diam(i)**2/4.
	   
C	*** calculation of the total density	   
	   dens=dens+diam(i)**3*prop(i)
	   
10	continue

        iota=4.*pi/long*k(0)
C	*** total density eq 2.2-23 p. 63 these S. Douté
	dens=6.*(1-compa)/pi/dens

        stot=dens*sects

        etot=dens*secte

        singalbc=stot
     &  /(etot+0.5*(1.+compa)*iota)
        if(singalbc.gt.0.9999) singalbc=0.9999

        tauc=(etot*(1.+0.5*(1.-compa))+iota)
     &  *rdepth*1.d3
 

C       *** calculation of the reflection and transmission functions
C           for a rough specularly reflecting matrix
C        case 1 : isotropic illumination

	call coefreflcp(n(0),k(0),se,si)

        ri=(1-dsqrt(1-singalbc))/(1+dsqrt(1-singalbc))

        y=dsqrt(1-singalbc)*tauc

        factexp1=dexp(-4.*y)
        factexp2=dexp(-2.*y)

        R0=(ri*(1-si*ri)+(si-ri)*factexp1)/
     &  ((1-si*ri)**2-(si-ri)**2*factexp1)*
     &  (1-se)*(1-si)+se

        T0=((1-si*ri)*factexp2+ri*(si-ri)*factexp2)/
     &  ((1-si*ri)**2-(si-ri)**2*factexp1)*
     &  (1-se)*(1-si)

	return
	end
C
C	**************************************************************

C	**************************************************************
C	*                                                            *
C	*  Calculation of the external or internal reflection        *
C	*  coefficients for an irregular equant inclusion            *
C       *  (version using Hapke formula to calculate Fresnel         *
C       *  coefficients)                                             *
C	*                                                            *
C	**************************************************************
C

	function coefreflc(indc,n1,k1,n2,k2,diam,long)

C
C	|------------------------------------------------------------|
C	|  Parameters :                                              |
C	|  - in   indc entrance indicator (se  indc=1 ; si  indc=2 ) |
C	|         n1c real index of the first medium                 |
C	|         k1c complex index of the first medium              |
C	|         n2c real index of the second medium                |
C	|         k2c complex index of the second medium             |
C	|         diam inclusion size (micron)                       |
C	|         long wavelength (micron)                           |
C	|                                                            |
C	|  - out  coefreflc reflection coefficient                   |
C	|                                                            |
C	|  Other quantities :                                        |
C	|         rpe perpendicular Fresnel coefficient              |
C	|         rpa parallel Fresnel coefficient                   |
C	|         zlegen vector of the Legendre Polynomial roots     |
C	|         a vector of the associated Christoffel numbers     |
C	|         eps threshold for considering a material as non    |
C       |         absorbent                                          |
C	|                                                            |
C	|  Modification hystory : 06/2014, F. Andrieu                |
C	|         correction of the expression of Hapke parallel and |
C	|         perpandicular reflexion coefficients               |
C	|                                                            |
C	|------------------------------------------------------------|

	implicit none
	integer*2 i,indc
	double precision n1,n2,k1,k2,n,k,mu
        double precision beta1, beta2, g1, g2
	double precision rpe,rpa
	double precision a(32),zlegen(32),long,pi,diam,coefreflc
        double precision  eps
        parameter(eps=1.d-15)
	common/legendre16/zlegen,a

C	*************************************************************
C	*** initialisation (Hapke 2012 p47 eq 4.5-4.6)
        
        n=(n1*n2+k1*k2)/(n1**2+k1**2)
        k=(n1*k2-n2*k1)/(n1**2+k1**2)
	

	coefreflc=0.
C	
	do 10 i=1,32
	
	mu=0.5*zlegen(i)+0.5
	
C	*** calculation of the Fresnel coefficients for an interface
C	between two non absorbent materials with the incidence
C	arcos(mu)	(Hapke 2012 sec 4.3.4 eq 4.37-4.38 pp55-58)
	if((k1.le.eps).and.(k2.le.eps)) then
	   if(n2.gt.(n1*dsqrt(1-mu**2))) then
	      beta1=n1*mu
	      beta2=dsqrt(n2**2-n1**2*dsqrt(1-mu**2)**2)
	      rpe=abs((beta1-beta2)/(beta1+beta2))
	      rpa=abs((beta1/n1**2-beta2/n2**2)/
     &        (beta1/n1**2+beta2/n2**2))
           else
C	* total reflection           
              rpe=1.
              rpa=1.
           endif 

        else
           
C	*** calculation of the Fresnel coefficients for an interface
C	between two materials with a least one absorbent material (Hapke 4.44-4.49 (2012) ou 4.33-4.38 (1993))
C PREVIOUS VERSION
C           g1=0.5*(n**2-k**2-dsqrt(1-mu**2)+
C     &     dsqrt((n**2-k**2-dsqrt(1-mu**2))**2+4.*n**2*k**2)) 
C           g2=0.5*(-(n**2-k**2-dsqrt(1-mu**2))+
C     &     dsqrt((n**2-k**2-dsqrt(1-mu**2))**2+4.*n**2*k**2))
C           rpe=((mu-g1)**2+g2**2)/((mu+g1)**2+g2**2)
C           rpe=dsqrt(rpe)
C           rpa=(((n**2-k**2)*mu-g1)**2+(2.*n*k*mu-g2)**2)/
C     &     (((n**2-k**2)*mu+g1)**2+(2.*n*k*mu-g2)**2)
C           rpa=dsqrt(rpa)

C CORRECTED VERSION
           g1=dsqrt(0.5*(n**2-k**2-(1-mu**2)+
     &     dsqrt((n**2-k**2-(1-mu**2))**2+4.*n**2*k**2))) 
           g2=dsqrt(0.5*(-(n**2-k**2-(1-mu**2))+
     &     dsqrt((n**2-k**2-(1-mu**2))**2+4.*n**2*k**2)))
           rpe=((mu-g1)**2+g2**2)/((mu+g1)**2+g2**2)
           rpe=dsqrt(rpe)
           rpa=(((n**2-k**2)*mu-g1)**2+(2.*n*k*mu-g2)**2)/
     &     (((n**2-k**2)*mu+g1)**2+(2.*n*k*mu-g2)**2)
           rpa=dsqrt(rpa)

   
     	endif
     	
    	if(indc.eq.1) then
C	* if se is calculated, a differential absorption factor  
C	is taken into account
     	   coefreflc=coefreflc+dexp(-2.*pi/long*k1*diam*(1-mu))
     &     *a(i)*mu*(rpe**2+rpa**2)*0.5	 
     	else
     	   coefreflc=coefreflc+a(i)*mu*(rpe**2+rpa**2)*0.5    !Integration using Gauss quadrature approximation
    	endif
    	
10	continue

	return
	end
C	
C	**************************************************************
