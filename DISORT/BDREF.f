      REAL FUNCTION  BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI )

c     Supplies surface bi-directional reflectivity.
c     
c     NOTE 1: Bidirectional reflectivity in DISORT is defined
c     by Eq. 39 in STWL.
c     NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c     angles) are positive.
c     
c     INPUT:
c     
c     WVNMLO : Lower wavenumber (inv cm) of spectral interval
c     
c     WVNMHI : Upper wavenumber (inv cm) of spectral interval
c     
c     MU     : Cosine of angle of reflection (positive)
c     
c     MUP    : Cosine of angle of incidence (positive)
c     
c     DPHI   : Difference of azimuth angles of incidence and reflection
c     (radians)
c     
c     LOCAL VARIABLES:
c     
c     B0     : empirical factor to account for the finite size of
c     particles in Hapke's BDR model
c     
c     B      : term that accounts for the opposition effect
c     (retroreflectance, hot spot) in Hapke's BDR model
c     
c     CTHETA : cosine of phase angle in Hapke's BDR model
c     
c     GAMMA  : albedo factor in Hapke's BDR model
c     
c     H0     : H( mu0 ) in Hapke's BDR model
c     
c     H      : H( mu ) in Hapke's BDR model
c     
c     HH     : angular width parameter of opposition effect in Hapke's
c     BDR model
c     
c     P      : scattering phase function in Hapke's BDR model
c     
c     THETA  : phase angle (radians); the angle between incidence and
c     reflection directions in Hapke's BDR model
c     
c     W      : single scattering albedo in Hapke's BDR model
c     
c     
c     Called by- DREF, SURFAC
c     +-------------------------------------------------------------------+
c     .. Scalar Arguments ..

      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
c     ..
c     .. Local Scalars ..

      REAL      B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W, i, e
      REAL      MUP_B, MU_B
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, CORRECT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC COS, SQRT
c     ..

c     ** Hapke's BRDF model (times Pi/Mu0)
c     ** (Hapke, B., Theory of reflectance
c     ** and emittance spectroscopy, Cambridge
c     ** University Press, 1993, Eq. 8.89 on
c     ** page 233. Parameters are from
c     ** Fig. 8.15 on page 231, expect for w.)

      W = 0.69
      R = 11.0
      BB = 0.241
      CC = 0.478
      HH   = 0.085
      B0   = 1.0
      R = 3.14159 * R / 180.

      i = ACOS(MUP)
      e = ACOS(MU)

c     Little trick to make the function work properly!
c     DPHI = ACOS(MU);

c$$$      CTHETA = MU*MUP+SIN(ACOS(MU))*SIN(ACOS(MUP))*COS(DPHI)
c$$$      THETA = ACOS( CTHETA )
      CTHETA = MU*MUP+SQRT(1.-MU**2)*SQRT(1.-MUP**2)*COS( DPHI )
      THETA = ACOS(CTHETA)

      CALL CORRECT( i, e, ABS(DPHI), R, MUP, MU, MUP_B, MU_B )

      P=(1.-CC)*(1.-BB**2)/((1.+2.*BB*CTHETA+BB**2)**(3./2.))
      P=P+CC*(1.-BB**2)/((1.-2.*BB*CTHETA+BB**2)**(3./2.))

      B    = B0 * HH / ( HH + TAN( THETA/2.) )

      GAMMA = SQRT( 1. - W )
      H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
      H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

      BDREF = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

c     S function not working yet!
      BDREF = BDREF * S( i, e, MUP, MU, MUP_B, MU_B, abs(DPHI), R ) 
     & * MUP / COS(i)

      MUP = COS(i)
      MU = COS(e)

      RETURN
      END

c     ------------------

      SUBROUTINE  CORRECT( i, e, DPHI, R, MUP, MU, MUP_B, MU_B )

      REAL      i, e, DPHI, R, MUP, MU, MUP_B, MU_B

      REAL      xidz, pi

      pi = 3.14159

      xidz = 1. / SQRT( 1. + pi * TAN(R)**2 )

      IF( i.LE.e ) THEN

         MUP=xidz*(COS(i)+SIN(i)*TAN(R)*
     &        (COS(DPHI)*E2(R,e)+SIN(DPHI/2.)**2*E2(R,i))
     &        /(2.-E1(R,e)-(DPHI/pi)*E1(R,i)))

         MU=xidz*(COS(e)+SIN(e)*TAN(R)*
     &        (E2(R,e)-SIN(DPHI/2.)**2*E2(R,i))
     &        /(2.-E1(R,e)-(DPHI/pi)*E1(R,i)))

      ELSE

         MUP=xidz*(COS(i)+SIN(i)*TAN(R)*
     &        (E2(R,i)-SIN(DPHI/2.)**2*E2(R,e))
     &        /(2.-E1(R,i)-(DPHI/pi)*E1(R,e)))

         MU=xidz*(COS(e)+SIN(e)*TAN(R)*
     &        (COS(DPHI)*E2(R,i)+SIN(DPHI/2.)**2*E2(R,e))
     &        /(2.-E1(R,i)-(DPHI/pi)*E1(R,e)))

      END IF

      MU_B=xidz*(COS(e)+SIN(e)*TAN(R)*(E2(R,e))/(2.-E1(R,e)))
      MUP_B=xidz*(COS(i)+SIN(i)*TAN(R)*(E2(R,i))/(2.-E1(R,i)))


      END

c     ------------------

      REAL FUNCTION  E1( R, X )

      REAL      R, X

      E1 = EXP( -2./3.14159 * cot(R) * cot(X) )

      RETURN
      END

c     ------------------

      REAL FUNCTION  E2( R, X )

      REAL      R, X

      E2 = EXP( -1./3.14159 * cot(R)**2 * cot(X)**2 )

      RETURN
      END

c     ------------------

      REAL FUNCTION  cot( X )

      REAL      X

      IF( X.EQ.0. ) THEN 
         cot = 10.**16 
      ELSE 
         cot = COS(X) / SIN(X)
      END IF

      RETURN
      END

c     ------------------

      REAL FUNCTION S( i, e, MUP, MU, MUP_B, MU_B, DPHI, R )

      REAL  i, e, MUP, MU, MUP_B, MU_B, DPHI, R 

      REAL  f, xidz, pi

      pi = 3.14159

      IF( DPHI.GT.3 ) THEN
         f=0.
      ELSE
         f=EXP(-2.*TAN(DPHI/2.))
      END IF

      xidz = 1. / SQRT(1. + pi * TAN(R)**2)

      IF( i.LE.e ) THEN
         S=MU_B*MUP_B*(1.-f+f*xidz*COS(i)/MUP_B)
         S=MU*COS(i)*xidz/S
      ELSE
         S=MU_B*MUP_B*(1.-f+f*xidz*COS(e)/MU_B)
         S=MU*COS(i)*xidz/S
      END IF

      RETURN
      END
