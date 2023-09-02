MODULE static   

IMPLICIT NONE

REAL(8), PARAMETER  :: DIEL_0=8.85418d-14 ! F/cm

REAL(8), PARAMETER  :: DIEL_METAL=0.0D0

REAL(8), PARAMETER  :: m0=5.6856D-16 !eV s^2 / cm^2   rest mass of electron
real(8), PARAMETER  :: m0ISO=9.10938D-28   ! gram   rest mass of electron
real(8), PARAMETER  :: gram2eV=m0/m0ISO    ! eV s^2 / cm^2 / gram
REAL(8), PARAMETER  :: hbar=6.58211899E-16 !eV s

REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
REAL(8), PARAMETER  :: PII=3.14159265d0
REAL(8), PARAMETER  :: ELCH=1.60217653d-19   ! C
real(8), PARAMETER  :: v_light=2.998d10      ! cm/s

REAL(8), PARAMETER :: htqm=(hbar)**2/m0		 !/ELCH
real(8), parameter :: hb2m=7.6305d-16     	 ! eV*cm^2

complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)

END MODULE static
