MODULE CHI2GAUSS
USE CONSTRAINTS
USE SPLINES
USE COVAR
IMPLICIT NONE

INTEGER gaussnmom, gaussndata, icalls, irun
REAL*8 :: gaussdomega, chi2alpha, chi2val, chi2reg
REAL*8, ALLOCATABLE :: gaussref(:), xref(:), yref(:), y2ref(:), curspec(:)
REAL*8, ALLOCATABLE :: gauss(:), dgauss(:), fgauss(:), refpars(:)
REAL*8, ALLOCATABLE :: xexp(:), yexp(:), y2exp(:)
INTEGER nexp, ndw
parameter(nexp=20,ndw=10)

CONTAINS
   SUBROUTINE init_chi2gauss(nmom, exact, ndata, prior, domega)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nmom, ndata
      REAL*8, INTENT(IN) :: exact(nmom), prior(ndata), domega
      INTEGER i
      
      gaussnmom=nmom
      ALLOCATE(gaussref(nmom))
      gaussref=exact
      gaussndata = ndata
      ALLOCATE(xref(ndata), yref(ndata), y2ref(ndata))
      yref = prior
      gaussdomega = domega
      DO i=1,ndata
         xref(i)=(i-1)*domega         
      ENDDO
      CALL spline(xref, yref, ndata, y2ref)
      
      allocate(xexp(nexp*ndw), yexp(nexp*ndw), y2exp(nexp*ndw))
      DO i=1,nexp*ndw
         xexp(i)=dble(i-1.0)/nexp
         yexp(i)=exp(-0.5*xexp(i)**2)         
      ENDDO
      CALL spline(xexp, yexp, nexp*ndw, y2exp)
      ALLOCATE(curspec(ndata))
      ALLOCATE(gauss(ndata), dgauss(ndata), fgauss(ndata))
      icalls = 0
      irun = 0
      chi2alpha = 0.0
   END SUBROUTINE

   REAL*8 FUNCTION lorentz(x,pars)
      REAL*8, INTENT(IN) :: pars(3)
      REAL*8, INTENT(IN) :: x
      lorentz = ((x-pars(2))/pars(3))**2
      lorentz = pars(1)**2/(1.0d0+lorentz)/(pars(3)*3.1415927)
   END FUNCTION
   
   REAL*8 FUNCTION chi2_gauss(pars, npars)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: npars
      REAL*8, INTENT(IN) :: pars(npars)
      
      REAL*8 :: gaussedprior(gaussndata), gaussedmoms(gaussnmom)
      REAL*8 specdist, specw2
      CHARACTER*10 srun
      
      INTEGER i, j
      
      ! computes the trial spectrum
      gaussedprior = 0.0d0
      DO i=1,gaussndata
         do j=1,npars,3
            gaussedprior(i) = gaussedprior(i) + lorentz(xref(i),pars(j:)) 
         enddo
      ENDDO
      
      ! gets the new moments
      CALL get_moments(gaussnmom, gaussndata, gaussedprior, gaussdomega, gaussedmoms)

      chi2_gauss=0.0d0
      DO i=1,gaussnmom
         IF (covariance(1,1).ne. 0.d0) THEN
         ! weight only by the diagonal for now -- the variance
            chi2_gauss = chi2_gauss + (gaussref(i) - gaussedmoms(i))**2/covariance(i,i)
         ELSE
            chi2_gauss = chi2_gauss + (gaussref(i) - gaussedmoms(i))**2
         ENDIF
      ENDDO
      
      specdist = 0.0d0
      specw2 = 0.0d0
      DO i=2,gaussndata-1  ! scalar product between the spectra
         specdist = specdist + yref(i)*gaussedprior(i)
         specw2 = specw2 + yref(i)*yref(i)
      ENDDO 
      specdist = specdist + 0.5*(yref(1)*gaussedprior(1)+yref(gaussndata)*gaussedprior(gaussndata))
      specw2 = specw2 + 0.5*(yref(1)*yref(1)+yref(gaussndata)*yref(gaussndata))
      !write(*,*) specdist, specw2, yref(100), gaussedprior(100)
      
      specdist = -dlog(specdist/specw2)      
      
      
      chi2reg = 0.0d0
      do i=1, npars
         chi2reg = chi2reg + (1.0-pars(i)/refpars(i))**2
      enddo 
      
      ! writes out diagnostics
      write(srun,'(I3)') irun
      IF (MOD(ICALLS,1000)==0) THEN
      if (irun ==0) THEN
         OPEN(666, FILE="GAUSS",status="unknown", action="write")
      ELSE
         
         OPEN(666, FILE=trim("GAUSS."//adjustl(srun)),status="unknown", action="write")
      ENDIF
      write(666,*) "# ", gaussndata, gaussdomega
      DO i=1,gaussndata
         write(666,*) (i-1)*gaussdomega/4.5563353e-06, gaussedprior(i), gauss(i)/4.5563353e-06, dgauss(i), fgauss(i)
      ENDDO
      CLOSE(666)
      if (irun ==0) THEN
         OPEN(666, FILE="GAUSSMOM",status="unknown", action="write")
      ELSE
         OPEN(666, FILE=trim("GAUSSMOM."//adjustl(srun)),status="unknown", action="write")
      ENDIF
      write(666,*) "# ", gaussnmom, chi2_gauss
      DO i=1,gaussnmom
         write(666,*) i, gaussref(i), gaussedmoms(i)
      ENDDO
      CLOSE(666)
      WRITE(*,*) pars
      write(*,*) "Area is ", sum(gaussref)-0.5*(gaussref(1)+gaussref(gaussnmom)), &
        sum(gaussedmoms)-0.5*(gaussedmoms(1)+gaussedmoms(gaussnmom))
        
      write(*,*) "Chi2 is ", chi2_gauss
      write(*,*) "Reg is ", chi2reg
      
      write(*,*) "K-L Divergence distance ", specdist   
      write(*,*) "****************************************"
      !write(*,*) "parameters are", pars
      END IF
      chi2val = chi2_gauss
      chi2_gauss = chi2_gauss + chi2reg*chi2alpha
      
      ICALLS = ICALLS + 1
      !if (ICALLS>10) call exit(0)
   END FUNCTION
END MODULE
