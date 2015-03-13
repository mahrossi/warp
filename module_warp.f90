MODULE CHI2WARP
USE CONSTRAINTS
USE SPLINES
USE COVAR
IMPLICIT NONE

INTEGER warpnmom, warpndata, icalls, irun
REAL*8 :: warpdomega, chi2alpha, chi2val, chi2reg
REAL*8, ALLOCATABLE :: warpref(:), xref(:), yref(:), y2ref(:), curspec(:)
REAL*8, ALLOCATABLE :: warp(:), dwarp(:), fwarp(:), refpars(:)
REAL*8, ALLOCATABLE :: xexp(:), yexp(:), y2exp(:)
INTEGER nexp, ndw
parameter(nexp=20,ndw=10)

CONTAINS
   SUBROUTINE init_chi2warp(nmom, exact, ndata, prior, domega)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nmom, ndata
      REAL*8, INTENT(IN) :: exact(nmom), prior(ndata), domega
      INTEGER i
      
      warpnmom=nmom
      ALLOCATE(warpref(nmom))
      warpref=exact
      warpndata = ndata
      ALLOCATE(xref(ndata), yref(ndata), y2ref(ndata))
      yref = prior
      warpdomega = domega
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
      ALLOCATE(warp(ndata), dwarp(ndata), fwarp(ndata))
      icalls = 0
      irun = 0
      chi2alpha = 0.0
   END SUBROUTINE

   SUBROUTINE compute_warp(params, ndim)
   ! This warp is \int dx exp(SUM(gaussian(ck, sigmak, x))),
   ! with means fixed to certain points in the grid.
   ! 1) It assumes the parameters are ordered by c1, sigma1, center1, c2, sigma2,
   ! center2, ... etc.
   ! 2) The integration assumes UNIFORM grid spacing
   ! OUTPUTS: dwarp, warp
   IMPLICIT NONE
   INTEGER, INTENT(IN):: ndim
   REAL*8, INTENT(IN) :: params(ndim)
   integer :: i
   real*8 :: lndwx(ndim),lndwy(ndim), lndwy2(ndim), wi, maxw, dw, rx

   maxw = warpdomega * warpndata
   dw = maxw/(ndim-1)
   do i=1,ndim
      lndwx(i)=(i-1)*dw
      lndwy(i)=params(i)
   enddo

   call spline(lndwx,lndwy,ndim,lndwy2)

   do i=1,warpndata
     call splint(lndwx,lndwy,lndwy2,ndim,xref(i),rx)
     dwarp(i)=dexp(rx)
   enddo
   
   ! Integrate through trapezium rule
   warp(1)=0.0d0   
   do i=2, warpndata
      warp(i) = warp(i-1) + 0.5d0*(dwarp(i-1)+dwarp(i))*warpdomega
   enddo
   END SUBROUTINE

   SUBROUTINE compute_fwarp(params, ndim)
   ! Warp weight as a sum of Gaussians on a grid
   IMPLICIT NONE
   INTEGER, INTENT(IN):: ndim
   REAL*8, INTENT(IN) :: params(ndim)
   integer :: i
   real*8 :: lndwx(ndim),lndwy(ndim), lndwy2(ndim), wi, maxw, dw, rx

   maxw = warpdomega * warpndata
   dw = maxw/(ndim-1)
   do i=1,ndim
      lndwx(i)=(i-1)*dw
      lndwy(i)=params(i)
   enddo

   call spline(lndwx,lndwy,ndim,lndwy2)

   do i=1,warpndata
     call splint(lndwx,lndwy,lndwy2,ndim,xref(i),rx)
     fwarp(i)=dexp(rx)
   enddo
   END SUBROUTINE
   
   REAL*8 FUNCTION chi2_warp(pars, npars)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: npars
      REAL*8, INTENT(IN) :: pars(npars)
      
      REAL*8 :: warpedprior(warpndata), warpedmoms(warpnmom)
      REAL*8 :: vec_diff(warpnmom), dummy(warpnmom)
      REAL*8 specdist, specw2
      CHARACTER*10 srun
      
      INTEGER i
      
      ! computes the warp
      !CALL compute_warp(pars(1:npars), npars)
      CALL compute_warp(pars(1:npars/2), npars/2)
      CALL compute_fwarp(pars(npars/2+1:npars), npars/2)
      
      
      ! warps the prior
      DO i=1,warpndata
         IF (warp(i)<=xref(warpndata)) THEN
            CALL splint(xref, yref, y2ref, warpndata, warp(i), warpedprior(i))
            
            warpedprior(i) = warpedprior(i) * dwarp(i) * fwarp(i)
         ELSE
            warpedprior(i) = 0.0d0
         ENDIF
      ENDDO
      
      ! gets the new moments
      CALL get_moments(warpnmom, warpndata, warpedprior, warpdomega, warpedmoms)

      chi2_warp=0.0d0
      IF (covariance(1,1).ne. 0.d0) THEN
         vec_diff=warpref - warpedmoms
         dummy=MATMUL(invcov, vec_diff)
         chi2_warp=dot_product(vec_diff, dummy)
      ELSE
         DO i=1,warpnmom
           chi2_warp = chi2_warp + (warpref(i) - warpedmoms(i))**2
         ENDDO
      ENDIF
      
      specdist = 0.0d0
      specw2 = 0.0d0
      DO i=2,warpndata-1  ! scalar product between the spectra
         specdist = specdist + yref(i)*warpedprior(i)
         specw2 = specw2 + yref(i)*yref(i)
      ENDDO 
      specdist = specdist + 0.5*(yref(1)*warpedprior(1)+yref(warpndata)*warpedprior(warpndata))
      specw2 = specw2 + 0.5*(yref(1)*yref(1)+yref(warpndata)*yref(warpndata))
      !write(*,*) specdist, specw2, yref(100), warpedprior(100)
      
      specdist = -dlog(specdist/specw2)      
      
      ! writes out diagnostics
      write(srun,'(I3)') irun
      IF (MOD(ICALLS,1000)==0) THEN
      if (irun ==0) THEN
         OPEN(666, FILE="WARP",status="unknown", action="write")
      ELSE
         
         OPEN(666, FILE=trim("WARP."//adjustl(srun)),status="unknown", action="write")
      ENDIF
      write(666,*) "# ", warpndata, warpdomega
      DO i=1,warpndata
         write(666,*) (i-1)*warpdomega/4.5563353e-06, warpedprior(i), warp(i)/4.5563353e-06, dwarp(i), fwarp(i)
      ENDDO
      CLOSE(666)
      if (irun ==0) THEN
         OPEN(666, FILE="WARPMOM",status="unknown", action="write")
      ELSE
         OPEN(666, FILE=trim("WARPMOM."//adjustl(srun)),status="unknown", action="write")
      ENDIF
      write(666,*) "# ", warpnmom, chi2_warp
      DO i=1,warpnmom
         write(666,*) i, warpref(i), warpedmoms(i)
      ENDDO
      CLOSE(666)
      write(*,*) "Area is ", sum(warpref)-0.5*(warpref(1)+warpref(warpnmom)), &
        sum(warpedmoms)-0.5*(warpedmoms(1)+warpedmoms(warpnmom))
      write(*,*) "Reg 1,2,3 of warp are ", reg1(warp,warpndata), reg2(dwarp,warpndata), reg3(dwarp,warpndata)
      write(*,*) "Reg 2,3 of fwarp are ", reg2(fwarp,warpndata), reg3(fwarp,warpndata)
      write(*,*) "EMD Warp distance is", regemw(warp,dwarp,yref,fwarp,warpndata)
      write(*,*) "EMD distance is", emd(warpedprior,yref,warpndata)
      write(*,*) "Chi2 is ", chi2_warp
      write(*,*) "K-L Divergence distance ", specdist   
      write(*,*) "****************************************"
      !write(*,*) "parameters are", pars
      END IF
      chi2val = chi2_warp
      chi2reg = emd(warpedprior,yref,warpndata)
      !chi2reg = regemw(warp,dwarp,yref,fwarp,warpndata)*1d-5 !(reg1(warp,warpndata) * 1d-6 + reg2(dwarp,warpndata) * 1d-6  + reg3(dwarp, warpndata) * 0d-17 & 
                !            + reg2(fwarp, warpndata) *1d-6 + reg3(fwarp, warpndata) * 0d-17)
      chi2_warp = (chi2_warp + chi2reg*chi2alpha)*1e10
      
      ICALLS = ICALLS + 1
      !if (ICALLS>10) call exit(0)
   END FUNCTION
   
   REAL*8 FUNCTION emd(f1, f2, ndim)
   ! Implements earth mover distance in 1D
      REAL*8 f1(ndim), f2(ndim), cf1, cf2
      INTEGER i, ndim
      
      cf1=0.5*f1(1)
      cf2=0.5*f2(1)
      emd=abs(cf1-cf2)*0.5
      do i=2, ndim-1 
         cf1 = cf1 + f1(i)
         cf2 = cf2 + f2(i)
         emd = emd + abs(cf1-cf2)
      enddo
      cf1 = cf1 + f1(ndim)*0.5
      cf2 = cf2 + f2(ndim)*0.5
      emd=emd+abs(cf1-cf2)*0.5
      
      emd = emd * warpdomega
   END FUNCTION
   
   
   REAL*8 FUNCTION regemw(fun,dfun,pfun,ffun,ndim)
   ! Implements regularization \int ((W[w]-w) P(w))^2  dw
      REAL*8 fun(ndim), dfun(ndim), pfun(ndim), ffun(ndim)
      INTEGER i, ndim
      
      regemw=0
      do i=2, warpndata-1
      !do i=warpndata/2, warpndata-1
         regemw=regemw+((dabs(fun(i)-xref(i))*pfun(i))**2+(dabs(1-dfun(i))**2)+(dabs(1-ffun(i))**2) )
      enddo
      regemw = regemw + 0.5*(( (dabs(fun(1)-xref(1))*pfun(1))**2+dabs(1-dfun(1))**2+(dabs(1-ffun(1))**2) ) &
      +((dabs(fun(warpndata)-xref(warpndata))*pfun(warpndata))**2 + &
      dabs(1-dfun(warpndata))**2 +dabs(1-ffun(warpndata))**2 ) )
      !reg1 = reg1 + 0.5*((1-fun(warpndata)/xref(warpndata))**2) ! fun(1) is probably zero
      regemw = regemw*warpdomega
   END FUNCTION
   
   
   REAL*8 FUNCTION reg1(fun,ndim)
   ! Implements regularization \int (1-W[w]/w)^2 dw
      REAL*8 fun(ndim)
      INTEGER i, ndim
      
      reg1=0
      do i=2, warpndata-1
      !do i=warpndata/2, warpndata-1
         reg1=reg1+(fun(i)-xref(i))**2.d0
      enddo
      reg1 = reg1 + 0.5*((fun(1)-xref(1))**2+(fun(warpndata)-xref(warpndata))**2)
      !reg1 = reg1 + 0.5*((1-fun(warpndata)/xref(warpndata))**2) ! fun(1) is probably zero
      reg1 = reg1*warpdomega
   END FUNCTION
   
   REAL*8 FUNCTION reg2(fun, ndim)
   ! Implements regularization \int (W'[w]-1)^2 dw
      REAL*8 fun(ndim)
      INTEGER i, ndim
      
      reg2=0
      do i=2, warpndata-1
         reg2=reg2+(fun(i)-1)**2.d0
      enddo
      reg2 = reg2 + 0.5*((fun(1)-1)**2+(fun(warpndata)-1)**2)
      reg2 = reg2 * warpdomega
   END FUNCTION

   REAL*8 FUNCTION reg3(fun, ndim)
      
   ! Implements regularization \int (W''[w])**2 dw
      REAL*8 fun(ndim)
      INTEGER i, ndim
      REAL*8 d2
      
      reg3=0
      do i=2, warpndata-1
         d2 = ((fun(i+1)-fun(i-1))/(2.0*warpdomega))**2.d0
         reg3=reg3+ d2 * xref(i)
      enddo
      reg3 = reg3 * warpdomega
   END FUNCTION
END MODULE
