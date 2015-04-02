PROGRAM moments
   USE MCIO
   USE CONSTRAINTS
!   USE CHI2GAUSS
   USE CHI2WARP
   USE MINIMIZE
   USE COVAR
   IMPLICIT NONE

   INTEGER :: nprior, nexact
   REAL*8 :: domega, dlambda, betah
   REAL*8, ALLOCATABLE :: prior(:), exact(:), moms(:)
   CHARACTER*256 :: priorfile, exactfile, covfile, parsfile
   REAL*8 :: omfac
   REAL*8, ALLOCATABLE :: fixpars(:), inpars(:)
   INTEGER i, j, k, niter, nalpha, ngauss, pl
   REAL*8 :: chi2, regalpha, falpha, temp
   REAL*8, ALLOCATABLE :: simplex(:,:), simplexval(:)
   LOGICAL :: exist
   parameter(omfac=4.5563353e-06) !invcm to atomic frequency
   NAMELIST /input/ priorfile, exactfile, covfile, parsfile, temp, ngauss, regalpha, falpha, nalpha
   
   ! DEFAULT INPUT PARAMETERS
   priorfile = 'prior.dat'
   exactfile = 'exact.dat'
   covfile = 'covar.dat' 
   parsfile = ''
   temp = 100
   ngauss = 40
   regalpha=1d-1
   nalpha = 10
   falpha = 0.5
   ! END INPUT PARAMETER
   READ(5,input)
   betah = 1./(temp*3.1668152e-06)
   

   ! Reads prior spectrum and exact imaginary-time moments (and covariance)
   OPEN(11, FILE=priorfile)
   CALL READGRID(11, prior, nprior, domega) 
   CLOSE(11)
   
   domega = domega*omfac 
   
   OPEN(11, FILE=exactfile)
   CALL READGRID(11, exact, nexact, dlambda) 
   CLOSE(11)

   ALLOCATE(covariance(nexact, nexact))
   ALLOCATE(invcov(nexact, nexact))
   covariance=0.d0
   invcov=0.d0
   inquire(file=covfile, exist=exist)
   if (exist) then
     OPEN(11, FILE=covfile)
     CALL read_covar(11, covariance, nexact) 
   ! Already calculates the inverse covariance
     CALL invertcov(covariance, nexact, invcov)
     CLOSE(11)
   end if
!   CHECK IF ALL IS FINE
!   CALL writecov(covariance, invcov, nexact)
   
   CALL init_fn(nexact, nprior, dlambda, domega, betah)   
   ALLOCATE(moms(nexact))
   
   ! get initial moments for the given prior
   CALL get_moments(nexact, nprior, prior, domega, moms)      
   DO i=1,nexact
      WRITE(*,*) (i-1)*dlambda, exact(i), moms(i)
   ENDDO
   
   
   CALL init_chi2warp(nexact, exact, nprior, prior, domega) 
   ! CALL init_chi2gauss(nexact, exact, nprior, prior, domega) 
   chi2alpha = regalpha
   
   
   ! try to compute warp
   ALLOCATE(fixpars(ngauss), inpars(ngauss), refpars(ngauss))
 
   ! initialize fixpars in case they were not initialized
   ! RANDOM INIT:  
   ! do i=1,ngauss
   !    fixpars(i)=(rand()-0.5)*0.01
   !    !fixpars(i)=(rand())*10
   ! enddo
   ! INIT WITH WARP 1 to 1 with the original grid:
    do i=1,ngauss/2
      pl=nprior/ngauss*i
      fixpars(i)=domega*(pl-1)
    enddo
   ! multiplying function keeps the old random init
    do i=(ngauss/2)+1,ngauss
      fixpars(i)=(rand()-0.5)*0.01
    enddo
   ! END INIT -- if there is a file it is read below 
   ! and it overwrites this
   write(*,*) 'fixi', fixpars
   
   inpars = 0.0d0
   refpars = 0.0d0
   IF (parsfile/='') THEN
      OPEN(11, FILE=parsfile)
      CALL READGRID(11, inpars, ngauss, temp) 
      CLOSE(11)
      fixpars=inpars
      refpars = inpars      
   ENDIF
   
   ALLOCATE(simplex(ngauss+1,ngauss), simplexval(ngauss+1))
   
   call srand(12345) 
!  COMPLETELY RANDOM INITIALIZATION
!   do i=1,ngauss+1
!      do j=1,ngauss
!        simplex(i,j)=fixpars(j)*(1+(rand()-0.5)*0.1)
!      enddo
!!      simplexval(i)=chi2_gauss(simplex(i,:),ngauss)      
!      simplexval(i)=chi2_warp(simplex(i,:),ngauss)      
!   enddo
! LESS RANDOM INITIALIZATION
   do i=1,ngauss+1
      simplex(i,:)=fixpars(:)
      if (i.gt.1) then
         simplex(i,i-1)=fixpars(i)+fixpars(i)*(1+(rand()-0.5)*0.1)
      endif
      simplexval(i)=chi2_warp(simplex(i,:),ngauss)
   enddo


   if (nalpha .gt. 0) THEN
      OPEN(1234, FILE="LCURVE",STATUS="unknown")   
      do k=1,nalpha
         chi2alpha = regalpha
         write (*,*) "Starting minimization with alpha = ", chi2alpha
! SIMPLEX WOULD HERE BE RANDOMIZED EVEN FURTHER     
!      do i=1,ngauss+1
!         do j=1,ngauss
!           !simplex(i,j)=simplex(i,j)*(1+(rand()-0.5)*0.1)
!           simplex(i,j)=simplex(i,j)*(1+(rand()-0.5)*0.5) ! randomizes a bit the simplex 
!           !simplex(i,j)=fixpars(i)*(1+(rand()-0.5))+fixpars(j)
!         enddo
!!         simplexval(i)=chi2_gauss(simplex(i,:),ngauss)      
!         simplexval(i)=chi2_warp(simplex(i,:),ngauss)      
!      enddo
! JUST UNCOMMENT ABOVE AND COMMENT LAST BLOCK TO GO BACK
         irun = k
!         CALL mebs_dhsimplex(chi2_gauss, ngauss, simplex, simplexval, 1d-7, niter)
         CALL mebs_dhsimplex(chi2_warp, ngauss, simplex, simplexval, 1d-10, niter)
         write(1234,*) chi2alpha, chi2val, chi2reg
         regalpha = regalpha * falpha
! RESTART AND RERANDOMIZE
         do i=1,ngauss+1
            if (i.gt.1) then
               simplex(i,i-1)=simplex(i, i-1)+simplex(i, i-1)*(1+(rand()-0.5)*0.1)
            endif
            simplexval(i)=chi2_warp(simplex(i,:),ngauss)
         enddo
! COMMENT THIS BLOCK AND UNCOMMENT ABOVE IF NEEDED
      enddo
      CLOSE(1234)
   ELSE
      chi2alpha = regalpha
      irun = 0
!      CALL mebs_dhsimplex(chi2_gauss, ngauss, simplex, simplexval, 1d-7, niter)
      CALL mebs_dhsimplex(chi2_warp, ngauss, simplex, simplexval, 1d-7, niter)
      write(*,*) "Final: ",chi2alpha, chi2val, chi2reg
   ENDIF
   !DO i=1,nprior
   !   WRITE(*,*) (i-1)*domega, warp(i), dwarp(i)
   !ENDDO
END PROGRAM
