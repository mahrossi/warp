program main
   use covar
   implicit none
   integer :: npoints, nmom
   double precision, allocatable :: p(:), pmodel(:), weight(:), temp(:) 
   double precision, allocatable :: momentfunc(:,:), w(:), lambda(:), Ttransp(:,:), U(:,:), S(:), Vtransp(:,:)
   double precision :: dw, betah, temperature, hbar, value, negentropy, alpha, chisquare, falpha
   integer :: ipoint, imom, i_code, i, j, nbeads, info, nalpha, i_alpha
   logical :: eof, exist
   character*40 :: trash, covfile, exactfile, priorfile

!!!!! system dependent parameters !!!!!
   temperature=25.0 ! Kelvin
   nbeads=24
   betah=7.63823305/temperature ! beta*hbar in these units, for temperature in Kelvin
   hbar=6.3507799345342777 ! in au*A^2/ps
   nalpha=20
   falpha=0.5
   alpha=2000.
   covfile='covar.dat' 
   exactfile='exact.dat'
   priorfile='prior.dat' 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   nmom=(nbeads/2)+1

   write(*,*) 'nmom', nmom

   ALLOCATE(lambda(nmom))

   do j=1, nmom
      lambda(j)=(j-1)*betah/nbeads
   enddo


! now lets read and store the prior spectrum
! first get dimensiont
! Get dimensions and read frequencies and model spectrum
   npoints=0
   eof=.false.
   open(10, file=priorfile)

   do while (.not.eof)
     read (10,*,iostat = i_code) trash, trash
     if (i_code.ne.0) then
        eof = .true.
     end if
     npoints = npoints+1
   end do!

   npoints = npoints-1

   write(*,*) npoints

   close(10)

   ALLOCATE(p(npoints), pmodel(npoints), weight(npoints), momentfunc(nmom,npoints), w(npoints), temp(nmom), &
            Ttransp(npoints, nmom), U(npoints, nmom), S(nmom), Vtransp(nmom, nmom))

   call readprior(pmodel, w, npoints, priorfile)
   dw=w(2)-w(1) ! ASSUMES EQUAL SPACING!!
   p(:)=0.d0

! Prior and frequencies are read, quantities are allocated. 
! Now get momentfunc matrix, which will suffer an SVD very soon
   do ipoint=1, npoints
      do imom=1, nmom
         call fvv(w(ipoint), betah, lambda(imom), value)
         momentfunc(imom,ipoint)=value
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first let's read and calculate the covariance matrix,
! the average of the imagacf and invert the covariance matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   call calculate_covar(nmom)
!   call invertcov(nmom)
!   call diagcov(nmom)
   OPEN(11, FILE=exactfile)
   inquire(file=covfile, exist=exist)
   if (exist) then
     OPEN(12, FILE=covfile)
     CALL read_covar(12, 11, nmom) 
   ! Already calculates the inverse covariance
     CALL invertcov(nmom)
     CLOSE(12)
   end if
   CLOSE(11)


   open(unit=20, file='imaginary.out', status='unknown')
        do j=1, nmom
           write(20, *) lambda(j)/betah, D(j), sqrt(covariance(j,j)/100)
        enddo
   close(unit=20, status='keep')

!
! test to see if imaginary data is fitted
!
      temp=matmul(momentfunc, dw*pmodel)
      open(unit=20, file='initial_imaginary.out', status='unknown')
           do j=1, nmom
              write(20, *) lambda(j)/betah, temp(j)
           enddo
      close(unit=20, status='keep')
!
! end test
!

! I will do SVD out of the maxent routine:
!
! Transpose momentfunc
   Ttransp=transpose(momentfunc)  

   call svd(Ttransp, U, S, Vtransp,nmom,npoints) 
!
! End of SVD
!
!
!  Start alpha loop and maxent calls
!
   open(unit=9, file='alpha_chisquare_entropy.out', status='unknown')
   do i_alpha=1, nalpha, 1
      negentropy=0.d0
      call maxent(p, pmodel, dw, npoints, momentfunc, nmom, alpha, chisquare, Ttransp, U, S, Vtransp)
!
! test to see if imaginary data is fitted
!
      temp=matmul(momentfunc, dw*p)
      open(unit=20, file='reconstructed_imaginary.out', status='unknown')
           do j=1, nmom
              write(20, *) lambda(j)/betah, temp(j)
           enddo
      close(unit=20, status='keep')
!
! end test
!
      ! calculate entropy here
      do ipoint=1, npoints
         if (pmodel(ipoint).gt.0.d0.and.p(ipoint).gt.0.d0) then
            negentropy=negentropy-dw*(p(ipoint)-pmodel(ipoint)-p(ipoint)*log(p(ipoint)/pmodel(ipoint)))
         endif
      enddo
      write(9, *) alpha, log10(negentropy), log10(chisquare)
      alpha=alpha*falpha
   enddo ! i_alpha
   close (unit=9, status='keep')

!
! End alpha loop and maxent calls
!

   open(unit=12, file='result.out', status='unknown')
   do ipoint=1, npoints
      write(12, *) w(ipoint), p(ipoint)
   enddo !ipoint
   close (unit=12, status='keep')
   return

end program main

subroutine maxent(p, pmodel, dw, npoints, momentfunc, nmom, alpha, chisquare, Ttransp, U, S, Vtransp)
! remember that pmodel is still not multiplied by dw
   use covar
   implicit none
   integer, intent(in) :: npoints, nmom
   double precision, intent(in) :: dw, pmodel(npoints), momentfunc(nmom, npoints), alpha
   double precision, intent(in) :: Ttransp(npoints, nmom), U(npoints, nmom), S(nmom), Vtransp(nmom,nmom)                        
   double precision, intent(inout) :: p(npoints), chisquare 
   double precision, allocatable :: F(:)
   double precision, allocatable :: uvec(:), FminusD(:), dummy2(:), gsing(:), sigma(:,:), Us(:,:), dgdu(:,:)
   double precision, allocatable :: g(:), M(:,:), K(:,:), pdiag(:,:), dummy3(:,:), dummy4(:,:), duvec(:),i1(:),i2(:)
   double precision :: mu, v1, v2, numerator, denominator, temp, error,  basis, wp(npoints), negentropy
   integer :: sing, niter
   integer, allocatable :: ipiv(:)
   integer :: ipoint, imom, i_code, i, j, info, i_iter
   logical :: accepted


   niter=20000


! Now let's allocate more stuff

   ALLOCATE(F(nmom), FminusD(nmom), &
            sigma(nmom,nmom), dummy2(nmom), g(nmom), M(nmom, nmom))


   ! initial criterion for singular space dimension
   sigma=0.d0
   sing=nmom
   do i=1, sing
      write(*,*) S(i)
!      if (S(i) .gt. maxval(S)*1d-15) then ! let's start with this criterion for the singular space.
!         sing = sing+1
         sigma(i,i)=S(i)
!      endif
   enddo
!   write(*,*) sigma

   ! Now let's allocate stuff that is in the singular space
  
   ALLOCATE(uvec(sing), pdiag(npoints, npoints), Us(npoints, sing), dummy3(npoints, sing), &
        duvec(sing), dummy4(sing, sing), ipiv(sing), i1(sing), i2(sing), gsing(sing), K(sing, sing), dgdu(sing, sing))

   uvec=0.d0
   duvec=0.d0
   Us(:,:)=U(:, 1:sing)

! START LOOP

   do i_iter=1, niter
   
      mu=1.0
   
      do ipoint=1, npoints
         basis=0.d0
         do i=1, sing
            basis=basis+Us(ipoint, i)*uvec(i)
         enddo
         p(ipoint)=pmodel(ipoint)*exp(basis) 
         wp(ipoint)=dw*p(ipoint)!pmodel(ipoint)*exp(basis)
   !      write(*,*) 'wp', wp(ipoint)
      enddo
   
   ! just to be perfect, let's treat correctly the endpoints
      wp(1)=wp(1)/2.d0
      wp(npoints)=wp(npoints)/2.d0
         
   
   !! ENTROPY ESTIMATE ALSO HERE
   
      negentropy=0.d0
      do ipoint=1, npoints
         if (pmodel(ipoint).gt.0.d0.and.p(ipoint).gt.0.d0) then
             negentropy=negentropy-dw*(p(ipoint)-pmodel(ipoint)-p(ipoint)*log(p(ipoint)/pmodel(ipoint)))
         endif
      enddo
      write(*,*) 's', negentropy
   !!
   
   
      F=matmul(momentfunc, wp)
   !   do i=1, nmom
   !      write(*,*) 'FD: ', F(i), D(i)
   !   enddo
      FminusD=F-D
      
      ! calculate chisquare
      dummy2=matmul(invcov, FminusD)
      chisquare=0.5*dot_product(FminusD, dummy2)
      write(*,*) 'alpha and chi', alpha, chisquare
   
   
      ! calculate g
      dummy2=matmul(FminusD, invcov)
      g=matmul(Vtransp, dummy2)
      g=matmul(sigma, g)
   
      gsing=g(1:sing)
     
      ! calculate dg/du
      M=matmul(transpose(Vtransp), sigma)
      M=matmul(invcov, M)
      M=matmul(Vtransp, M)
      M=matmul(sigma, M)
   !   write(*,*) 'M', M
   
      pdiag=0.d0
      do i=1, npoints
        pdiag(i,i)=wp(i) 
      enddo
      
   !  write(*,*) Us
   
      dummy3=matmul(pdiag, Us)
      K=matmul(transpose(Us), dummy3)
   !   write(*,*) 'K', K
   
   ! HERE WE CAN CALCULATE THE ERROR AS DEFINED BY BRYAN
      if (i_iter .gt. 1) then
         i1=matmul(K, uvec)
         i1=-alpha*i1
         denominator=dot_product(i1, i1)
         denominator=sqrt(denominator)
         i2=matmul(K, gsing)
         temp=dot_product(i2,i2)
         denominator=(denominator+sqrt(temp))**2
         i1=i1-i2
         numerator=2.d0*dot_product(i1, i1)
         error=numerator/denominator
         write(*,*) 'Iter =', i_iter, ' Error =', error
         if (error .lt. 1d-10) then
             write(*,*) 'Converged!'
             ! DEALLOCATE BEFORE RETURNING, JUST IN CASE
             DEALLOCATE( F, FminusD, sigma, dummy2, g, M, &
                 uvec, pdiag, Us, dummy3, duvec, dummy4, ipiv, i1, i2, gsing, K, dgdu)
             return
         endif
      endif
   ! END CALCULATION OF CRAZY ERROR
   
    
      dgdu=matmul(M(1:sing, 1:sing),K)
   
      accepted=.false.
   
      do while (accepted.eqv..false.) 
         dummy4=0.d0
         do i=1, sing
           dummy4(i,i)=(alpha+mu)
         enddo  
   
         dummy4=(dummy4+dgdu)
   
         duvec=(-alpha*uvec-gsing)
   
   
   !   ! now must solve a system that is dummy4*duvec=duvec
         call dgesv(sing, 1, dummy4, sing, ipiv, duvec, sing, info )
         if (info.ne.0) stop 'OPS'
   
   !! CHECK MU
        i1=matmul(K, duvec)
        v1=dot_product(duvec, i1)
        v2=0.d0
        do i=1, npoints
           v2=v2+pmodel(i)
        enddo
        v2=v2
   !! reject step if v1 more than v2
        if (v1.gt.0.2*v2) then
           write(*,*) ' Step rejected. Mu too small: ', mu, v1, v2
           mu=mu*10
        else
           accepted=.true.
        endif
     enddo !accepted
   
     uvec=uvec+duvec
   
  enddo !i_iter

! If it does not converge in niter, stop... 
  write(*,*) 'Could not converge iter, for alpha=', alpha

! DEALLOCATE STUFF
   DEALLOCATE(F, FminusD, sigma, dummy2, g, M, &
              uvec, pdiag, Us, dummy3, duvec, dummy4, ipiv, i1, i2, gsing, K, dgdu)

  return

end subroutine

subroutine readprior(pmodel, w, npoints, priorfile)
   implicit none
   real*8 :: pmodel(npoints), w(npoints)
   integer :: npoints
   integer :: i_code, ipoint, i
   logical :: eof
   character*40 :: trash, priorfile

   pmodel(:)=0.d0
   w(:)=0.d0

   eof = .false.
   open(11, FILE=priorfile)

   ! reads the file now to keep frequency and spectral information
   i=1
   do while (.not.eof)
   read (11,*,iostat = i_code) w(i), pmodel(i)
      if (i_code.ne.0) then
        eof = .true.
      end if
      i = i+1
   end do
   close(11)
! clean pmodel
   do ipoint=1, npoints
      if (pmodel(ipoint).lt.maxval(pmodel)/5000.) then
         pmodel(ipoint)=0.d0
      endif
   enddo !ipoint
! pmodel is clean!

end subroutine

subroutine fvv(omega, betah, lambda, value)
! vv moment function at any beta
   implicit none
   double precision, intent(in) :: omega, betah, lambda
   double precision, intent(out) :: value
   double precision :: hbar, pi
   pi=4.d0*atan(1.d0)
   if (omega.eq.0.d0) then
      value=1.d0/(pi)
   else
      value=betah*omega*cosh(betah*omega/2.-lambda*omega)/(2.*pi*sinh(betah*omega/2.))
   endif
   return
end subroutine

subroutine svd(Ttransp, U, S, Vtransp,nmom, npoints)
!
!     SVD of Ttransp using lapack
!
          implicit none
          INTEGER            INFO,  nmom, npoints, LWORK
          DOUBLE PRECISION   Ttransp(npoints,nmom), U(npoints,nmom), S(nmom), Vtransp(nmom,nmom)
          DOUBLE PRECISION   DUMMY(npoints,nmom)
          DOUBLE PRECISION, allocatable :: WORK(:) 
          DUMMY=Ttransp
          LWORK=6*nmom+2*npoints
          allocate(WORK(LWORK))

          call DGESVD('S','S',npoints, nmom, DUMMY, npoints, S, U, npoints, Vtransp, nmom, WORK, LWORK, INFO  )
          if (INFO .ne.0) stop 'no luck!'

          deallocate(WORK)
end subroutine

