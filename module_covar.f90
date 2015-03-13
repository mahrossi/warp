      module covar
      implicit none
      REAL*8, allocatable :: covariance(:,:), invcov(:,:)
      contains

      subroutine read_covar(unit, cov, nmom)
!
!     Reads the covariance matrix -- each line is one column
!
           implicit none
           integer, intent(in) :: nmom, unit
           logical :: eof
           integer :: i_code, imom, imeasure, imom2, i
           double precision :: trash, averages(nmom)
           double precision :: cov(nmom,nmom)

           do imom=1, nmom
              READ(unit,*,iostat = i_code) (cov(imom, i), i=1, nmom)
           enddo

      end subroutine 

      subroutine diagcov(cov, eigen, nmom)
!
!     Inverts the covariance matrix using lapack
!
          implicit none
          integer, intent(in) :: nmom
          DOUBLE PRECISION, intent(in) :: cov(nmom,nmom)
          CHARACTER          JOBZ, UPLO
          INTEGER            INFO, LDA, LWORK, N
          DOUBLE PRECISION   DUMMY( nmom, nmom ), WORK(3*nmom), eigen(nmom)


          DUMMY=cov
          call DSYEV('V', 'U', nmom, DUMMY, nmom, eigen, WORK, 3*nmom, INFO )
          if (INFO .ne.0) stop 'smth wrong with diag'

          write(*,*) 'eigenvalues', eigen

      end subroutine

      subroutine invertcov(cov, nmom, inv)
!
!     Inverts the covariance matrix using lapack
!
          implicit none
          integer, intent(in) :: nmom
          DOUBLE PRECISION, intent(in) :: cov(nmom,nmom)
          INTEGER            INFO,  M, N
          INTEGER            IPIV( nmom )
          DOUBLE PRECISION   DUMMY( nmom, nmom ), WORK(nmom)
          DOUBLE PRECISION :: inv(nmom, nmom)

          DUMMY=cov
          call DGETRF( nmom, nmom, DUMMY, nmom, IPIV, INFO )
          if (INFO .ne.0) stop 'argh!'
          call DGETRI( nmom, DUMMY, nmom, IPIV, WORK, nmom, INFO )
          if (INFO .ne.0) stop 'argh 2!'
          inv=DUMMY
      end subroutine

      subroutine writecov(cov, inv, nmom)
!
!     Writes the covariance matrix and its inverse, if asked
!
          implicit none
          DOUBLE PRECISION, intent(in) :: cov(nmom,nmom), inv(nmom, nmom)          
          integer :: i, j, nmom
!
!          if (.not.allocated(covariance)) stop 'must calculate covar first'
!          if (.not.allocated(invcov)) stop 'must invert covar first'
!
          write(*,*) 'COVAR'
          do i=1,nmom
             write(*,*) ( cov(i,j), j=1, nmom )
          enddo
          write(*,*) 'MAX ', maxval(abs(cov))
          write(*,*) 'MIN ', minval(abs(cov))
          write(*,*) 'COVAR INVERSE'
          do i=1,nmom
             write(*,*) ( invcov(i,j), j=1, nmom )
          enddo
          write(*,*) 'MAX ', maxval(abs(inv))
          write(*,*) 'MIN ', minval(abs(inv))
      end subroutine

      end module
