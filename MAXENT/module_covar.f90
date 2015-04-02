      module covar
      implicit none
      double precision, allocatable :: covariance(:,:), invcov(:,:), D(:)
      integer*8 :: nmeasure
      contains

      subroutine read_covar(unitcov, unitexact, nmom)
!
!     Reads the covariance matrix -- each line is one column
!
           implicit none
           integer, intent(in) :: nmom, unitcov, unitexact
           logical :: eof
           integer :: i_code, imom, imeasure, imom2, i
           double precision :: trash, averages(nmom), dummy

           allocate(covariance(nmom,nmom), D(nmom))
           write(*,*) covariance

           do imom=1, nmom
              READ(unitcov,*,iostat = i_code) (covariance(imom, i), i=1, nmom)
           enddo

           do imom=1, nmom
              READ(unitexact,*,iostat = i_code) dummy, D(imom)
           enddo

      end subroutine 

      subroutine calculate_covar(nmom)
!
!     Calculates the covariance Matrix - also spills out the average of the imagacf
!     Only non natural name: D is the average of all measurements, here the imagacf accounted for
!     File has to be called 'onlyvvdata_imag.dat' for now, can change later
!
           implicit none
           integer, intent(in) :: nmom
           logical :: eof
           integer :: i_code, imom, imeasure, imom2, i
           double precision, allocatable :: alldata(:)
           double precision :: trash, averages(nmom)

           allocate(covariance(nmom,nmom))

           eof=.false.
           open(10, file='onlyvvdata_imag.dat')
           read (10,*,iostat = i_code) nmeasure
           if (i_code.ne.0) then
             eof = .true.
           end if
           

           allocate(alldata(nmom*nmeasure))
      
           i=1
           do while (.not.eof)
              read (10,*,iostat = i_code) alldata(i)
              if (i_code.ne.0) then
                 eof = .true.
              end if
              i = i+1
           end do
           close(10)

!          calculate averages
           averages(:)=0.d0
           do imom=1, nmom
              do imeasure=1,nmeasure 
                 averages(imom)=averages(imom)+alldata(imeasure+(imom-1)*nmeasure)/nmeasure
              enddo ! imeasure
           enddo !imom

!          calculate full covariance matrix in a very non optimized way - must be able to
!          do through calculating means and outer products on the fly...
           covariance(:,:)=0.d0
           do imom=1, nmom
              do imom2=imom, nmom
                 do imeasure=1, nmeasure
                    covariance(imom, imom2)=covariance(imom, imom2)+(averages(imom)-alldata(imeasure+(imom-1)*nmeasure)) &
                                            *(averages(imom2)-alldata(imeasure+(imom2-1)*nmeasure))
                 enddo !imeasure
              enddo !imom2
           enddo !imom
           ! it is symmetric, so now just copy to the other part
           do imom=1, nmom
              do imom2=imom, nmom
                 covariance(imom2, imom)=covariance(imom, imom2)
              enddo
           enddo

!          and normalize by nmeasure*(nmeasure-1)
           covariance=covariance/(nmeasure*(nmeasure-1))
!          set the D array in Bryan's notation:
           D=averages

      end subroutine 

      subroutine diagcov(nmom)
!
!     Inverts the covariance matrix using lapack
!
          implicit none
          integer, intent(in) :: nmom
          CHARACTER          JOBZ, UPLO
          INTEGER            INFO, LDA, LWORK, N
          DOUBLE PRECISION   DUMMY( nmom, nmom ), WORK(3*nmom), eigen(nmom)

          DUMMY=covariance
          call DSYEV('V', 'U', nmom, DUMMY, nmom, eigen, WORK, 3*nmom, INFO )
          if (INFO .ne.0) stop 'smth wrong with diag'

          write(*,*) 'eigenvalues', eigen

      end subroutine

      subroutine invertcov(nmom)
!
!     Inverts the covariance matrix using lapack
!
          implicit none
          integer, intent(in) :: nmom
          INTEGER            INFO,  M, N
          INTEGER            IPIV( nmom )
          DOUBLE PRECISION   DUMMY( nmom, nmom ), WORK(nmom)
          allocate(invcov(nmom,nmom))
          DUMMY=covariance
          call DGETRF( nmom, nmom, DUMMY, nmom, IPIV, INFO )
          if (INFO .ne.0) stop 'argh!'
          call DGETRI( nmom, DUMMY, nmom, IPIV, WORK, nmom, INFO )
          if (INFO .ne.0) stop 'argh 2!'
          invcov=DUMMY
      end subroutine

      subroutine writecov(nmom)
!
!     Writes the covariance matrix and its inverse, if asked
!
          implicit none
          integer :: i, j, nmom

          if (.not.allocated(covariance)) stop 'must calculate covar first'
          if (.not.allocated(invcov)) stop 'must invert covar first'

          write(*,*) 'COVAR'
          do i=1,nmom
             write(*,*) ( covariance(i,j), j=1, nmom )
          enddo
          write(*,*) 'MAX ', maxval(abs(covariance))
          write(*,*) 'MIN ', minval(abs(covariance))
          write(*,*) 'COVAR INVERSE'
          do i=1,nmom
             write(*,*) ( invcov(i,j), j=1, nmom )
          enddo
          write(*,*) 'MAX ', maxval(abs(invcov))
          write(*,*) 'MIN ', minval(abs(invcov))
      end subroutine

      end module
