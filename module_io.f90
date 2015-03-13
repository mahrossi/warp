MODULE MCIO

IMPLICIT NONE


CONTAINS
   SUBROUTINE readgrid(nunit, ydata, ndata, dx)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nunit
      INTEGER, INTENT(out) :: ndata
      REAL*8, INTENT(out) :: dx
      REAL*8, INTENT(out), ALLOCATABLE :: ydata(:)
      INTEGER :: i
      CHARACTER*20 :: dummy
      
      READ(nunit,*) dummy, ndata, dx
      IF (ALLOCATED(ydata)) DEALLOCATE(ydata)
      ALLOCATE(ydata(ndata))
      
      DO i=1,ndata
         READ(nunit,*) dummy, ydata(i)
      ENDDO   
   END SUBROUTINE

   SUBROUTINE writegrid(nunit, ydata, ndata, dx)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nunit
      INTEGER, INTENT(IN) :: ndata
      REAL*8, INTENT(IN) :: dx
      REAL*8, INTENT(IN) :: ydata(ndata)
      INTEGER :: i
      
      WRITE(nunit,*) "# ", ndata, dx
      DO i=1,ndata
         WRITE(nunit,*) (i-1)*dx, ydata(i)
      ENDDO   
   END SUBROUTINE

END MODULE
