MODULE CONSTRAINTS
   IMPLICIT NONE
   
   REAL*8 :: pi
   REAL*8, ALLOCATABLE :: fnlambdas(:,:)   
   PARAMETER (pi=4.0d0*atan(1.0d0))
   CONTAINS
      
   SUBROUTINE init_fn(nmom, ndata, dlambda, domega, betah)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndata, nmom
      REAL*8, INTENT(IN) :: domega, dlambda, betah
      INTEGER i,j
      
      IF (ALLOCATED(fnlambdas)) DEALLOCATE(fnlambdas)
      ALLOCATE(fnlambdas(nmom,ndata))
      DO i=1,nmom
         fnlambdas(i,1) = 1.d0/(pi)
         DO j=2,ndata
            fnlambdas(i,j) =betah*(j-1)*domega*cosh(betah*((j-1)* &
               domega/2.-dlambda*domega*(i-1)*(j-1)))/(2.*pi*sinh(betah*(j-1)*domega/2.))
         ENDDO
      ENDDO
   END SUBROUTINE
   
   SUBROUTINE get_moments(nmom, ndata, fomega, domega, mlambda)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndata, nmom
      REAL*8, INTENT(IN) :: domega, fomega(ndata)
      REAL*8, INTENT(OUT) :: mlambda(nmom)
      
      INTEGER i,j
      DO i=1,nmom
         mlambda(i) = (fnlambdas(i,1) * fomega(1)+fnlambdas(i,ndata) * fomega(ndata))*0.5d0
      ENDDO      
      DO j=2,ndata-1
         DO i=1,nmom
            mlambda(i) = mlambda(i) + fnlambdas(i,j) * fomega(j)
         ENDDO
      ENDDO
      mlambda = mlambda * domega      
   END SUBROUTINE
END MODULE
