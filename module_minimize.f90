MODULE MINIMIZE
IMPLICIT NONE
CONTAINS

SUBROUTINE mebs_dhsimplex(fun, ndim, vertices, fvals, ftol, iter)
IMPLICIT NONE
INTEGER :: iter,ndim
real*8, external :: fun
REAL*8 :: ftol,vertices(ndim+1,ndim),fvals(ndim+1),funk,TINY, dummy
PARAMETER (TINY=1.e-11) 
! EXTERNAL funk
! MR: I will transform funk into a call to an external routine get_fn, with
! the function chosen by 'choice', which is a global variable in the fn module
! and is set in the input
! Multidimensional minimization of the function funk(x) where x(1:ndim) is a 
! vector in ndim dimensions, by the downhill simplex method of Nelder and Mead. 
! The matrix p->vertices(1:ndim+1,1:ndim) is input. 
! Its ndim+1 rows are ndim-dimensional vectors which are the vertices 
! of the starting simplex. Also input is the vector y(1:ndim+1)->fvals, 
! whose components must be pre-initialized to the values of 
! funk evaluated at the ndim+1 vertices (rows) of vertices; 
! and ftol the fractional convergence tolerance to be achieved in the 
! function value (n.b.!). On output, p and y will have been reset to ndim+1 
! new points all within ftol of a minimum function value, and iter gives the 
! number of function evaluations taken.
INTEGER i,ihi,ilo,inhi,j,m,n
REAL*8 rtol,sum,swap,ysave,ytry,psum(ndim)

iter=0
1 do n=1,ndim 
   sum=0.
      do m=1,ndim+1
         sum=sum+vertices(m,n)
      enddo
   psum(n)=sum 
enddo

2 ilo=1
if (fvals(1).gt.fvals(2)) then
   ihi=1
   inhi=2 
else
   ihi=2
   inhi=1
endif

do i=1,ndim+1 
   if(fvals(i).le.fvals(ilo)) ilo=i 
   if(fvals(i).gt.fvals(ihi)) then
      inhi=ihi
      ihi=i
   else if(fvals(i).gt.fvals(inhi)) then
      if(i.ne.ihi) inhi=i
   endif 
enddo

rtol=2.*abs(fvals(ihi)-fvals(ilo))/(abs(fvals(ihi))+abs(fvals(ilo))+TINY)

if (rtol.lt.ftol) then 
   swap=fvals(1)
   fvals(1)=fvals(ilo)
   fvals(ilo)=swap 
   do n=1,ndim
      swap=vertices(1,n)
      vertices(1,n)=vertices(ilo,n)
      vertices(ilo,n)=swap
   enddo 
!   write(*,*) 'MADE IT!!'
! Here the code found the minimum
   dummy=fun(vertices(1,:), ndim) !, ndim, dummy, .TRUE.)
   return
endif

!write(*,*) 'I am in iteration ', iter
! if (iter.ge.ITMAX) stop 'ITMAX exceeded in mebs_dhsimplex'

iter=iter+2

call mebstry(fun,vertices,fvals,psum,ndim,ihi,-1.d0, ytry) 

if (ytry.le.fvals(ilo)) then
   call mebstry(fun,vertices,fvals,psum,ndim,ihi,2.d0, ytry)
else if (ytry.ge.fvals(inhi)) then
   ysave=fvals(ihi)
   call mebstry(fun,vertices,fvals,psum,ndim,ihi,0.5d0, ytry)
   if (ytry.ge.ysave) then
      do i=1,ndim+1
         if(i.ne.ilo)then
            do j=1,ndim 
               psum(j)=0.5*(vertices(i,j)+vertices(ilo,j))
               vertices(i,j)=psum(j)
            enddo 
            funk=fun(psum, ndim)!, ndim, funk, .FALSE.)
!            write(*,*) 'function value', funk
            fvals(i)=funk
         endif
      enddo 
      iter=iter+ndim 
      goto 1
   endif 
else
   iter=iter-1 
endif
goto 2 

END SUBROUTINE mebs_dhsimplex

SUBROUTINE mebstry(fun,vertices,y,psum,ndim,ihi,fac,output) 
IMPLICIT NONE
REAL*8 y(ndim+1)
INTEGER ihi,ndim
real*8, external :: fun
REAL*8 output,fac,vertices(ndim+1,ndim),psum(ndim)
! Extrapolates by a factor fac through the face of the simplex across 
! from the high point, tries it, and replaces the high point if the new 
! point is better.
INTEGER j
REAL*8 fac1,fac2,ytry,ptry(ndim), funk 

fac1=(1.-fac)/ndim 
fac2=fac1-fac

do j=1,ndim
   ptry(j)=psum(j)*fac1-vertices(ihi,j)*fac2 
enddo

funk=fun(ptry,ndim) !, ndim, funk, .FALSE.)

ytry=funk

if (ytry.lt.y(ihi)) then
   y(ihi)=ytry 
   do j=1,ndim
      psum(j)=psum(j)-vertices(ihi,j)+ptry(j)
      vertices(ihi,j)=ptry(j) 
   enddo 
endif 
output=ytry 
return

END SUBROUTINE mebstry

END MODULE MINIMIZE
