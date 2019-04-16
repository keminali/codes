!-- File  verify.f90  ------------------------------------------------------
!
!  Performs several verification calculations given a file of grid spacings 
!  and some observed quantity corresponding to each grid spacing.
!
!  Computes:
!  - order of convergence
!  - Richardson extrapolation to zero grid spacing
!  - grid convergence indices (GCI)
!
!--------------------------------------------------------------------------

   PROGRAM verify

   IMPLICIT NONE

   INTEGER :: n
   INTEGER :: nd

   INTEGER, PARAMETER :: ndmax = 10

   REAL :: ratio
   REAL :: fexact
   REAL :: fsafe
   REAL :: p

   REAL, DIMENSION(ndmax) :: f
   REAL, DIMENSION(ndmax) :: gcif
   REAL, DIMENSION(ndmax) :: r
   REAL, DIMENSION(ndmax) :: x

!--------------------------------------------------------------------------

!..Write Header.

   WRITE(*,*) ' '
   WRITE(*,*) '--- VERIFY: Performs verification calculations ---'

!..Read in the file, determine the number of data points, and output.

   DO n = 1, ndmax
     READ(*,*,end=10) x(n), f(n)
   ENDDO

10 CONTINUE

   nd = n - 1

   WRITE(*,*) ' '
   WRITE(*,*) 'Number of data sets read = ', nd

   WRITE(*,*) ' '
   WRITE(*,*) '     Grid Size     Quantity'
   WRITE(*,*) ' '

   DO n = 1, nd
     WRITE(*,'(2(F14.6))') x(n), f(n)
   ENDDO

!..Compute the grid refinement ratio, r, between each pair.

   DO n = 1, nd-1
     r(n) = x(n+1) / x(n)
   ENDDO

!..Estimate the order of convergence using the first three data pairs 
!..and assuming that the grid refinement ratio is constant, r(1) = r(2).
!..This is done using Eqn. 5.10.6.1 of (Roache, 1998).

   p = log( ( f(3) - f(2) ) / ( f(2) - f(1) ) )  /  log( r(1) )

   WRITE(*,*) ' '
   WRITE(*,*) 'Order of convergence using first three finest grid '
   WRITE(*,*) 'and assuming constant grid refinement (Eqn. 5.10.6.1) '
   WRITE(*,*) 'Order of Convergence, p = ', p

!..Perform Richardson extrapolation to estimate a zero grid value of f.

   fexact = f(1) + ( f(1) - f(2) ) / ( r(1)**p - 1.0 )

   WRITE(*,*) ' '
   WRITE(*,*) 'Richardson Extrapolation: Use above order of convergence'
   WRITE(*,*) 'and first and second finest grids (Eqn. 5.4.1) '
   WRITE(*,*) 'Estimate to zero grid value, f_exact = ', fexact

!..Compute Grid Convergence Index (GCI) for each fine grid using Eqn. 5.6.1
!..from Roache's book. Use factor of safety as recommended on page 123.

   IF ( nd .ge. 3 ) then
     fsafe = 1.25
   ELSE
     fsafe = 3.0
   ENDIF

   DO n = 1, nd-1
     gcif(n) = fsafe * ( abs( f(n+1) - f(n) ) / f(n) ) / ( r(n)**p - 1.0 )
   ENDDO

   WRITE(*,*) ' '
   WRITE(*,*) 'Grid Convergence Index on fine grids. Uses p from above.'
   WRITE(*,*) 'Factor of Safety = ', fsafe
   WRITE(*,*) ' '
   WRITE(*,*) '  Grid     Refinement            '
   WRITE(*,*) '  Step      Ratio, r       GCI(%)'
   DO n = 1, nd-1
     WRITE(*,'(2x,i2,i3,2(f14.6))') n, n+1, r(n), gcif(n)*100.0
   ENDDO

!..Examine if asymptotic range has been achieved by examining ratio 
!..of Eqn. 5.10.5.2 of Roache's book is one.

   IF ( nd .ge. 3 ) then

     WRITE(*,*) ' '
     WRITE(*,*) 'Checking for asymptotic range using Eqn. 5.10.5.2.'
     WRITE(*,*) 'A ratio of 1.0 indicates asymptotic range.'
     WRITE(*,*) ' '
     WRITE(*,*) ' Grid Range    Ratio'
     DO n = 1, nd-2
       ratio = r(n)**p * gcif(n) / gcif(n+1)
       WRITE(*,'(2x,i1,i1,i2,i1,f14.6)') n, n+1, n+1, n+2, ratio
     ENDDO

   ENDIF


!..Write Trailer.

   WRITE(*,*) ' '
   WRITE(*,*) '--- End of VERIFY ---'
   WRITE(*,*) ' '
 
!--------------------------------------------------------------------------

    END PROGRAM verify

