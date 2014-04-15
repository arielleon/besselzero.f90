program besselzero
    !
    !  This program is to print the first k values of zeros of &
    !  Bessel function j_n(x) (kind = 1) or y_n(x) (kind = 2)" for n >= 0. 
    ! 
    

    use const
    implicit none
    integer :: n, k, kind
    integer :: kk
    real(kind=8), allocatable, dimension(:) :: x
    interface
        subroutine sub_besselzero(n, k, kind, x)
             integer, intent(in) :: n, k, kind
             real(kind=8), dimension(:), intent(out) :: x    
        end subroutine sub_besselzero
    end interface
 
    print *, "This program is to print the first k values of zeros of &
             Bessel function j_n(x) (kind = 1) or y_n(x) (kind = 2)"
    print *, "Please type in input:"
    print *, "The integer n= "
    read *, n
    print *, "The integer k="
    read *, k
    print *, "kind = (Either 1 or 2)"
    read *, kind

    
    allocate(x(k))
    call sub_besselzero(n, k, kind,  x)

    do kk = 1, k
       print *, kk, x(kk)
    end do
    print *,"END"
    
end program besselzero



subroutine sub_besselzero(n, k, kind, x)
    use const
    implicit none

    integer, intent(in) :: n, k, kind
    real(kind=8), dimension(:), intent(out) :: x
    real(kind=8), external :: findzero
  
    INTERFACE
       SUBROUTINE  Sort(x, Size)
           REAL(KIND=8), DIMENSION(1:), INTENT(INOUT) :: x
           INTEGER, INTENT(IN)                   :: Size
       END SUBROUTINE SORT
    END INTERFACE

    integer :: k3
    real(kind=8) :: x0
    real(kind=8), allocatable, dimension(:) :: all_x
    real(kind=8), allocatable, dimension(:) :: all_dx
    real(kind=8), allocatable, dimension(:) :: all_xtemp
    real(kind=8), allocatable, dimension(:) :: all_xout
    
    integer :: kk
    integer :: j

    k3 = 3*k

    ! Same as all_xmat = zeros(k3, 1)
    allocate (all_x(k3))
    do kk =1,k3

       ! Initial guess of zeros
       x0 = 1+sqrt(2.0)+(kk-1)*pi+n+n**0.4
       ! Do Halley's method
       all_x(kk) = findzero(n, x0, kind)
       if ( dabs(all_x(kk)+1) .le. 1e-5) then
         print *, "Bad guess."
       end if      
    end do

    ! all_x = sort(all_x)
    
    call sort(all_x, k3)

    ! dx = [1; abs(diff(x))]
    ! x = x(dx>1e-8)
    ! x = x(1:k)
    
    allocate (all_dx(k3))
   
    all_dx(1) = 1
    do kk = 2, k3
       all_dx(kk) = dabs(all_x(kk) - all_x(kk-1))
    end do
    
    allocate (all_xtemp(k3))
    j = 1
    do kk = 1, k3
       if (all_dx(kk) .gt. 1e-8) then
          all_xtemp(j) = all_x(kk)
          j = j+1
       end if
    end do

    allocate(all_xout(k))
    
    do kk = 1, k
       all_xout(kk) = all_xtemp(kk)
    end do
    
    x = all_xout
    
    
    deallocate (all_x, all_dx, all_xtemp, all_xout)
   

    

end subroutine sub_besselzero

real(kind=8)  function findzero(n, x0, kind)
     use const
     implicit none
   
     integer, intent(in) :: n, kind
     real(kind=8), intent(in) :: x0
    
     real(kind=8) :: x1 
     integer :: n1, n2
     real(kind=8) :: tol = 1e-12     ! Tolerance
     integer :: maxit = 100;     ! Maximum number of times to iterate
     real(kind=8) :: err         ! Initial error
     integer :: iter
     real(kind=8) :: a, b, x02, x
   
     n1 = n+1
     n2 = n*n     
    
     !!! Att: always do the initial explicitly
     iter = 0
     err = 1
     x1 = x0
     do while (dabs(err) .gt. tol .and. iter .le. maxit)
        select case (kind)
        case (1)
             a = bessel_jn(n, x1)
             b = bessel_jn(n1, x1)
        case (2) 
             a = bessel_yn(n, x1)
             b = bessel_yn(n1, x1)
        end select
     
        x02 = x1*x1
       
        err = 2*a*x1*(n*a-b*x1)/(2*b*b*x02-a*b*x1*(4*n+1)+(n*n1+x02)*a*a)

        x = x1-err
        x1 = x
        iter = iter+1
                  
     end do

     if ( iter .gt. (maxit-1)) then
        x = -1
     end if

     findzero = x

end function findzero

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      REAL(KIND=8), DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      REAL(KIND=8)                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, End		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location        	! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      REAL(KIND=8), INTENT(INOUT) :: a, b
      REAL(KIND=8)                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      REAL(KIND=8), DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location
      INTERFACE
            INTEGER FUNCTION  FindMinimum(x, Start, End)
                    REAL(KIND=8), DIMENSION(1:), INTENT(IN) :: x
                    INTEGER, INTENT(IN)                :: Start, End
                    REAL(KIND=8)                            :: Minimum
                    INTEGER                            :: Location
                    INTEGER                            :: i
            END FUNCTION FindMinimum
      END INTERFACE
      DO i = 1, Size-1			! except for the last
         Location = FindMinimum(x, i, Size)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
