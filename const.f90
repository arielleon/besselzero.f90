module const
    implicit none

    ! Define a real kind type q with at least 6 decimal 
    ! digits and an exponent range from 10**30 to 10**(-30)
    integer, parameter :: q = selected_real_kind(p = 6, r = 30)

    ! Define pi
    real(kind=q), parameter :: pi = 3.1415926536_q

end module const
