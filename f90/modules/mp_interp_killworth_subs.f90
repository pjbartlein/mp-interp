module mp_interp_killworth_subs
! subroutines for implementing the Killworth (1996) "adjusted-input linear" mean-preserving interpolation
    
! Killworth, P.D. (1996) Time interpolation of forcing fields in ocean models
! J. Physical Oceanography 26:136-143
    
implicit none

    integer             :: debug_unit=10
    logical             :: debug_write = .false.  
    
contains
    
subroutine mp_interp_killworth(n_outer, n_inner, nctrl, ym, yfill, x_ctrl, nsubint, &
    npad, no_negatives, match_mean, tol, ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

    use pseudo_daily_interp_subs

    implicit none
    
    integer(4), intent(in)      :: nctrl, ntargs                    ! number of control and target points
    integer(4), intent(in)      :: n_outer, n_inner                 ! number of outer and inner intervals (e.g. months, years)
    real(8), intent(in)         :: ym(nctrl)                        ! input
    real(8), intent(in)         :: yfill                            ! fill value
    real(8), intent(in)         :: x_ctrl(nctrl)                    ! x_cntrl (used for debugging)
    integer(4), intent(in)      :: nsubint(nctrl)                   ! number of subintervals
    integer(4), intent(in)      :: npad                             ! number of padding each end
    logical, intent(in)         :: no_negatives, match_mean         ! logical control variables
    real(8), intent(in)         :: tol                              ! tolerance for enforce_mean()
    real(8), intent(in)         :: x_targ(ntargs)                   ! x_targ (for debuggin)
    integer(4), intent(in)      :: max_nctrl_in                     ! max number of "inner" intervals in an "outer" interval
    integer(4), intent(in)      :: max_ntargs_in                    ! max number of subintervals in an "outer" interval 
    
    real(8), intent(inout)        :: y_int(ntargs)                  ! interpolated values
    real(8), intent(inout)        :: ym_int(nctrl)                  ! mean of interpolated values
                                                                    ! inner interval = e.g. months, outer interval = e.g. year
    real(8)             :: ym_in(max_nctrl_in), x_ctrl_in(max_nctrl_in)
    real(8)             :: x_targ_in(max_ntargs_in)
    real(8)             :: y_int_out(max_ntargs_in), ym_int_out(max_ntargs_in)
    real(8)             :: rmse
    
    integer(4)          :: nsubint_in(max_nctrl_in)
    
    integer(4)          :: n, nn
    integer(4)          :: beg_inner, end_inner, beg_ctrl, end_ctrl, nctrl_in, beg_targ, end_targ, ntarg_in, beg_int, end_int, nint_out
    
    if (debug_write) write (debug_unit,'(a)') "In mp_interp_killworth()"
    
    ! loop over number of years
    nn = 1 ! starting pseudo-daily value
    do n = 1, n_outer
        
        ! range of inner intervals (e.g. months)
        beg_inner = (n - 1) * n_inner + 1 
        end_inner = n * n_inner 
        
        ! range of control values (padded)
        beg_ctrl = beg_inner - npad
        end_ctrl = end_inner + npad

        if (n .eq. 1) beg_ctrl = 1
        if (n .eq. n_outer) end_ctrl = n * n_inner
    
        nctrl_in = end_ctrl - beg_ctrl + 1
        
        ! range of target values, including padding
        beg_targ = sum(nsubint(1:(beg_ctrl - 1))) + 1
        end_targ = beg_targ + sum(nsubint(beg_ctrl:end_ctrl)) - 1
        
        ntarg_in = end_targ - beg_targ + 1
        
        ! range of interpolated values to save
        beg_int = nn 
        end_int = nn + sum(nsubint(beg_inner:end_inner)) - 1
        
        nint_out = end_int - beg_int + 1
        
        ! update nn for next year
        nn = end_int + 1
        !write (*, '(5i6)') n, beg_int, end_int, nint_out, nn
        
        if (debug_write) write (debug_unit, '(a)') " "
        if (debug_write) write (debug_unit, '("n,beg_inner,end_inner,n_inner): ",4i9)') &
            n,beg_inner, end_inner, end_inner - beg_inner + 1
        if (debug_write) write (debug_unit, '("  beg_ctrl,end_ctrl,nctrl_in,x_ctrl(beg_ctrl),x_ctrl(end_ctrl): ",3i9,2f16.9)') &
                beg_ctrl,end_ctrl,nctrl_in,x_ctrl(beg_ctrl),x_ctrl(end_ctrl)
        if (debug_write) write (debug_unit, '("  beg_targ,end_targ,ntarg_in,x_targ(beg_targ),x_targ(end_targ): ",3i9,2f16.9)') &
                beg_targ,end_targ,ntarg_in,x_targ(beg_targ),x_targ(end_targ)
        if (debug_write) write (debug_unit, '("  beg_int,end_int,nint_in,x_targ(beg_int),x_targ(end_int): ",3i9,2f16.9)') &
                beg_int,end_int,nint_out,x_targ(beg_int),x_targ(end_int)
    
        ! load temporary variables
        
        ym_in = 0.0d0; x_ctrl_in = 0.0d0; nsubint_in = 0; y_int_out = 0.0

        ym_in(1:nctrl_in) = ym(beg_ctrl:end_ctrl)
        x_ctrl_in(1:nctrl_in) = x_ctrl(beg_ctrl:end_ctrl)
        nsubint_in(1:nctrl_in) = nsubint(beg_ctrl:end_ctrl)   
        x_targ_in(1:ntarg_in) = x_targ(beg_targ:end_targ)
        
        if (debug_write) write (debug_unit, '(a)') " ym_in(1:nctrl_in)"
        if (debug_write) write (debug_unit, '(16f16.9)') ym_in(1:nctrl_in)
        if (debug_write) write (debug_unit, '(a)') " x_ctrl_in(1:nctrl_in)"
        if (debug_write) write (debug_unit, '(16f16.9)') x_ctrl_in(1:nctrl_in)
        if (debug_write) write (debug_unit, '(a)') " nsubint_in(1:nctrl_in)"
        if (debug_write) write (debug_unit, '(16i16)') nsubint_in(1:nctrl_in)
        if (debug_write) write (debug_unit, '(a)') " x_targ_in(1), x_targ_in(ntarg_in)"
        if (debug_write) write (debug_unit, '(16f16.9)') x_targ_in(1), x_targ_in(ntarg_in)
        
        ! if any missing values, set output to yfill, otherwise call hz_int()
        if (any(ym_in(1:nctrl_in) .eq. yfill)) then
        
            if (debug_write) write (debug_unit, '("missing values, interval ", i8)') n
            y_int(beg_targ:end_targ) = yfill
            
        else
        
            if (debug_write) write (debug_unit, '("call hz_int(), interval ", i8)') n
            
            call ainv_int(nctrl_in, ym_in, yfill, x_ctrl_in, nsubint_in, ntarg_in, x_targ_in, y_int_out, ym_int_out)
            
            if (debug_write) write (debug_unit, '(a)') "back from hz_int()"
        
            if (debug_write) write (debug_unit,'(10f10.6)') y_int_out
            if (debug_write) write (debug_unit,'(16f10.6)') ym_in(1:nctrl_in)
            if (debug_write) write (debug_unit,'(16f10.6)') ym_int_out(1:nctrl_in)
            if (debug_write) write (debug_unit,'(16f10.6)') ym_in(1:nctrl_in) - ym_int_out(1:nctrl_in)
            if (debug_write) then
                call interp_stat(nctrl_in, ym_in, ym_int_out, rmse)
                write (debug_unit,'("rmse: ", g12.4)') rmse
            end if
        
            y_int(beg_int:end_int) = y_int_out((beg_int - beg_targ + 1):(beg_int - beg_targ + nint_out))

            if (debug_write) then
                write (debug_unit, '(a)') "n,beg_targ, beg_int, end_int, nint_out, beg_int - beg_targ + 1, beg_int - beg_targ + nint_out"
                write (debug_unit, '(8i9)') n,beg_targ, beg_int, end_int, nint_out, beg_int - beg_targ + 1, beg_int - beg_targ + nint_out
                write (debug_unit, '(a)') "y_int_out((beg_int - beg_targ + 1):(beg_int - beg_targ + nint_out))"
                write (debug_unit,'(10f10.6)') y_int_out((beg_int - beg_targ + 1):(beg_int - beg_targ + nint_out))
                write (debug_unit, '(a)') "y_int(beg_int:end_int)"
                write (debug_unit,'(10f10.6)') y_int(beg_int:end_int)
            end if
                
        end if

    end do
    
    ! no negatives (also check for NaN's)
    if (no_negatives) then
        do n = 1, ntargs
            if (y_int(n) .lt. 0.0d0) y_int(n) = 0.0d0
            if (isnan(y_int(n))) y_int(n) = 0.0d0
        end do
    end if
        
    if (debug_write) write (debug_unit,'(a)') "y_int_before_enforce_mean"
    if (debug_write) write (debug_unit,'(10f8.2)') y_int
        
    if (debug_write) write (debug_unit,*) nctrl, ntargs
    if (match_mean) then
        if (debug_write) write (debug_unit, '(a)') "call enforce_mean() "
        call enforce_mean(nctrl, ntargs, nsubint, tol, ym, y_int, yfill)
    end if
    
    call interval_mean(nctrl, nsubint, ntargs, y_int, yfill, ym_int)
    
    if (debug_write) write (debug_unit,'(a)') "y_int_after_enforce_mean"
    if (debug_write) write (debug_unit,'(10f8.2)') y_int
    if (debug_write) write (debug_unit,'(12f8.3)') ym_int

    
end subroutine mp_interp_killworth
    
subroutine ainv_int(nctrl, ym, ymiss, x_ctrl, nsubint, ntargs, x_targ, y_int, ym_int)
! manages Killworth-style "A-inv" pseudo-daily interpolation

    use pseudo_daily_interp_subs
    use lapack95    ! Intel MKL on Windows
    !use f95_lapack  ! https://masuday.github.io/fortran_tutorial/external.html#blas-and-lapack + lapack from Homebrew on mac
    
    implicit none
    
    integer(4), intent(in)          :: nctrl, ntargs
    integer(4), intent(in)          :: nsubint(nctrl)                   ! interval lengths
    real(8), intent(in)             :: ym(nctrl)                        ! input y
    real(8), intent(in)             :: ymiss                            ! missing / fill value
    real(8), intent(in)             :: x_ctrl(nctrl)                    ! input x
    real(8), intent(in)             :: x_targ(ntargs)                   ! target values
    real(8), intent(out)            :: y_int(ntargs)                    ! interpolated values
    real(8), intent(out)            :: ym_int(nctrl)                    ! means of interpolated values
    
    ! local variables
    real(8)                         :: A(nctrl, nctrl)                  ! matrix A
    real(8)                         :: A_inv(nctrl,nctrl)
    real(8)                         :: e(nctrl), f(nctrl), g(nctrl)     ! A elements
    real(8)                         :: am_int(nctrl), ymp(nctrl)
    integer(4)                      :: ipiv(nctrl)
    real(8)                         :: work(nctrl)
    
    integer(4)                      :: m, m1, info
    
    if (debug_write) write (debug_unit, '(a)') "In ainv_interp"
    if (debug_write) write (debug_unit, '(12f9.2)') ym
    if (debug_write) write (debug_unit, '(12f9.2)') x_ctrl
    if (debug_write) write (debug_unit, '(a)') "x_targ"
    if (debug_write) write (debug_unit, '(12f9.2)') x_targ
    
    call dayinterp(nctrl, ntargs, nsubint, ym, y_int)
    call interval_mean(nctrl, nsubint, ntargs, y_int, ymiss, ym_int)
    
    if (debug_write) write (debug_unit, '(a)') "ym, ym_int, ym_int - ym"
    if (debug_write) write (debug_unit,'(12f9.4)') ym
    if (debug_write) write (debug_unit,'(12f9.4)') ym_int
    if (debug_write) write (debug_unit,'(12f9.4)') ym_int - ym
    
    ! elements
    e = 0.0d0; f = 0.0d0; g = 0.0d0
        e(1) = nsubint(1) / (4.0d0 * (nsubint(nctrl) + nsubint(1)))
    do m = 2, nctrl
        e(m) = nsubint(m) / (4.0d0 * (nsubint(m - 1) + nsubint(m)))
    end do
        g(nctrl) = nsubint(nctrl) / (4.0d0 * (nsubint(nctrl) + nsubint(1)))
    do m = 1, nctrl - 1
        g(m) = nsubint(m) / (4.0d0 * (nsubint(m) + nsubint(m+1)))
    end do
    do m = 1, nctrl
        f(m) = 1.0d0 - e(m) - g(m)
    end do
    if (debug_write) write (debug_unit, '(a)') "e, f, g"
    if (debug_write) write (debug_unit, '(12f9.4)') e, f, g

    ! define A
    A = 0.0d0
    do m = 2, nctrl
        A(m, m - 1) = e(m)
        A(m - 1, m) = g(m - 1)
        A(m, m) = f(m)
    end do
        A(1, 1) = f(1)
        A(nctrl, 1) = g(nctrl)
        A(1, nctrl) = e(1)
    if (debug_write) write (debug_unit, '(a)') "A"
    do m1 = 1, nctrl
        if (debug_write) write (debug_unit, '(12f9.4)') (A(m1, m), m = 1, nctrl)
    end do
    
    am_int = matmul(A, ym)
    if (debug_write) write (debug_unit, '(a)') "am_int"
    if (debug_write) write (debug_unit, '(12f9.4)') am_int
    
    ! inverse
    A_inv = A
    call DGETRF(nctrl, nctrl, A_inv, nctrl, ipiv, info)
    call DGETRI(nctrl, A_inv, nctrl, ipiv, work, nctrl, info)
    
    if (debug_write) write (debug_unit, '(a)') "A_inv"
    do m1 = 1, nctrl
        if (debug_write) write (debug_unit, '(12f9.4)') (A_inv(m1, m), m = 1, nctrl)
    end do
    
    !! test the inverse
    !Id = 0.0d0
    !Id = matmul(A_inv, A)
    !if (debug_write) write (debug_unit, '(a)') "I"
    !do m1 = 1, nctrl
    !    if (debug_write) write (debug_unit, '(12f9.4)') (Id(m1, m), m = 1, nctrl)
    !end do
    
    ymp = matmul(A_inv, ym)
    if (debug_write) write (debug_unit, '(a)') "ymp"
    if (debug_write) write (debug_unit, '(12f9.4)') ymp
    
    call dayinterp(nctrl, ntargs, nsubint, ymp, y_int)
    call interval_mean(nctrl, nsubint, ntargs, y_int, ymiss, ym_int)
    
    if (debug_write) write (debug_unit, '(a)') "ym, ymp_int, ymp_int - ym"
    if (debug_write) write (debug_unit,'(12f9.4)') ym
    if (debug_write) write (debug_unit,'(12f9.4)') ym_int
    if (debug_write) write (debug_unit,'(12f9.4)') ym_int - ym
    
end subroutine ainv_int     
    
    
end module mp_interp_killworth_subs
