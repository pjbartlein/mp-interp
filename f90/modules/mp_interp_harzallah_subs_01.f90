module mp_interp_harzallah_subs
! subroutines for implementing the Harzallah (1995) "iterative spline" mean-preserving interpolation
    
! Harzallah, A. (1995) The interpolation of data series using a constrained iterating technique
! Monthly Weather Review 123:2251-2254.
    
! This implementation provides three versions of spline fitting.
    
implicit none

    integer             :: debug_unit=10
    logical             :: debug_write = .true.  
    
    contains
    
 !subroutine month_to_day_ts_hz(spline_case, n_outer, nm, nctrl, ym, yfill, x_ctrl, nsubint, &
 !   no_negatives, match_mean, tol, ntargs, x_targ, y_int, ym_int)   
    
subroutine mp_interp_harzallah(n_outer, n_inner, nctrl, ym, yfill, x_ctrl, nsubint, &
    spline_case, npad, no_negatives, match_mean, tol, ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

    use mean_preserving_subs

    implicit none
    
    integer(4), intent(in)      :: nctrl, ntargs                    ! number of control and target points
    integer(4), intent(in)      :: n_outer, n_inner                 ! number of outer and inner intervals (e.g. months, years)
    real(8), intent(in)         :: ym(nctrl)                        ! input
    real(8), intent(in)         :: yfill                            ! fill value
    real(8), intent(in)         :: x_ctrl(nctrl)                    ! x_cntrl (used for debugging)
    integer(4), intent(in)      :: nsubint(nctrl)                   ! number of subintervals
    integer(4), intent(in)      :: spline_case                      ! selects spline-fitting procedure
    integer(4), intent(in)      :: npad                             ! number of padding each end
    logical, intent(in)         :: no_negatives, match_mean         ! logical control variables
    real(8), intent(in)         :: tol                              ! tolerance for enforce_mean()
    real(8), intent(in)         :: x_targ(ntargs)                   ! x_targ (for debuggin)
    integer(4), intent(in)      :: max_nctrl_in                     ! max number of "inner" intervals in an "outer" interval
    integer(4), intent(in)      :: max_ntargs_in                    ! max number of subintervals in an "outer" interval 
    
    real(8), intent(out)        :: y_int(ntargs)                    ! interpolated values
    real(8), intent(out)        :: ym_int(nctrl)                    ! mean of interpolated values
                                                                    ! inner interval = e.g. months, outer interval = e.g. year
    real(8)             :: ym_in(max_nctrl_in), x_ctrl_in(max_nctrl_in)
    real(8)             :: x_targ_in(max_ntargs_in)
    real(8)             :: y_int_out(max_ntargs_in), ym_int_out(max_ntargs_in)
    real(8)             :: rmse
    integer(4)          :: nsubint_in(max_nctrl_in)
    
    integer(4)          :: i, n
    integer(4)          :: nn, beg_inner, end_inner, beg_subint, end_subint, n_inner_in, ntargs_out
    integer(4)          :: nbeg, nend, nsubint_padded, extra_nsubint
    
    if (debug_write) write (debug_unit,'(a)') "In mp_interp_harzallah()"

    ! loop over number of years
    nn = 1
    do n = 1, n_outer
    
        nbeg = npad; nend = npad
    
        if (n .eq. 1) nbeg = 0
        if (n .eq. n_outer) nend = 0
    
        beg_inner = (n - 1) * n_inner + 1
        end_inner = beg_inner + (n_inner - 1)
        ntargs_out = sum(nsubint((beg_inner - nbeg):(end_inner + nend)))
        beg_subint = nn
        end_subint = beg_subint - 1 + sum(nsubint(beg_inner:end_inner))
        nn = end_subint + 1
        nsubint_padded = (end_subint - beg_subint + 1)
    
        n_inner_in = (end_inner + nend) - (beg_inner - nbeg) + 1
        if (n .gt. 1) then
            extra_nsubint = sum(nsubint((beg_inner - nbeg):(beg_inner))) - nsubint(beg_inner)
        else 
            extra_nsubint = 0
        end if

        if (debug_write) write (debug_unit, '("n,beg_inner,end_inner,nbeg,nend,beg_inner-nbeg,end_inner+nend,: ", 7i5)') & 
            n, beg_inner, end_inner, nbeg, nend, beg_inner - nbeg, end_inner + nend
        if (debug_write) write (debug_unit, '("n_inner_in, ntargs_out, beg_subint, end_subint, nsubint_padded, extra_nsubint: ", 2i5, 2i8, 2i4)') &
            n_inner_in, ntargs_out, beg_subint, end_subint, nsubint_padded, extra_nsubint
    
        ! load temporary variables
    
        ym_in = 0.0d0; x_ctrl_in = 0.0d0; nsubint_in = 0
        ym_in(1:n_inner_in) = ym((beg_inner - nbeg):(end_inner + nend))
        x_ctrl_in(1:n_inner_in) = x_ctrl((beg_inner - nbeg):(end_inner + nend))
        nsubint_in(1:n_inner_in) = nsubint((beg_inner - nbeg):(end_inner + nend))

        do i = 1, ntargs_out
            x_targ_in(i) = dble(beg_subint - extra_nsubint) + dble(i - 1)
            !write (10,*) i, x_targ_in(i)
        end do
    
        if (debug_write) write (debug_unit, '(16g14.6)') ym_in(1:n_inner_in)
        if (debug_write) write (debug_unit, '(16f14.2)') x_ctrl_in(1:n_inner_in)
        if (debug_write) write (debug_unit, '(16i14)') nsubint_in(1:n_inner_in)
        
        ! if any missing values, set output to yfill, otherwise call hz_int()
        if (any(ym_in(1:n_inner_in) .eq. yfill)) then
        
            !if (debug_write) write (debug_unit, '("missing values, interval ", i8)') n
            ym_int(beg_inner:end_inner) = yfill
            y_int(beg_subint:end_subint) = yfill
            
        else

            !if (debug_write) write (debug_unit, '("call hz_int(), interval ", i8)') n
            call hz_int(spline_case, n_inner_in, ym_in, yfill, x_ctrl_in, nsubint_in, ntargs_out, x_targ_in, y_int_out, ym_int_out)
    
            ym_int(beg_inner:end_inner) = ym_int_out((nbeg + 1):(ntargs-nend))
            y_int(beg_subint:end_subint) = y_int_out((extra_nsubint+1):(extra_nsubint + nsubint_padded - 1))          
                
        end if

    end do
    
    ! no negatives
    if (no_negatives) then
        do n = 1, ntargs
            if (y_int(n) .lt. 0.0d0) y_int(n) = 0.0d0
        end do
    end if
    
    ! check for NaNs
    if (no_negatives) then
        do n = 1, ntargs
            !if (debug_write) write (debug_unit,'("n, y_int(n): ", i6, g12.4, 1x, l1)') n, y_int(n), isnan(y_int(n))
            if (isnan(y_int(n))) y_int(n) = 0.0d0
        end do
    end if
        
    if (debug_write) write (debug_unit,'(a)') "y_int_before_enforce_mean"
    if (debug_write) write (debug_unit,'(10f8.2)') y_int
        
    if (debug_write) write (debug_unit,*) nctrl, ntargs
    if (match_mean) then
        if (debug_write) write (debug_unit, '(a)') "call enforce_mean() "
        call enforce_mean(no_negatives, nctrl, ntargs, nsubint, tol, ym, y_int, yfill)
    end if
    
    call interval_mean(nctrl, nsubint, ntargs, y_int, yfill, ym_int)
    
    if (debug_write) write (debug_unit,'(a)') "y_int_after_enforce_mean"
    if (debug_write) write (debug_unit,'(10f8.2)') y_int
    if (debug_write) write (debug_unit,'(12f8.3)') ym_int
    
end subroutine mp_interp_harzallah
    
subroutine hz_int(spline_case, nctrl, ym, ymiss, x_ctrl, nsubint, ntargs, x_targ, y_int, ym_int)
! manages Harzalla spline pseudo-daily interpolation

    use mean_preserving_subs
    
    implicit none
    
    integer(4), intent(in)          :: nctrl, ntargs, spline_case
    integer(4), intent(in)          :: nsubint(nctrl)                   ! interval lengths
    real(8), intent(in)             :: ym(nctrl)                        ! input y
    real(8), intent(in)             :: ymiss                            ! missing / fill value
    real(8), intent(in)             :: x_ctrl(nctrl)                    ! input x
    real(8), intent(in)             :: x_targ(ntargs)                   ! target values
    real(8), intent(out)            :: y_int(ntargs)                    ! interpolated values
    real(8), intent(out)            :: ym_int(nctrl)                    ! means of interpolated values
    
    ! local variables
    integer(4)                      :: ibcbeg = 0, ibcend = 0           ! boundary condition flags
    real(8)                         :: ybcbeg = 0.0d0, ybcend = 0.0d0   ! boundary conditions
    real(8)                         :: ypp(nctrl)                       ! derivatives
    real(8)                         :: ypval(ntargs), yppval(ntargs)    ! first and second derivatives
    
    integer(4), parameter           :: max_iter = 100
    real(8)                         :: tol = 0.01                       ! tolerance
    real(8)                         :: resid(nctrl)
    real(8)                         :: ym_mat(max_iter + 1, nctrl)       ! matrix
    real(8)                         :: y_int_mat(max_iter, ntargs)
    
    integer(4)                      :: i, k
    integer(4)                      :: np = 3, err
    
    if (debug_write) write (debug_unit, '(a)') "In hz_interp ============="
    if (debug_write) write (debug_unit, *) nctrl, ntargs, spline_case
    if (debug_write) write (debug_unit, '(12f8.2)') ym
    if (debug_write) write (debug_unit, '(12f11.2)') x_ctrl
    if (debug_write) write (debug_unit, '(12f11.2)') x_targ
    
    ym_mat = 0.0d0; y_int_mat = 0.0d0
    ym_mat(1, :) = ym(:)
    
    do k = 1, max_iter

        if (debug_write) write (debug_unit, '("ym_mat", i4)') k
        if (debug_write) write (debug_unit,'(12g12.4)') ym_mat(k, :)
        if (debug_write) write (debug_unit,*) spline_case
        
        select case(spline_case)
            case (1)
                ! Burkhardt
                call spline_cubic_set ( nctrl, x_ctrl, ym_mat(k, :), ibcbeg, ybcbeg, ibcend, ybcend, ypp )
                if (debug_write) write (debug_unit,'(12g16.8)') ypp
                do i = 1, ntargs
                    call spline_cubic_val ( nctrl, x_ctrl, ym_mat(k, :), ypp, x_targ(i), y_int_mat(k,i), ypval, yppval )
                    if (debug_write) write (debug_unit,'(i4, 4g16.8)') i, x_targ(i), y_int_mat(k,i) !, ypval, yppval
                end do     
            case (2)
                ! Burkhardt
                call spline_pchip_set ( nctrl, x_ctrl, ym_mat(k, :), ypp )
                if (debug_write) write (debug_unit,'(12g16.8)') ypp      
                call spline_pchip_val ( nctrl, x_ctrl, ym_mat(k, :), ypp, ntargs, x_targ, y_int_mat(k, :) )
            case (3)
                ! Akima ACM 697   
                call uvip3p (np, nctrl, x_ctrl, ym_mat(k, :), ntargs, x_targ, y_int_mat(k, :), err)
            case default
                stop "spline_case"
        end select
        
        if (debug_write) write (debug_unit, '("y int_mat", i4)') k
        if (debug_write) write (debug_unit,'(10g12.4)') y_int_mat(k, :)
        
        ! interval mean values
        if (debug_write) write (debug_unit, '("calling interval_mean() from hz_int(), k, nctrl, ntargs: ", 3i8)') k, nctrl, ntargs
        call interval_mean(nctrl, nsubint, ntargs, y_int_mat(k, :), ymiss, ym_int)
        if (debug_write) write (debug_unit, '("ym_int", i4)') k
        if (debug_write) write (debug_unit,'(12g12.4)') ym_int
     
        ! residuals
        resid = ym_mat(k, :) - ym_int
        if (debug_write) write (debug_unit, '("residuals", i4)') k
        if (debug_write) write (debug_unit,'(12g12.4)') resid
              
        if (debug_write) write (debug_unit,'("max resid: ", g12.4)') maxval(resid)
        if (maxval(resid) .lt. tol) exit
        
        ym_mat(k + 1, :) = resid(:)
        if (debug_write) write (debug_unit, '("ym_mat new", i4)') k
        if (debug_write) write (debug_unit,'(12g12.4)') ym_mat(k+1, :)
     
    end do
    
    if (debug_write) write (debug_unit,'("k (out): ", i4)') k
    
    ! sum over iterations  -- replace by updating?
    y_int = 0.0d0
    do k = 1, max_iter
        do i = 1, ntargs
            y_int(i) = y_int(i) + y_int_mat(k, i)
        end do
    end do
    
    if (debug_write) write (debug_unit, '("calling interval_mean() from hz_int() end, nctrl, ntargs: ", 3i8)') nctrl, ntargs    
    call interval_mean(nctrl, nsubint, ntargs, y_int, ymiss, ym_int)

    if (debug_write) write (debug_unit,'(12f10.2)') ym
    if (debug_write) write (debug_unit,'(12f10.2)') ym_int
    if (debug_write) write (debug_unit,'(12f10.5)') ym_int - ym
    
end subroutine hz_int

end module mp_interp_harzallah_subs