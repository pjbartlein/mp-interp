module mp_interp_rymes_and_meyer_subs
! subroutines for implementing the Rymes abd Meyer (2001) "iterative smoothing" mean-preserving interpolation
    
! Rymes, M.D. and D.R. Meyers (2001) Mean preserving alogrithm for smoothing interpolating averaged data.
! Solar Energy 71:225-231.
    
implicit none

    integer             :: debug_unit=10
    logical             :: debug_write = .true.  
    
    contains
    
!subroutine month_to_day_ts_rm(ny, nm, nctrl, ym, x_ctrl, nsubint, lowerbound, lower, upperbound, upper, &
!    no_negatives, match_mean, ntargs, x_targ, y_int, ym_int)
    
subroutine mp_interp_rymes_and_meyer(n_outer, n_inner, nctrl, ym, yfill, x_ctrl, nsubint,  &
    lowerbound, lower, upperbound, upper, npad, no_negatives, match_mean, tol, & 
    ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

    use mean_preserving_subs

    implicit none
    
    integer(4), intent(in)      :: nctrl, ntargs                    ! number of control and target points
    integer(4), intent(in)      :: n_outer, n_inner                 ! number of outer and inner intervals (e.g. months, years)
    real(8), intent(in)         :: ym(nctrl)                        ! input
    real(8), intent(in)         :: yfill                            ! fill value
    real(8), intent(in)         :: x_ctrl(nctrl)                    ! x_cntrl (used for debugging)
    integer(4), intent(in)      :: nsubint(nctrl)                   ! number of subintervals
    logical, intent(in)         :: lowerbound, upperbound           ! contrain interpolation?
    real(8), intent(in)         :: lower(nctrl), upper(nctrl)       ! lower and upper bounds
    integer(4), intent(in)      :: npad                             ! number of padding each end
    logical, intent(in)         :: no_negatives, match_mean         ! logical control variables
    real(8), intent(in)         :: tol                              ! tolerance for enforce_mean()
    real(8), intent(in)         :: x_targ(ntargs)                   ! x_targ 
    integer(4), intent(in)      :: max_nctrl_in                     ! max number of "inner" intervals in an "outer" interval
    integer(4), intent(in)      :: max_ntargs_in                    ! max number of subintervals in an "outer" interval 
    
    real(8), intent(out)        :: y_int(ntargs)                    ! interpolated values
    real(8), intent(out)        :: ym_int(nctrl)                    ! mean of interpolated values
                                                                    ! inner interval = e.g. months, outer interval = e.g. year
    
    real(8)                     :: ybeg, yend                       ! initial endpoint values
    integer(4), parameter       :: maxit = 50                       ! max number of iterations
    integer(4), parameter       :: init_type = 2                    ! initialization type, 1 = spread, 2 = linear interp
    
    real(8)             :: ym_in(max_nctrl_in), x_ctrl_in(max_nctrl_in)
    real(8)             :: x_targ_in(max_ntargs_in)
    real(8)             :: y_int_out(max_ntargs_in), ym_int_out(max_ntargs_in)
    real(8)             :: lower_in(max_nctrl_in), upper_in(max_nctrl_in)
    real(8)             :: rmse
    integer(4)          :: nsubint_in(max_nctrl_in)
    
    integer(4)          :: i, n
    integer(4)          :: nn, beg_inner, end_inner, beg_subint, end_subint, n_inner_in, ntargs_out
    integer(4)          :: nbeg, nend, nsubint_padded, extra_nsubint
    
    if (debug_write) write (debug_unit,'(a)') "In mp_interp_rymes-and-meyer()"

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
    
        ! initial values for beginning and ending y_int's 
    
        ybeg = ym(end_inner + nend) 
        yend = ym(beg_inner - nbeg) 
        if (debug_write) write (debug_unit, '(2f12.6)') ybeg, yend 
    
        ! load temporary variables
    
        ym_in = 0.0d0; x_ctrl_in = 0.0d0; nsubint_in = 0; lower_in = 0.0d0; upper_in = 0.0d0
        ym_in(1:n_inner_in) = ym((beg_inner - nbeg):(end_inner + nend))
        x_ctrl_in(1:n_inner_in) = x_ctrl((beg_inner - nbeg):(end_inner + nend))
        nsubint_in(1:n_inner_in) = nsubint((beg_inner - nbeg):(end_inner + nend))
        lower_in(1:n_inner_in) = lower((beg_inner - nbeg):(end_inner + nend))
        upper_in(1:n_inner_in) = upper((beg_inner - nbeg):(end_inner + nend))
    
        do i = 1, ntargs_out
            x_targ_in(i) = dble(beg_subint - extra_nsubint) + dble(i - 1)
            !write (10,*) i, x_targ_in(i)
        end do
    
        if (debug_write) write (debug_unit, '(16f8.2)') ym_in(1:n_inner_in)
        if (debug_write) write (debug_unit, '(16f8.2)') x_ctrl_in(1:n_inner_in)
        if (debug_write) write (debug_unit, '(16i8)') nsubint_in(1:n_inner_in)
        
        ! if any missing values, set output to yfill, otherwise call hz_int()
        if (any(ym_in(1:n_inner_in) .eq. yfill)) then
        
            !if (debug_write) write (debug_unit, '("missing values, interval ", i8)') n
            ym_int(beg_inner:end_inner) = yfill
            y_int(beg_subint:end_subint) = yfill
            
        else
            
            !if (debug_write) write (debug_unit, '("call rm_int(), interval ", i8)') n
            call rm_int(n_inner_in, ym_in, yfill, x_ctrl_in, nsubint_in, ybeg, yend, lowerbound, lower_in, upperbound, upper_in, & 
                init_type, ntargs_out, maxit, x_targ_in, y_int_out, ym_int_out)
    
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
    
end subroutine mp_interp_rymes_and_meyer
    
    
subroutine rm_int(nctrl, ym, ymiss, x_ctrl, nsubint, ybeg_in, yend_in, lowerbound, lower, upperbound, upper, &
    init_type, ntargs, maxit, x_targ, y_int, ym_int) 
! pseudo-daily interpolated values
! From Rymes, M.D. and D.R. Myers, 2001. Solar Energy 71:225-231
    
    use mean_preserving_subs

    implicit none

    integer(4), intent(in)              :: nctrl, ntargs
    real(8), intent(in)                 :: ym(nctrl)                    ! input AVG[]
    real(8), intent(in)                 :: ymiss                        ! missing / fill value
    real(8), intent(in)                 :: x_ctrl(nctrl)                ! x_cntrl (used for debugging)
    integer(4), intent(in)              :: nsubint(nctrl)               ! number of subintervals
    real(8), intent(in)                 :: ybeg_in, yend_in             ! initial endpoint values
    logical, intent(in)                 :: lowerbound, upperbound       ! contrain interpolation?
    real(8), intent(in)                 :: lower(nctrl), upper(nctrl)   ! lower and upper bounds
    integer(4), intent(in)              :: init_type                    ! initialization type, 1 = spread, 2 = linear interp
    integer(4), intent(in)              :: maxit                        ! max number of iterations
    real(8), intent(in)                 :: x_targ(ntargs)               ! x_targ (used for debugging)
    real(8), intent(out)                :: y_int(ntargs)                ! interpolated values MN[]
    real(8), intent(out)                :: ym_int(nctrl)                ! mean of interpolated values

    ! internal variables
    real(8), parameter  :: three = 3.0d0
    real(8)             :: ck(nctrl), fk(nctrl), fk_numer(nctrl), fk_denom(nctrl), rmse, ym_current(nctrl)
    real(8)             :: ybeg, yend
    real(8)             :: ck_out(ntargs, nctrl), ym_int_out1(ntargs, nctrl), y_int_out1(0:ntargs, ntargs), rmse_out(0:ntargs)
    integer(4)          :: ibeg(nctrl), iend(nctrl), lflag(nctrl), uflag(nctrl)
    
    integer(4)          :: i, j, k, l, m, iter ! i indexes targets, k indexes controls 
    
    if (debug_write) write (debug_unit, '(a)') "In rm_int"
    if (debug_write) write (debug_unit, *) nctrl, ntargs, maxit
    if (debug_write) write (debug_unit, '(14f9.2)') ym
    if (debug_write) write (debug_unit, '(14i9)') nsubint
    if (debug_write) write (debug_unit,*) lowerbound, upperbound

    !  get beginning and ending subinterval values for each interval
    i = 0
    do k = 1, nctrl
        ibeg(k) = i + 1
        iend(k) = i + nsubint(k)
        i = i + nsubint(k)
    end do
    
    ! alternative initializations of y_int
    select case (init_type)
        case (1)
            ! spread values
            call dayspread(nctrl, ntargs, nsubint, ym, y_int)
        case (2) 
            ! linear interpolation
            call dayinterp(nctrl, ntargs, nsubint, ym, y_int)
        case default
            stop "initialize y_int"
    end select
    y_int_out1(0, :) = y_int(:)
    
    if (debug_write) write (debug_unit,'(a)') "initialize:  y_int"
    if (debug_write) write (debug_unit,'(10f9.2)') y_int
    if (debug_write) write (debug_unit,'(a)') "ibeg, iend"
    if (debug_write) write (debug_unit,'(10i9)') ibeg
    if (debug_write) write (debug_unit,'(10i9)') iend
    
    ! main loop
    
    if (debug_write) write (debug_unit,'("ym: ", 12f12.6)') ym
    
    do iter = 1, maxit !ntargs 
        
        if (debug_write) write (debug_unit, '("iteration: ", i5)') iter
    
        ! update y_int, and get three-point moving avarages of current interpolated values 
    
        do i = 2, (ntargs - 1)
            y_int(i) = (y_int(i-1) + y_int(i) + y_int(i+1)) / three   ! Eqn. 1
        end do
    
        y_int(1)  = (ybeg   + y_int(1)  +  y_int(2)) / three    ! Eqn. 2
        y_int(ntargs) = (y_int(ntargs - 1) + y_int(ntargs) + yend) / three
        
        ! correction factors
        
        ck=0.0d0
        do k = 1, nctrl                                               
            do i = ibeg(k), iend(k)
                ck(k) = ck(k) + (ym(k) - y_int(i))  ! Eqn. 4
            end do
            
            ck(k) = ck(k) / dble(ntargs) 
            
            do i = ibeg(k), iend(k)
                y_int(i) = y_int(i) + ck(k)                        ! apply corrections
            end do

            if (debug_write) write (debug_unit,'("iter,k,ibeg,iend, ck: ", 4i4, g14.6)') iter, k, ibeg(k), iend(k), ck(k)
            ck_out(iter,:) = ck(:)
            
        end do
        
        call interval_mean(nctrl, nsubint, ntargs, y_int, ymiss, ym_current)
          
        ! bounded interpolation
        
        if (lowerbound) then
             
            ym_int = 0.0d0; fk = 0.5d0; fk = 0.0d0; fk_numer = 0.0d0; fk_denom = 0.0d0
            lflag = 0
            
            do k = 1, nctrl                                               
                do i = ibeg(k), iend(k)
                    if (y_int(i) .lt. lower(k)) then 
                        y_int(i) = lower(k)
                        lflag(k) = 1
                    end if
                    ym_int(k) = ym_int(k) + y_int(i)
                    fk_numer(k) = fk_numer(k) + (y_int(i) - ym(k))  ! Eqn. 7
                    fk_denom(k) = fk_denom(k) + (y_int(i) - lower(k))
                end do
                ym_int(k) = ym_int(k) / dble(nsubint(k))    
                
                fk(k) = fk_numer(k)/fk_denom(k)
                if (debug_write) write (debug_unit,'("k,lflag,ym,ym_int,diff,fk: ", 2i5,4g12.4)') &
                    k,lflag(k),ym(k),ym_int(k),ym(k)-ym_int(k),fk(k)
                
                ! if mean interpolated value greater than original mean, adjust downwards
                if (lflag(k) .eq. 1 .and. ym_int(k) .gt. ym(k)) then  
                    do i = ibeg(k), iend(k)
                        y_int(i) = y_int(i) - fk(k) * (ym_int(k) - ym(k))    
                        if (y_int(i) .lt. lower(k)) y_int(i) = lower(k)
                    end do
                end if      
            end do
            
            call interval_mean(nctrl, nsubint, ntargs, y_int, ymiss, ym_current)
            
        end if
        
        if (upperbound) then
             
            ym_int = 0.0d0; fk = 0.5d0; fk = 0.0d0; fk_numer = 0.0d0; fk_denom = 0.0d0
            
            uflag = 0
            
            do k = 1, nctrl                                               
                do i = ibeg(k), iend(k)
                    if (y_int(i) .gt. upper(k)) then 
                        y_int(i) = upper(k)
                        uflag(k) = 1
                    end if
                    ym_int(k) = ym_int(k) + y_int(i)
                    fk_numer(k) = fk_numer(k) + (y_int(i) - ym(k))  ! Eqn. 7
                    fk_denom(k) = fk_denom(k) + (y_int(i) - upper(k))
                end do
                ym_int(k) = ym_int(k) / dble(nsubint(k))

                fk(k) = fk_numer(k)/fk_denom(k)
                if (debug_write) write (debug_unit,'("k,uflag,ym,ym_int,diff,fk: ", 2i5,4g12.4)') & 
                    k,uflag(k),ym(k),ym_int(k),ym(k)-ym_int(k),fk(k)
                  
                ! if mean interpolated value less than original mean, adjust upwards 
                if (uflag(k) .eq. 1 .and. ym_int(k) .lt. ym(k)) then
                    do i = ibeg(k), iend(k)
                        y_int(i) = y_int(i) + fk(k) * (ym(k) - ym_int(k))  
                        if (y_int(i) .gt. upper(k)) y_int(i) = upper(k)
                    end do
                end if  
            end do
            
            call interval_mean(nctrl, nsubint, ntargs, y_int, ymiss, ym_current)
            
        end if
         
        ! update beginning and ending values      
        
        ybeg = y_int(ntargs)
        if (ym(1) .eq. 0.0d0) ybeg = 0.0d0
        yend = y_int(1)
        if (ym(nctrl) .eq. 0.0d0) yend = 0.0d0


        ! means of updated data, and rmse of differences
        
        ym_int = 0.0d0; rmse = 0.0d0
        i = 0
        do k = 1, nctrl
            do j = 1, nsubint(k)
                i = i + 1
                ym_int(k) = ym_int(k) + y_int(i)
                y_int_out1(iter, i) = y_int(i)
            end do
            ym_int(k) = ym_int(k) / dble(nsubint(k))
            rmse = rmse + (ym(k) - ym_int(k))*(ym(k) - ym_int(k))
        end do
        rmse = sqrt(rmse / (nctrl - 1))
        
        ym_int_out1(iter, :) = ym_int(:)
        
        rmse_out(iter) = rmse
        
        !write (*, '("iteration, rmse: ",i4, g14.6)') iter, rmse
        if (debug_write) write (debug_unit, '("iteration, rmse: ",i4, g14.6)') iter, rmse
        if (debug_write) write (debug_unit, '(10f9.2)') y_int             
        if (debug_write) write (debug_unit, '(a)') "ym_int"
        if (debug_write) write (debug_unit, '(12f9.4)') ym_int
        if (debug_write) write (debug_unit, '(a)') "ym"
        if (debug_write) write (debug_unit, '(12f9.4)') ym
        if (debug_write) write (debug_unit, '(a)') "y_int_out1"
        if (debug_write) write (debug_unit, '(12f9.4)') y_int_out1(iter, :)
        
    end do
   
end subroutine rm_int
    
end module mp_interp_rymes_and_meyer_subs