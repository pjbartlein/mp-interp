module mp_interp_rymes_and_meyers_subs
! subroutines for implementing the Rymes abd Meyer (2001) "iterative smoothing" mean-preserving interpolation
    
! Rymes, M.D. and D.R. Meyers (2001) Mean preserving alogrithm for smoothing interpolating averaged data.
! Solar Energy 71:225-231.
    
implicit none

    integer             :: debug_unit=10
    logical             :: debug_write = .false.  
    
contains

subroutine mp_interp_rymes_and_meyers(n_outer, n_inner, nctrl, ym, yfill, x_ctrl, nsubint,  &
    lowerbound, lower, upperbound, upper, npad, no_negatives, match_mean, tol, & 
    ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

    use pseudo_daily_interp_subs

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
    
    real(8), intent(inout)        :: y_int(ntargs)                    ! interpolated values
    real(8), intent(inout)        :: ym_int(nctrl)                    ! mean of interpolated values
                                                                    ! inner interval = e.g. months, outer interval = e.g. year
    
    integer(4), parameter       :: maxit = 50                       ! max number of iterations
    integer(4), parameter       :: init_type = 2                    ! initialization type, 1 = spread, 2 = linear interp
    
    real(8)             :: ym_in(max_nctrl_in), x_ctrl_in(max_nctrl_in)
    real(8)             :: x_targ_in(max_ntargs_in)
    real(8)             :: y_int_out(max_ntargs_in), ym_int_out(max_ntargs_in)
    real(8)             :: lower_in(max_nctrl_in), upper_in(max_nctrl_in)
    real(8)             :: rmse
    
    integer(4)          :: nsubint_in(max_nctrl_in)
    
    integer(4)          :: n, nn
    integer(4)          :: beg_inner, end_inner, beg_ctrl, end_ctrl, nctrl_in, beg_targ, end_targ, ntarg_in, beg_int, end_int, nint_out
    
    if (debug_write) write (debug_unit,'(a)') "In mp_interp_rymes-and-meyer()"

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
        lower_in(1:nctrl_in) = lower(beg_ctrl:end_ctrl)
        upper_in(1:nctrl_in) = upper(beg_ctrl:end_ctrl)
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
    
            call rm_int(nctrl_in, ym_in, yfill, x_ctrl_in, nsubint_in, lowerbound, lower_in, upperbound, upper_in, &
                init_type, ntarg_in, maxit, x_targ_in, y_int_out, ym_int_out)
            
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
    
end subroutine mp_interp_rymes_and_meyers
    
    
subroutine rm_int(nctrl, ym, ymiss, x_ctrl, nsubint, lowerbound, lower, upperbound, upper, &
    init_type, ntargs, maxit, x_targ, y_int, ym_int) 
! pseudo-daily interpolated values
! From Rymes, M.D. and D.R. Myers, 2001. Solar Energy 71:225-231
    
    use pseudo_daily_interp_subs

    implicit none

    integer(4), intent(in)              :: nctrl, ntargs
    real(8), intent(in)                 :: ym(nctrl)                    ! input AVG[]
    real(8), intent(in)                 :: ymiss                        ! missing / fill value
    real(8), intent(in)                 :: x_ctrl(nctrl)                ! x_cntrl (used for debugging)
    integer(4), intent(in)              :: nsubint(nctrl)               ! number of subintervals
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
    
    integer(4)          :: i, j, k, iter ! i indexes targets, k indexes controls
    
    if (debug_write) write (debug_unit, '(a)') "In rm_int"
    if (debug_write) write (debug_unit, *) nctrl, ntargs, maxit
    if (debug_write) write (debug_unit, '(12f12.6)') ym
    if (debug_write) write (debug_unit, '(12i12)') nsubint
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
    ybeg = y_int(ntargs)
    yend = y_int(1)
    
    if (debug_write) write (debug_unit,'(a)') "initialize:  y_int"
    if (debug_write) write (debug_unit,'(10f12.6)') y_int
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
                if (debug_write) write (debug_unit, '("k,i,ck(k),ym(k),y_int(i),ym(k)-y_int(i): ", 2i6, 4f12.6)') &
                    k, i, ck(k), ym(k), y_int(i), ym(k)-y_int(i)
            end do
            
            ck(k) = ck(k) / dble(ntargs) 
            
            do i = ibeg(k), iend(k)
                y_int(i) = y_int(i) + ck(k)                        ! apply corrections
            end do

            if (debug_write) write (debug_unit,'("iter,k,ibeg,iend, ck: ", 4i4, f12.6)') iter, k, ibeg(k), iend(k), ck(k)
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
                if (debug_write) write (debug_unit,'("k,lflag,ym,ym_int,diff,fk: ", 2i5,4f12.6)') &
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
                if (debug_write) write (debug_unit,'("k,uflag,ym,ym_int,diff,fk: ", 2i5,4f12.6)') &
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
        
        if (debug_write) write (debug_unit, '("iteration, rmse: ",i4, g14.6)') iter, rmse
        if (debug_write) write (debug_unit, '(10f9.2)') y_int             
        if (debug_write) write (debug_unit, '(a)') "ym_int"
        if (debug_write) write (debug_unit, '(12f12.6)') ym_int
        if (debug_write) write (debug_unit, '(a)') "ym"
        if (debug_write) write (debug_unit, '(12f12.6)') ym
        if (debug_write) write (debug_unit, '(a)') "y_int_out1"
        if (debug_write) write (debug_unit, '(12f12.6)') y_int_out1(iter, :)
        
    end do
   
end subroutine rm_int
    
end module mp_interp_rymes_and_meyers_subs
