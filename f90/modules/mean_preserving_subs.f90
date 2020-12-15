module mean_preserving_subs
! subroutines for mean-preserving interpolation
    
    implicit none

    integer             :: debug_unit = 10
    integer             :: out_unit = 6
    logical             :: debug_em = .false.

contains
    
subroutine enforce_mean(nctrl, ntargs, nsubint, tol, xm_targ, xd, xmiss)
! adjustes daily values to match target mean

    implicit none
    
    integer(4), parameter   :: maxtargs = 7305
    integer(4), intent(in)  :: nctrl
    integer(4), intent(in)  :: ntargs
    integer(4), intent(in)  :: nsubint(nctrl)
    real(8), intent(in)     :: tol
    real(8), intent(in)     :: xm_targ(nctrl)
    real(8), intent(inout)  :: xd(ntargs)
    real(8), intent(in)     :: xmiss
    !
    integer(4), parameter   :: maxwgts = 31, maxiter=100
    real(8)                 :: xdm(nctrl)
    real(8)                 :: nonzero_sum(nctrl), xd_adjust_sum(nctrl), xd_next_sum(nctrl)
    real(8)                 :: xd_next_mean(nctrl), new_diff(nctrl)
    real(8)                 :: diff(nctrl), xd_adjust 
    
    integer(4)              :: ii, ib, ie, j, jj, k, niter(nctrl), ntotiter
    
    if (debug_em) write (debug_unit,'(a)') "In enforce_mean"
    if (debug_em) write (debug_unit,'("tol: ", f8.5)') tol
    if (debug_em) write (debug_unit,'("nctrl, ntargs: ", 2i8/)') nctrl, ntargs
        
    ! get means of current daily values
    
    call interval_mean(nctrl, nsubint, ntargs, xd, xmiss, xdm)
    
    ii = 0
    ntotiter = 1
    xd_adjust_sum = 0.0d0;  xd_next_sum = 0.0d0
    do k = 1, nctrl
        niter(k) = 0
        if (debug_em) write (debug_unit, '(/"Interval ", i8)') k 
        ib = ii + 1; ie = ib + nsubint(k) - 1
        
        if (debug_em) write (debug_unit, '("xm_targ,xdm,xm_targ-xdm: ",3f12.6)') xm_targ(k), xdm(k), xm_targ(k) - xdm(k)
        
        ! if xm(k) = 0.0, set all days to 0.0
        if (xm_targ(k) .eq. 0.0d0) then 
            if (debug_em) write(debug_unit, '(a)') "       k      ib      ie   ndays xm_targ"
            if (debug_em) write(debug_unit, '(4i8, f12.6)')  k, ib, ie, ie - ib + 1, xm_targ(k)
            if (debug_em) write (debug_unit, '(a)') "Setting all subinterval values to 0.0"
            xd(ib:ie) = 0.0d0
            if (debug_em) write (debug_unit, '(10g14.6)') xd(ib:ie) 
            
        else if (xm_targ(k) .eq. xmiss) then
            if (debug_em) write(debug_unit, '(a)') "       k      ib      ie   ndays xm_targ"
            if (debug_em) write(debug_unit, '(4i8, f12.6)')  k, ib, ie, ie - ib + 1, xm_targ(k)
            if (debug_em) write (debug_unit, '(a)') "Setting all subinterval values to missing"
            xd(ib:ie) = xmiss
            if (debug_em) write (debug_unit, '(10g14.6)') xd(ib:ie) 
        else    
            
            ! check for differences in actual mean and mean of interpolated values
            diff(k) = xm_targ(k) - xdm(k)
            if (debug_em) write(debug_unit, '(a)') "       k      ib      ie   ndays xm_targ"// &
                "         xdm            diff         tol  dabs(diff)"
            if (debug_em) write(debug_unit, '( 4i8, 5f12.6)')  & 
                k, ib, ie, ie - ib + 1, xm_targ(k), xdm(k), diff(k), tol, dabs(diff(k))
            
            ! if absolute value of difference is greater than tol, adjust subinterval values            
            if (dabs(diff(k)) .gt. tol) then
                
                if (debug_em) write (debug_unit, '(a)') "Adjusting subinterval values..."
                
                ! get sum of nonzero values
                if (debug_em) write (debug_unit, '(a)') "       j      jj          xd nonzero_sum"
                nonzero_sum = 0.0d0     
                jj = 0
                do j = ib, ie
                    jj = jj + 1
                    if (xd(j) .ne. 0.0d0)  nonzero_sum(k) = nonzero_sum(k) + xd(j)   
                    if (debug_em) write (debug_unit, '(2i8, 2f12.6)') j, jj, xd(j), nonzero_sum(k)
                end do 
                    
                ! apply porportions of differences
                if (debug_em) write (debug_unit, '(a)') "       j      jj          xd   xd_adjust"
                jj = 0
                do j = ib, ie
                    jj = jj + 1
                    if (xd(j) .ne. 0.0d0) then 
                        xd_adjust = ((xd(j) / nonzero_sum(k)) * diff(k)) * nsubint(k) 
                        xd(j) = xd(j) + xd_adjust
                        xd_next_sum(k) = xd_next_sum(k) + xd(j)
                        xd_adjust_sum(k) = xd_adjust_sum(k) + xd_adjust
                        if (debug_em) write (debug_unit, '(2i8, 3f12.6)') j, jj, xd(j), xd_adjust
                    else
                        if (debug_em) write (debug_unit, '(2i8, 3f12.6)') j, jj, xd(j)
                    end if   
                end do
                xd_next_mean(k) = xd_next_sum(k) / nsubint(k)
                new_diff(k) = xm_targ(k) - xd_next_mean(k)
                if (debug_em) write (debug_unit, '(a)') "xd_adjust_sum xd_next_sum xd_next_mean  new_diff"
                if (debug_em) write (debug_unit, '(4f12.6)') xd_adjust_sum(k), xd_next_sum(k), xd_next_mean(k), new_diff(k)
            else  
            
                if (debug_em) write (debug_unit, '(a)') "Adopting subinterval values..."
            
            end if
            
        end if
        ii = ii + nsubint(k)

    end do
    
end subroutine enforce_mean

subroutine dayinterp(nm,nd,monlen,zm,zd)
! Interpolate pseudo-daily values of monthly data.  Not mean-preserving.

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: zm(nm)
    real(8), intent(out)    :: zd(nd)
    
    integer(4)              :: nm1
    !integer(4)              :: midmon(nm),midmon2(0:nm+1)
    real(8)                 :: midmon(nm),midmon2(0:nm+1)
    real(8)                 :: zm2(0:nm+1)
    integer                 :: i,m
    
    !call midmonth_int(nm,monlen,midmon)
    call midmonth_real(nm,dble(monlen),midmon)
    
    ! pad data at beginning (m=0) and end (m=13)
    nm1=nm+1
    zm2(1:nm)=zm
    zm2(0)=zm(nm)
    zm2(nm1)=zm(1)
    midmon2(1:nm)=midmon
    midmon2(0)=1-(nd-midmon(nm))-2
    midmon2(nm1)=nd+midmon(1)   
    
    ! linear pseudo-daily interpolation
    do i=1,nd
        ! find month day i lies in
        do m=1,nm+1
            if (i.gt.midmon2(m-1) .and. i.le.midmon2(m)) exit
        end do   
        zd(i)=(dble(i-midmon2(m-1))/dble(midmon2(m)-midmon2(m-1)))*(zm2(m)-zm2(m-1))+zm2(m-1)       
    end do
    
end subroutine dayinterp

subroutine dayspread(nm,nd,monlen,zm,zd)
! block fill daily values from monthly means

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: zm(nm)
    real(8), intent(out)    :: zd(nd)
    
    integer                 :: i,j,m
    
    i=0
    do m=1,nm
        do j=1,monlen(m)
            i=i+1
            zd(i)=zm(m)
        end do
    end do
    
end subroutine dayspread

subroutine midmonth_int(nm,monlen,midmon)
! gets mid-month day number 

    implicit none
    
    integer(4), intent(in)  :: nm
    integer(4), intent(in)  :: monlen(nm)
    integer(4), intent(out) :: midmon(nm)
    
    integer(4)              :: m,endday(nm)

    ! midmonth day numbers
    m=1
    midmon(m)=ceiling(dble(monlen(m))/2.0d0)
    endday(m)=monlen(m)
    do m=2,nm
        midmon(m)=ceiling(dble(monlen(m))/2.0d0)+endday(m-1)
        endday(m)=endday(m-1)+monlen(m)
    end do
    !write (*,'(a)') "midmon, endday"
    !write (*,'(12i9)') midmon
    !write (*,'(12i9)') endday
    
end subroutine midmonth_int

subroutine midmonth_real(nm,monlen,midmon)
! gets mid-month day number 

    implicit none
    
    integer(4), intent(in)  :: nm
    real(8), intent(in)     :: monlen(nm)
    real(8), intent(out)    :: midmon(nm)
    
    real(8)              :: endday(nm)
    integer(4)           :: m

    ! midmonth day numbers
    m=1
    midmon(m)=monlen(m)/2.0d0
    endday(m)=monlen(m)
    do m=2,nm
        midmon(m)=(monlen(m)/2.0d0)+endday(m-1)
        endday(m)=endday(m-1)+monlen(m)
    end do
    !write (*,'(a)') "midmon, endday"
    !write (*,'(12f9.2)') midmon
    !write (*,'(12f9.2)') endday
    
end subroutine midmonth_real

subroutine dzero(nm,nd,monlen,xm,xd0)
! enforces 0.0 values of interpolated daily data when the monthly mean is 0.0

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: xm(nm)
    real(8), intent(inout)  :: xd0(nd)
    
    integer(4)              :: i,m,j

    ! zero daily values in months were xm=0.0
    i=0
    do m=1,nm
        do j=1,monlen(m)
            i=i+1
            if (xm(m).eq.0.0) xd0(i)=0.0          
        end do
    end do
    
    ! zero all other daily values
    do i=1,nd
        if (xd0(i).le.0.0) xd0(i)=0.0
    end do

end subroutine dzero

subroutine interval_mean(ninterval, nsubint, n, v_int, vmiss, vm_int)
    
    implicit none
        
    integer(4), intent(in)          :: ninterval, n
    integer(4), intent(in)          :: nsubint(ninterval)                  ! interval lengths
    real(8), intent(in)             :: v_int(n)                    ! input values
    real(8), intent(in)             :: vmiss                            ! missing/fill value
    real(8), intent(out)            :: vm_int(ninterval)                    ! output means
    
    integer(4)                      :: i, j, k, npts
    
    vm_int = 0.0d0
    i = 0
    do k = 1, ninterval
        npts = 0
        do j = 1, nsubint(k)
            i = i + 1 
            if (v_int(i) .ne. vmiss) then
                vm_int(k) = vm_int(k) + v_int(i)
                npts = npts + 1 
            end if
        end do

        if (npts .gt. 0) then     
            vm_int(k) = vm_int(k) / dble(npts)
        else 
            vm_int(k) = vmiss
        end if
    end do
    
end subroutine interval_mean

subroutine interp_stat(nctrl, ym, ym_int, rmse)

    implicit none
    
    integer(4), intent(in)      :: nctrl
    real(8), intent(in)         :: ym(nctrl), ym_int(nctrl)
    real(8), intent(out)        :: rmse
    
    integer(4)                  :: k
    
    rmse = 0.0d0

    do k = 1, nctrl
        rmse = rmse + (ym(k) - ym_int(k))*(ym(k) - ym_int(k))
    end do
    rmse = sqrt(rmse / (nctrl - 1))
    
end subroutine interp_stat

subroutine step_plot(nctrl, ntarg, nsubint, zm, zx, zint)
! block fill daily values from monthly means

    implicit none
    
    integer(4), intent(in)  :: nctrl, ntarg
    integer(4), intent(in)  :: nsubint(nctrl)
    real(8), intent(in)     :: zm(nctrl)
    real(8), intent(out)    :: zx(ntarg)   ! adjusted abcissas to make nice step plot
    real(8), intent(out)    :: zint(ntarg)
    
    real(8)                 :: halfstep
    
    integer                 :: i, j, k
    
    i = 0
    do k = 1, nctrl
        do j = 1, nsubint(k)
            i = i + 1
            zx(i) = dble(i)
            zint(i) = zm(k)
        end do
    end do
    
    ! adjust abcissas
    
    halfstep = (zx(2) - zx(1)) / 2.0d0
    i = 0
    do k = 1, nctrl - 1
        i = i + nsubint(k)
        zx(i) = zx(i) + halfstep
        zx(i + 1) = zx(i) 
    end do
   
end subroutine step_plot

end module mean_preserving_subs