program test_epstein
! tests Epstein (1991) harmonic pseudo-daily interpolation
    
use pseudo_daily_interp_subs
use mp_interp_epstein_subs
    
implicit none

! The data are structured as follows:
! there are n_outer outer intervals (e.g. years), n_inner "inner" intervals (e.g. months), and 
! nsubint subintervals withing each inner intervals (e.g. days), for transparency in this program, ny, nm, and nd
! are used to represent these intervals and subinterfals
! In this example, 2 m air temperature (tas) and precipitation-rate data (pr) for a representative Mediterranean- and 
! monsoon-climate grid point from NCEP-DOEv2 Reanalysis II data set are used to illusatrate the approach and graphically 
! compare the results to "observed" (i.e. reanalysis) daily data.

integer(4), parameter  :: nvars = 1                 ! number of demonstration time series
integer(4), parameter  :: ny = 1, nm = 12           ! number of years and months in a year
integer(4), parameter  :: nctrl = ny * nm           ! total number of control points
integer(4), parameter  :: maxtarg = 365             ! total number of days
integer(4), parameter  :: max_nctrl_in = 12         ! maximum number of inner intervals (months) in an outer interval (years)
integer(4), parameter  :: max_ntargs_in = 365       ! maximum number of subintervals (days) in an outer interval (years)
real(8)             :: ym(nctrl)                    ! control data (e.g. observed monthly means)            
real(8)             :: x_ctrl(nctrl)                ! time values of control points (e.g. mid-month day number in the range 1 to maxtarg)
real(8)             :: x_targ(maxtarg)              ! time values of target points (e.g. day number in the range 1 to maxtarg
real(8)             :: yrmn(nctrl), yrmndy(maxtarg) ! decimal year values, e.g. of months and days
integer(4)          :: iy_mon(nctrl), im_mon(nctrl) ! integer year and month of input  data (e.g. observed monthly means)
integer(4)          :: ctrl_beg(nctrl)              ! time values (e.g. days) of the beginning of each innter interval (e.g. month)
integer(4)          :: ctrl_mid(nctrl)              ! time values (e.g. days) of the middle of each innter interval (e.g. month)
integer(4)          :: nsubint(nctrl)               ! number of subintervals (e.g. days) in each inner interval

real(8)             :: y_int(maxtarg)               ! interpolated value for each subinterval (e.g. days)  
real(8)             :: ym_int(nctrl)                ! inner interval (e.g. monthsly) means of the subinterval values (e.g. days)

real(8)             :: y_obs(maxtarg)               ! observed subinterval (e.g. days) values for comparison
integer(4)          :: iy(maxtarg), im(maxtarg), id(maxtarg)    ! year, month and day indices of observed daily values
real(8)             :: step_x(maxtarg), step_y(maxtarg)         ! x- and y-values for constructing a step plot of observed monthly values

real(8)             :: tol                          ! tolerance value for enforce_mean()
real(8)             :: yfill = 1e32                 ! missing/fill value

integer(4)          :: i, m, n, mm, nn, ivar        ! various indices
integer(4)          :: ntargs                       ! total number of target points

character(1)        :: header
character(16)       :: methodname = "epstein_"
character(1024)     :: sourcepath, interppath, dataname, infile

real(4)             :: total_secs, loop_secs

! control variables
logical             :: no_negatives = .false.       ! no negative interpolated values (e.g. precip)
logical             :: match_mean = .false.         ! use enfore_mean() to match observed means
logical             :: smooth = .true.              ! smooth across outer intervals

! source and output pathss, modify as necessary
sourcepath = "e:\Projects\MeanPreserving\data\test_methods\"
interppath = "e:\Projects\MeanPreserving\data\test_methods\"

total_secs = secnds(0.0)

do ivar = 1, nvars
    
    loop_secs = secnds(0.0)
    
    select case(ivar)
    case (1)
        dataname = "epstein_fig01"
        no_negatives = .false.
        match_mean = .true.
        tol = 0.01
        smooth = .false.
    case default
        stop "ivar"
    end select
        
    infile = trim(dataname)//".csv"
    open (1, file = trim(sourcepath)//trim(infile))
    read (1,'(a)') header
    
    open (2, file = trim(sourcepath)//"YrMnDy.csv")
    read (2,'(a)') header
    
    open (3, file = trim(interppath)//trim(dataname)//"_mon.csv")
    open (4, file = trim(interppath)//trim(dataname)//"_day.csv")

    ! read monthly data to interpolate, and observed daily data for comparison
    nn = 0
    do n = 1, ny
        do m = 1, nm
            nn = nn + 1
            read (1,*) iy_mon(nn), im_mon(nn), yrmn(nn), ctrl_beg(nn), ctrl_mid(nn), nsubint(nn), ym(nn)
!            write (10,'(3i5, f14.8, 2i8, i3, f12.5)') nn, iy_mon(nn), im_mon(nn),  &
!                yrmn(nn), ctrl_beg(nn), ctrl_mid(nn), nsubint(nn), ym(nn)
        end do
    end do
    close (1)
    
    do nn = 1, maxtarg
        read (2,*) iy(nn), im(nn), id(nn), yrmndy(nn)
!        write (10,'(4i5, f14.8, 3f12.5)') nn,iy(nn), im(nn), id(nn), yrmndy(nn)
    end do
    ntargs = maxtarg
    close(2)

    x_ctrl = dble(ctrl_mid)
    !write (10,'(12f8.2)') x_ctrl

    do i = 1, ntargs
        x_targ(i) = ctrl_beg(1) + dble(i-1)
    end do
    !write (10,'(10f12.5)') x_targ
    
    ! mean-preserving interpolation

    call mp_interp_epstein(ny, nm, nctrl, ym, yfill, x_ctrl, nsubint, &
        smooth, no_negatives, match_mean, tol, ntargs, max_nctrl_in, max_ntargs_in, y_int, ym_int)

    ! write out monthly time series
    write (3, '(a)') "Year, Month, YrMn, MonBeg, MidMonth, ndays, ym, ym_int"
    nn = 0
    do n = 1, ny
        do m = 1, nm
            nn = nn + 1
            write (3, '(2(i4, ", "), f13.8, ", ", 2(i9, ", "), i3, 2(", ", g14.6))') & 
                iy_mon(nn), im_mon(nn), yrmn(nn), ctrl_beg(nn), ctrl_mid(nn), nsubint(nn), ym(nn), ym_int(nn)
            !write (10,'(3i5, f14.8, 2i8, i3, f12.5)') & 
            !    nn, iy(nn), im(nn), yrmn(nn), ctrl_beg(nn), ctrl_ ym_int(nn)
        end do
    end do
    
    ! write out daily time series, including data to make a "step plot" of 
    write (4, '(a)') "Year, Month, Day, YrMnDy, y_obs, y_int, step_x, step_y"
    nn = 0; mm = 0
    do n = 1, ny   
        do m = 1, nm
            mm = mm + 1
            do i = 1, nsubint(mm)
                nn = nn + 1       
                write (4, '(3(i4, ", "), f13.8, 2(", ", g14.6), ", ",f13.8, ", ", g14.6)') & 
                    iy(nn), im(nn), id(nn), yrmndy(nn), y_obs(nn), y_int(nn), step_x(nn), step_y(nn)
                !write (10,'(5i5, f14.8, 4f12.5)') n,nn, iy(nn), im(nn), id(nn), yrmndy(nn), & 
                !    y_obs(nn), y_int(nn), step_x(nn), step_y(nn)
            end do
        end do
    end do

    close(1); close(2); close(3)



    close(1); close(3)

    write (*,'("loop, seconds: ", i4, f8.4)') ivar, secnds(loop_secs)

end do

write (*,'("total seconds:     ", f8.4)') secnds(total_secs)

end
