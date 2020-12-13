program demo_1yr_rymes_and_meyer_01
! demonstrates Rhymes and Meyer (2001) iterative smoothing pseudo-daily interpolation  -- single year of different data sets
    
use mean_preserving_subs
use mp_interp_rymes_and_meyer_subs
    
implicit none

! The data are structured as follows:
! there are n_outer outer intervals (e.g. years), n_inner "inner" intervals (e.g. months), and 
! nsubint subintervals withing each inner intervals (e.g. days), for transparency in this program, ny, nm, and nd
! are used to represent these intervals and subinterfals.
! In this example, a single year's worth of data is used.

integer(4), parameter  :: ny = 1, nm = 12           ! number of years and months in a year
integer(4), parameter  :: nctrl = ny * nm           ! total number of control points
integer(4), parameter  :: maxtarg = 365             ! total number of days
integer(4), parameter  :: max_nctrl_in = 12         ! maximum number of inner intervals (months) in an outer interval (years)
integer(4), parameter  :: max_ntargs_in = 366       ! maximum number of subintervals (days) in an outer interval (years)
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

! control variables
logical             :: no_negatives = .false.       ! no negative interpolated values (e.g. precip)
logical             :: match_mean = .false.         ! use enfore_mean() to match observed means
logical             :: lowerbound = .false., upperbound = .false.  ! set lower and upper bounds?
real(8)             :: lower(nctrl), upper(nctrl)   ! upper and lower bounds

integer(4)          :: npad                         ! number of months to pad at each end of an inner interval

integer(4)          :: demo_ctrl_beg(nm), demo_ctrl_mid(nm), demo_nsubint(nm)
real(4)             :: demo_tas(nm), demo_pr(nm)

data demo_ctrl_beg /     1,    32,    60,    91,   121,   152,   182,   213,   244,   274,   305,   335 /
data demo_ctrl_mid /    15,    43,    74,   104,   135,   165,   196,   227,   257,   288,   318,   349 /
data demo_nsubint  /    31,    28,    31,    30,    31,    30,    31,    31,    30,    31,    30,    31 /
data demo_tas      / -6.77, -4.33,  1.52,  8.35, 14.74, 20.29, 23.22, 22.17, 16.94, 10.35,  3.08,  4.29 /
data demo_pr       /  20.0,   0.0,  30.0,   0.0,   0.5,  10.0,   0.0,   0.0,   5.0,   5.0,   5.0,   0.0 /

! all examples
tol = 0.0001
x_ctrl = dble(demo_ctrl_mid)
ntargs = maxtarg
do i = 1, ntargs
    x_targ(i) = ctrl_beg(1) + dble(i-1)
end do
nsubint = demo_nsubint

npad = 2

! example 1
write (*,'(a)') "Example 1:  temperature, no adjustment of means"
ym = demo_tas
no_negatives = .false.; match_mean = .false. 
lowerbound = .false.
upperbound = .false.
lower = 0.0d0 / 0.0d0
upper = 0.0d0 / 0.0d0

call mp_interp_rymes_and_meyer(ny, nm, nctrl, ym, yfill, x_ctrl, nsubint, & 
    lowerbound, lower, upperbound, upper, npad, no_negatives, match_mean, tol, & 
    ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

write (*, '(a)') " "
write (*, '(" input: ", 12f8.2)') ym
write (*, '(a)') " y_int: "
write (*, '(15f8.2)') y_int
write (*, '(a)') " "
write (*, '("     ym: ", 12f9.3)') ym
write (*, '(" ym_int: ", 12f9.3)') ym_int
write (*, '("   diff: ", 12f9.3)') ym - ym_int
write (*, '(a)') " "
write (*, '(a)') "There are amall differences beteen the means of the input data and the monthly means of the"
write (*, '(a)') "interpolated data.  Next, enforce the equality of the mean values by adjusting the interpolated ones."

! example 2
write (*, '(a)') " "
write (*,'(a)') "Example 2:  temperature, enforce equality of input and interpolateed mean values"
ym = demo_tas
no_negatives = .false.; match_mean = .true. 
lowerbound = .false.
upperbound = .false.
lower = 0.0d0 / 0.0d0
upper = 0.0d0 / 0.0d0

call mp_interp_rymes_and_meyer(ny, nm, nctrl, ym, yfill, x_ctrl, nsubint, & 
    lowerbound, lower, upperbound, upper, npad, no_negatives, match_mean, tol, & 
    ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

write (*, '(a)') " "
write (*, '(" input: ", 12f8.2)') ym
write (*, '(a)') " y_int: "
write (*, '(15f8.2)') y_int
write (*, '(a)') " "
write (*, '("     ym: ", 12f9.3)') ym
write (*, '(" ym_int: ", 12f9.3)') ym_int
write (*, '("   diff: ", 12f9.3)') ym - ym_int
write (*, '(a)') " "
write (*, '(a)') "The adjustment of the mean values with tol = 0.0001 reduces the difference."
write (*, '(a)') "Note that there is little practical difference between the interpolated values above."

! example 3
write (*, '(a)') " "
write (*,'(a)') "Example 3:  precipitation rate, no_negetives and lowerbound = .true., but no adjustment of means"
ym = demo_pr
no_negatives = .true.; match_mean = .false. 
lowerbound = .true.
upperbound = .false.
lower = 0.0d0 
upper = 0.0d0 / 0.0d0

call mp_interp_rymes_and_meyer(ny, nm, nctrl, ym, yfill, x_ctrl, nsubint, & 
    lowerbound, lower, upperbound, upper, npad, no_negatives, match_mean, tol, & 
    ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

write (*, '(a)') " "
write (*, '(" input: ", 12f8.2)') ym
write (*, '(a)') " "
write (*, '(a)') " y_int: "
write (*, '(15f8.2)') y_int
write (*, '(a)') " "
write (*, '("     ym: ", 12f9.3)') ym
write (*, '(" ym_int: ", 12f9.3)') ym_int
write (*, '("   diff: ", 12f9.3)') ym - ym_int
write (*, '(a)') " "
write (*, '(a)') "There are some relatively large differences beteen the means of the input data and the monthly means of the"
write (*, '(a)') "interpolated data.  Next, enforce the equality of the mean values by adjusting the interpolated ones."

! example 4
write (*, '(a)') " "
write (*,'(a)') "Example 4:  precipitation rate, no_negetivesl and lower bound = .true., enforce equality of input and interpolateed mean values"
ym = demo_pr
no_negatives = .true.; match_mean = .true. 
lowerbound = .true.
upperbound = .false.
lower = 0.0d0 
upper = 0.0d0 / 0.0d0

call mp_interp_rymes_and_meyer(ny, nm, nctrl, ym, yfill, x_ctrl, nsubint, & 
    lowerbound, lower, upperbound, upper, npad, no_negatives, match_mean, tol, & 
    ntargs, x_targ, max_nctrl_in, max_ntargs_in, y_int, ym_int)

write (*, '(a)') " "
write (*, '(" input: ", 12f8.2)') ym
write (*, '(a)') " y_int: "
write (*, '(15f8.2)') y_int
write (*, '(a)') " "
write (*, '("     ym: ", 12f9.3)') ym
write (*, '(" ym_int: ", 12f9.3)') ym_int
write (*, '("   diff: ", 12f9.3)') ym - ym_int
write (*, '(a)') " "
write (*, '(a)') "The adjustment of the mean values with tol = 0.0001 reduces the difference."
write (*, '(a)') " "
end