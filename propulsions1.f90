program orbits_n_thrust
use dislin
implicit none

! declare variables
integer :: h, i, j, k, hr, minute, b, c, d, p, q, z, l
integer, parameter :: wp = selected_real_kind(p=12)
real(wp) :: F, F1, F2, Isp, Isp1, Isp2, dt, temp1, temp2, temp3, tol, time, s, dv1, dv2, dv, m1, m2, m
real(wp), parameter :: pi = 3.141592654_wp, u = 3.9860044e5_wp, g = 9.806_wp
real(wp), parameter :: a_leo = 8530._wp, a_meo = 13200._wp
real(wp) :: r_perigee, r_apogee
real(wp), dimension(2) :: e, a
real(wp), dimension(5) :: xo, xn
real(wp), allocatable, dimension(:,:) :: burn1, coast, burn2
character :: method, cont

cont = 'y'

do while (cont == 'y' .or. cont == 'Y')

! initialize variables
xn(3) = a_leo
xn(4) = 0._wp
xn(1) = 0._wp
xn(2) = sqrt(u / a_leo)
xo = 0._wp
time = 0._wp

call prompt_user('Enter value for initial mass (kg): ',xn(5),35)
m = xn(5)
call prompt_user('Enter value for thrust (N) for first motor: ',F1,44)
F = F1
call prompt_user('Enter value for Isp (s) for first motor: ',Isp1,41)
Isp = Isp1
call prompt_user('Enter value for thrust (N) for second motor: ',F2,45)
call prompt_user('Enter value for Isp (s) for second motor: ',Isp2,42)
call prompt_user('Enter value for time step (s): ',dt,31)
temp1 = dt
call prompt_user('Enter tolerance value for final orbit radius (km): ',tol,51)
method = 'a'
do while(method/='t'.and.method/='r'.and.method/='T'.and.method/= 'R')
	write(*,*)
	write(*,*) 'Use Trapezoidal Rule or Runge-Kutta (t / r): '
	read(*,*) method
end do
write(*,*)
write(*,*) 'Enter graph update frequency (sec): '
read(*,*) temp2
b = floor(temp2 / dt)

! initialize graphics
call metafl('XWIN')
call window(640,1,1280,1000)
call page(5000,4000)
call disini
call axslen(2000,2000)
call axspos(500,2500)
call labtyp('VERT','X')
call graf(-15000., 15000., -15000., 3000., -15000., 15000., -15000., 3000.)

! perform calculations
j = 0
r_apogee = a_leo

open(unit = 10, status = 'scratch')

! First burn
do while (r_apogee < a_meo - tol .or. r_apogee > a_meo + tol)
	
	do i = 1, 5
		xo(i) = xn(i)
	end do
	
	if (method == 't' .or. method == 'T') then
		call trapezoidal_rule(xo, xn, dt, F, Isp, u)
	else
		call runge_kutta(xo, xn, dt, F, Isp, u)
	end if
	
	call instant_orbit(xn(3), xn(1), xn(2), u, a(2), e(2), r_perigee, r_apogee)
	if (r_apogee > a_meo + tol) then
		do i = 1, 5
			xn(i) = xo(i)
		end do
		dt = dt / 2._wp
	else
		time = time + dt
		j = j + 1
		write(10,*) xn(1), xn(2), xn(3), xn(4), xn(5)
	end if
	
	if (xn(5) <= 1000._wp) then
		write(*,*)
		write(*,*) ' * * * *   You ran out of fuel!   * * * * '
		exit
	end if
end do

allocate(burn1(j,5))
rewind(10)
do i = 1, j
	read(10,*) burn1(i,1), burn1(i,2), burn1(i,3), burn1(i,4), burn1(i,5)
end do
close(10)
m1 = m - burn1(j,5)
call delta_v(m, burn1(j,5), Isp, g, dv1)
dv2 = sqrt(u/a_meo) - sqrt(u*(2._wp / r_apogee - 1._wp / a(2)))
m2 = 1000._wp * (2.7182818_wp ** (dv2 / g / Isp2 * 1000._wp) - 1._wp)
temp3 = m2 * g * Isp2 / F2
open(unit = 20, status = 'scratch')
k = 0
dt = temp1
l = floor(temp3 / dt / 2._wp)

p = 500
q = 2700
z = 4100
call messag('Initial mass: ',q,p)
call number(real(m),5,z,p)
p = p + 100
p = p + 100
call messag('* * * *   First Burn   * * * *',q,p)
p = p + 100
call messag('Mass used during first burn (kg):',q,p)
call number(real(m1),5,z,p)
p = p + 100
call messag('Mass after first burn (kg):',q,p)
call number(real(burn1(j,5)),5,z,p)
p = p + 100
call messag('Delta V form first burn (km/s):',q,p)
call number(real(dv1),5,z,p)
p = p + 100
call messag('Predicted apogee after first burn (km):',q,p)
call number(real(r_apogee),4,z,p)
p = p + 100
call messag('Obrit eccentricity right after burn:',q,p)
call number(real(e(2)),7,z,p)
p = p + 100
call messag('Duration of first burn:',q,p)
hr = floor(time / 3600._wp)
call number(real(hr),-1,z,p)
call messag(':',999,999)
minute = floor(mod(time,3600._wp) / 60._wp)
call number(real(minute),-1,999,999)
call messag(':',999,999)
call number(real(mod(time,60._wp)),2,999,999)
call sendbf

temp2 = time
time = 0._wp
F = 0._wp

! After burn coast into meo orbit
no_thrust: do
	
	do i = 1, 5
		xo(i) = xn(i)
	end do
	
	if (method == 't' .or. method == 'T') then
		call trapezoidal_rule(xo, xn, dt, F, Isp, u)
	else
		call runge_kutta(xo, xn, dt, F, Isp, u)
	end if
	
	if (xn(3) < xo(3)) then
		if (dt <= temp1 * .125_wp) then
			exit no_thrust
		else
			do i = 1, 5
				xn(i) = xo(i)
			end do
			dt = dt / 2._wp
		end if
	else
		time = time + dt
		k = k + 1
		write(20,*) xn(1), xn(2), xn(3), xn(4), xn(5)
	end if
	
	call instant_orbit(xn(3), xn(1), xn(2), u, a(1), e(1), r_perigee, r_apogee)
	
end do no_thrust


allocate(coast(k-l,5))
rewind(20)
do i = 1, k-l
	read(20,*) coast(i,1), coast(i,2), coast(i,3), coast(i,4), coast(i,5)
end do
close(20)
dt = temp1
temp3 = time

p = p + 100
p = p + 100
call messag('* * * *   Coast   * * * *',q,p)
p = p + 100
hr = floor(time / 3600._wp)
minute = floor(mod(time,3600._wp) / 60._wp)
call messag('Duration of orbital coast:',q,p)
call number(real(hr),-1,z,p)
call messag(':',999,999)
call number(real(minute),-1,999,999)
call messag(':',999,999)
call number(real(mod(time,60._wp)),2,999,999)
time = 0._wp
p = p + 100
call messag('Current orbit eccentricity:',q,p)
call number(real(e(1)),5,z,p)
p = p + 100
call messag('Change in eccentricity through coast:',q,p)
call number(real(abs(e(2) - e(1))),7,z,p)
call sendbf

! Second burn
do i = 1, 5
	xo(i) = coast(k-l,i)
	xn(i) = xo(i)
end do
F = F2
Isp = Isp2
call instant_orbit(xo(3), xo(1), xo(2), u, a(2), e(2), r_perigee, r_apogee)
h = 0
open(unit = 30, status = 'scratch')

do while (a(2) < a_meo - tol .or. a(2) > a_meo + tol)
	
	do i = 1, 5
		xo(i) = xn(i)
	end do
	
	if (method == 't' .or. method == 'T') then
		call trapezoidal_rule(xo, xn, dt, F, Isp, u)
	else
		call runge_kutta(xo, xn, dt, F, Isp, u)
	end if
	
	call instant_orbit(xn(3), xn(1), xn(2), u, a(2), e(2), r_perigee, r_apogee)
	if (a(2) > a_meo + tol) then
		do i = 1, 5
			xn(i) = xo(i)
		end do
		dt = dt / 2._wp
	else
		time = time + dt
		h = h + 1
		write(30,*) xn(1), xn(2), xn(3), xn(4), xn(5)
	end if
	
	if (xn(5) <= 1000._wp) then
		write(*,*)
		write(*,*) ' * * * *   You ran out of fuel!   * * * * '
		exit
	end if
end do

allocate(burn2(h,5))
rewind(30)
do i = 1, h
	read(30,*) burn2(i,1), burn2(i,2), burn2(i,3), burn2(i,4), burn2(i,5)
end do
close(30)
dt = temp1
m2 = coast(k-l,5) - burn2(h,5)
call delta_v(coast(k-l,5),burn2(h,5),Isp,g,dv2)
m = burn2(h,5)
dv = dv1 + dv2

p = p + 100
p = p + 100
call messag('* * * *   Second Burn   * * * *',q,p)
p = p + 100
call messag('Mass used during second burn (kg):',q,p)
call number(real(m2),5,z,p)
p = p + 100
call messag('Mass after second burn (kg):',q,p)
call number(real(burn2(j,5)),5,z,p)
p = p + 100
call messag('Delta V form second burn (km/s):',q,p)
call number(real(dv2),5,z,p)
p = p + 100
call messag('Obrit eccentricity right after burn:',q,p)
call number(real(e(2)),7,z,p)
p = p + 100
call messag('Duration of second burn:',q,p)
hr = floor(time / 3600._wp)
call number(real(hr),-1,z,p)
call messag(':',999,999)
minute = floor(mod(time,3600._wp) / 60._wp)
call number(real(minute),-1,999,999)
call messag(':',999,999)
call number(real(mod(time,60._wp)),2,999,999)
call sendbf





hr = floor((time + temp2 + temp3) / 3600._wp)
minute = floor(mod(time + temp2 + temp3,3600._wp) / 60._wp)
p = p + 100
p = p + 100
call messag('* * * *   Totals   * * * *',q,p)
p = p + 100
call messag('Total flight time:',q,p)
call number(real(hr),-1,z,p)
call messag(':',999,999)
call number(real(minute),-1,999,999)
call messag(':',999,999)
call number(real(mod(time + temp2 + temp3,60._wp)),2,999,999)
p = p + 100
call messag('Final mass: ',q,p)
call number(real(m),5,z,p)
p = p + 100
call messag('Total delta V used: ', q, p)
call number(real(dv),5,z,p)
call sendbf


! Display simulation
a(1) = a_leo
a(2) = a(1)
e(1) = 0._wp
e(2) = e(1)

d = floor(pi * a_meo ** (3._wp/2._wp) / sqrt(u) / dt / 2._wp)
! b = (j + k) / 250










do i = 1, j
	if (i < j .and. i > 2 .and. mod(real(i),real(b)) == 0.) then
		call instant_orbit(burn1(i,3), burn1(i,1), burn1(i,2), u, a(2), e(2), r_perigee, r_apogee)
		call orbits(a(1), a(2), e(1), e(2), burn1(i-b,3), burn1(i,3), burn1(i-b,4), burn1(i,4))
		if (i > d) then
			call polar2cart_graf(burn1(i-d:i,3), burn1(i-d:i,4), d)
		else
			call polar2cart_graf(burn1(1:i,3), burn1(1:i,4), i)
		end if
		call sendbf
		a(1) = a(2)
		e(1) = e(2)
		c = 0
	end if
	c = c + 1
	if (i == j) then
		call instant_orbit(burn1(i,3), burn1(i,1), burn1(i,2), u, a(2), e(2), r_perigee, r_apogee)
		call orbits(a(1), a(2), e(1), e(2), burn1(i-c,3), burn1(i,3), burn1(i-c,4), burn1(i,4))
		call polar2cart_graf(burn1(1:i,3), burn1(1:i,4), i)
		call sendbf
		a(1) = a(2)
		e(1) = e(2)
	end if
end do
call color('yellow')
do i = 1, k - l
	if (i > b .and. mod(real(i),real(b)) == 0.) then
		call polar2cart_graf(coast(i-b:i,3), coast(i-b:i,4),b)
		call sendbf
	end if
	if (i == k - l) then
		call polar2cart_graf(coast(1:i,3), coast(1:i,4),i)
		call sendbf
	end if
end do
call color('fore')
do i = 1, h
	if (i < h .and. i > 2 .and. mod(real(i),real(b)) == 0.) then
		call instant_orbit(burn2(i,3), burn2(i,1), burn2(i,2), u, a(2), e(2), r_perigee, r_apogee)
		call orbits(a(1), a(2), e(1), e(2), burn2(i-b,3), burn2(i,3), burn2(i-b,4), burn2(i,4))
		if (i > d) then
			call polar2cart_graf(burn2(i-d:i,3), burn2(i-d:i,4), d)
		else
			call polar2cart_graf(burn2(1:i,3), burn2(1:i,4), i)
		end if
		call color('yellow')
		call polar2cart_graf(coast(:,3), coast(:,4), k-l)
		call color('fore')
		call polar2cart_graf(burn1(:,3), burn1(:,4), j)
		call sendbf
		a(1) = a(2)
		e(1) = e(2)
		c = 0
	end if
	c = c + 1
	if (i == h) then
		call instant_orbit(burn2(i,3), burn2(i,1), burn2(i,2), u, a(2), e(2), r_perigee, r_apogee)
		call orbits(a(1), a(2), e(1), e(2), burn2(i-c,3), burn2(i,3), burn2(i-c,4), burn2(i,4))
		call polar2cart_graf(burn2(1:i,3), burn2(1:i,4), i)
		call color('yellow')
		call polar2cart_graf(coast(1:k-l,3), coast(1:k-l,4), k-l)
		call color('fore')
		call polar2cart_graf(burn1(:,3), burn1(:,4), j)
		call sendbf
		a(1) = a(2)
		e(1) = e(2)
	end if
end do

write(*,*)
write(*,*) ' * * *   Simulation Complete   * * *'

call disfin

write(*,*)
write(*,*) 'Would you like to repeat the program (y/n): '
read(*,*) cont

deallocate(burn1, coast, burn2)
end do

end program

subroutine polar2cart_graf(r,theta,n)
integer, intent(in) :: n
integer, parameter :: wp = selected_real_kind(p=12)
real(wp), intent(in), dimension(n) :: r, theta
real, dimension(n) :: x, y
integer :: i

do i = 1, n
	x(i) = real(r(i) * cos(theta(i)))
	y(i) = real(r(i) * sin(theta(i)))
end do

call curve(x,y,n)

end subroutine polar2cart_graf

subroutine orbits(a1, a2, e1, e2, r1, r2, v_inert1, v_inert2)
integer, parameter :: wp = selected_real_kind(p=12)
real(wp), intent(in) :: a1, a2, e1, e2, r1, r2, v_inert1, v_inert2
real(wp) :: v1, v2
real(wp), parameter :: a_leo = 8530._wp, a_meo = 13200._wp, r_earth = 6371._wp	! km
real(wp), parameter :: pi = 3.141592654_wp
real(wp), dimension(361) :: earth_r, earth_v, leo_r, meo_r, r_1, r_2, v_1, v_2
integer :: i
if (e1 == 0. .or. e2 == 0.) then
	v1 = 0._wp
	v2 = 0._wp
else
	v1 = acos((a1/r1*(1._wp-e1**2)-1._wp)/e1)
	v2 = acos((a2/r2*(1._wp-e2**2)-1._wp)/e2)
end if
earth_r = r_earth
leo_r = a_leo
meo_r = a_meo
do i = 1, 361
	earth_v(i) = dble(i-1) * pi / 180._wp
	v_1(i) = dble(i-1) * pi / 180._wp
	v_2(i) = dble(i-1) * pi / 180._wp
	r_1(i) = a1*(1._wp-e1**2)/(1._wp+e1*cos(v_1(i) + v1 - v_inert1))
	r_2(i) = a2*(1._wp-e2**2)/(1._wp+e2*cos(v_2(i) + v2 - v_inert2))
end do
call color('GREEN')
call polar2cart_graf(earth_r, earth_v, 361)
call color('red')
call dashl
call polar2cart_graf(leo_r, earth_v, 361)
call color('blue')
call dashl
call polar2cart_graf(meo_r, earth_v, 361)
call color('back')
call solid
call polar2cart_graf(r_1, v_1, 361)
call color('orange')
call polar2cart_graf(r_2, v_2, 361)
call color('fore')
end subroutine orbits

subroutine prompt_user(string,var,length)
integer, intent(in) :: length
integer, parameter :: wp = selected_real_kind(p=12)
character(len=length), intent(in) :: string
real(wp), intent(out) :: var
write(*,*)
write(*,*) string
read(*,*) var
end subroutine prompt_user

subroutine instant_orbit(r, Vr, Vv, u, a, e, rp, ra)
integer, parameter :: wp = selected_real_kind(p=12)
real(wp), intent(in) :: r, Vr, Vv, u
real(wp), intent(out) :: a, e, rp, ra
a = u / (2._wp * u / r - Vr**2 - Vv**2)
e = r / u * sqrt((Vv**2 - u / r)**2 + (Vr * Vv)**2)
rp = a * (1._wp - e)
ra = a * (1._wp + e)
end subroutine instant_orbit

subroutine trapezoidal_rule(xo, xn, dt, F, Isp, u)
integer, parameter :: wp = selected_real_kind(p=12)
real(wp), intent(in), dimension(5) :: xo
real(wp), intent(in) :: dt, u, F, Isp
real(wp), intent(out), dimension(5) :: xn
real(wp), dimension(5) :: x_temp, fo, fn
real(wp) :: gamma
real(wp), parameter :: g = 9.806_wp
integer :: i

! x variables
! 1 is radial velocity
! 2 is tangetial velocity
! 3 is radial position
! 4 is the true anomoly
! 5 is mass
! f variable are the time derivitives of x

! f values
gamma = atan(xo(1)/xo(2))
fo(1) = xo(2)**2 / xo(3) + F * sin(gamma) / xo(5) / 1000._wp - u / xo(3)**2
fo(2) =  F * cos(gamma) / xo(5) / 1000._wp - xo(1) * xo(2) / xo(3)
fo(3) = xo(1)
fo(4) = xo(2) / xo(3)
fo(5) = -1._wp * F / g / Isp

! prediction step
do i = 1, 5
	x_temp(i) = xo(i) + dt * fo(i)
end do

! new f values
gamma = atan(x_temp(1)/x_temp(2))
fn(1) = x_temp(2)**2 / x_temp(3) + F * sin(gamma) / x_temp(5) / 1000._wp - u / x_temp(3)**2
fn(2) =  F * cos(gamma) / x_temp(5) / 1000._wp - x_temp(1) * x_temp(2) / x_temp(3)
fn(3) = x_temp(1)
fn(4) = x_temp(2) / x_temp(3)
fn(5) = -1._wp * F / g / Isp

! correction step
do i = 1, 5
	xn(i) = xo(i) + dt / 2._wp * (fo(i) + fn(i))
end do

end subroutine trapezoidal_rule

subroutine runge_kutta(xo, xn, dt, F, Isp, u)
integer, parameter :: wp = selected_real_kind(p=12)
real(wp), intent(in), dimension(5) :: xo
real(wp), intent(in) :: dt, u, F, Isp
real(wp), intent(out), dimension(5) :: xn
real(wp), dimension(5) :: k1, k2, k3, k4, x_temp
real(wp) :: gamma
real(wp), parameter :: g = 9.806_wp
integer :: i

! k1
gamma = atan(xo(1)/xo(2))
k1(1) = xo(2)**2 / xo(3) + F * sin(gamma) / xo(5) / 1000._wp - u / xo(3)**2
k1(2) =  F * cos(gamma) / xo(5) / 1000._wp - xo(1) * xo(2) / xo(3)
k1(3) = xo(1)
k1(4) = xo(2) / xo(3)
k1(5) = -1._wp * F / g / Isp

! x by k1
do i = 1, 5
	x_temp(i) = xo(i) + dt / 2._wp * k1(i)
end do

! k2
gamma = atan(x_temp(1)/x_temp(2))
k2(1) = x_temp(2)**2 / x_temp(3) + F * sin(gamma) / x_temp(5) / 1000._wp - u / x_temp(3)**2
k2(2) =  F * cos(gamma) / x_temp(5) / 1000._wp - x_temp(1) * x_temp(2) / x_temp(3)
k2(3) = x_temp(1)
k2(4) = x_temp(2) / x_temp(3)
k2(5) = -1._wp * F / g / Isp

! x by k2
do i = 1, 5
	x_temp(i) = xo(i) + dt / 2._wp * k2(i)
end do

! k3
gamma = atan(x_temp(1)/x_temp(2))
k3(1) = x_temp(2)**2 / x_temp(3) + F * sin(gamma) / x_temp(5) / 1000._wp - u / x_temp(3)**2
k3(2) =  F * cos(gamma) / x_temp(5) / 1000._wp - x_temp(1) * x_temp(2) / x_temp(3)
k3(3) = x_temp(1)
k3(4) = x_temp(2) / x_temp(3)
k3(5) = -1._wp * F / g / Isp

! x by k3
do i = 1, 5
	x_temp(i) = xo(i) + dt * k3(i)
end do

! k4
gamma = atan(x_temp(1)/x_temp(2))
k4(1) = x_temp(2)**2 / x_temp(3) + F * sin(gamma) / x_temp(5) / 1000._wp - u / x_temp(3)**2
k4(2) =  F * cos(gamma) / x_temp(5) / 1000._wp - x_temp(1) * x_temp(2) / x_temp(3)
k4(3) = x_temp(1)
k4(4) = x_temp(2) / x_temp(3)
k4(5) = -1._wp * F / g / Isp

! special runge_kutta equation
do i = 1, 5
	xn(i) = xo(i) + dt / 6._wp * (k1(i) + k2(i) + k3(i) + k4(i))
end do

end subroutine runge_kutta

subroutine delta_v(m_initial, m_final, Isp, g0, dv)
integer, parameter :: wp = selected_real_kind(p=12)
real(wp), intent(in) :: m_initial, m_final, Isp, g0
real(wp), intent(out) :: dv
dv = g0 * Isp * log( m_initial / m_final) / 1000._wp
end subroutine delta_v
