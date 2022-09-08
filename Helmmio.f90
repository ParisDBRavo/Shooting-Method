Module numbers
!arreglos
Real(kind=8), allocatable, dimension(:) :: t, rhs
Real(kind=8), allocatable, dimension(:, :) :: u
Real(kind=8), allocatable, dimension(:) :: k1, k2, k3, k4

! Parameters
Real(kind=8) x0, y0

! Method and code variables
Real(kind=8) t0, tmax, dt
Integer :: N, resolution_label, NE

! Shooting parameters
Real(kind=8) x1_p, epsilon, kk0, dkk, kk
Integer k_it_max, num_of_nodes, zeroes

End module 

! - - - - - - - 

Program rk4_Helmholtz

Use numbers
Implicit none

! Counters defined locally
Integer i, j, l

Resolution_label =1
NE = 4
N = 100
T0 = 0.0d0
Tmax = 1.0d0
N = 2**(resolution_label-1) *N
! Number of cells for the domain 
X0 = 0.0d0
Y0 = -1.1d0
kk0 = -0.5
dkk=0.1
k_it_max = 10000
epsilon = 1.e-8
num_of_nodes = 0.0d0

!Allocate memory for the various arrays
!t:time
!u(1) : x
!u(2) : y
allocate(t(0:N), u(1:NE,0:N), rhs(1:NE),k1(1:NE),k2(1:NE), k3(1:NE),k4(1:NE))
!------PART A------
dt = (tmax-t0) /dble(N)

do i=0, N
    t(i) = t0+ dt*dble(i)
end do

!------PART B------
open(1,file = 'Helmholtz.dat')

u(1,0) = x01
u(2,0) = y01
u(3,0) = x02
u(4,0) = y02

kk = kk0
zeroes = 0

do l=1, k_it_max
  
  do i=1, N
    do j=1, 4
        if(j.eq.1)then 
            call calcrhs(t(i-1), u(:,i-1))
            k1 = rhs        
        else if(j.eq.2)then 
            call calcrhs(t(i-1)+0.5d0*dt, u(:,i-1)+0.5d0*k1(:)*dt)
            k2 = rhs
        else if(j.eq.3)then 
            call calcrhs(t(i-1)+0.5d0*dt, u(:,i-1)+0.5d0*k2(:)* dt)
            k3 = rhs
        else if(j.eq.4)then 
            call calcrhs(t(i-1)+dt, u(:,i-1)+k3(:)* dt)
            k4 = rhs    
            u(:,i) = u(:,i-1)+(1.0d0/6.0d0)*(k1+2.0d0*k2+2.0d0*k3+k4)*dt
        end if
    end do
    if(y0.ge.0)then
    if((u(1,i)*u(1,i-1).le.0.0d0).and.(i.ge.2).and.(i.le.N-2))then
        zeroes = zeroes +1 
    
    end if





    end if

    
  end do
!Changing k
x1_p= U(1,N)
print *,l-1,kk, U(1,N)



if(y0.ge.0)then
if((x1_p.lt.0.0).and.(zeroes.eq.num_of_nodes)) then
    kk= kk-dkk
    dkk = 0.5d0*dkk
else
    kk = kk+dkk
end if
if(abs(u(1,N)).lt.epsilon) then
    print *,'Best k =', kk
    do i=0, N, 2**resolution_label
        write(1,*) t(i), u(1,i), u(2,i)    

    end do
    write(1,*)
    write(1,*)
    exit
end if


!Modificación mia
else if(y0.lt.0)then

if((x1_p.gt.0.0).and.(zeroes.eq.num_of_nodes)) then

    kk= kk-dkk
    dkk = 0.5d0*dkk
else

    kk = kk+dkk
end if

if(abs(u(1,N)).lt.epsilon) then
    print *,'Best k =', kk
    do i=0, N, 2**resolution_label
        write(1,*) t(i), u(1,i), u(2,i)    

    end do
    write(1,*)
    write(1,*)
    exit
end if

end if


end do

 close(1)
end program


subroutine calcrhs(my_t,my_u)
use numbers
implicit none

real(kind=8), intent(in) :: my_t
real(kind=8), dimension(NE), intent(in) ::  my_u

    rhs(1) =my_u(2)
    rhs(2) =-kk**2 * my_u(1)
    rhs(3) =
    rhs(4) =

end subroutine calcrhs
