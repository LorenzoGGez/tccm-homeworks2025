! function calculating the energy potential derivate
double precision function dpot(eps, sig, r) result(u)
        implicit none
        double precision, intent(in) :: eps, sig
        double precision, intent(in) :: r
        double precision :: d,l,p

d=24.D0*eps/r
l=2.D0*(sig/r)**12
p=(sig/r)**6
u=d*(p-l)
end function dpot



double precision function V(eps, sig, n, dis)
        implicit none
        integer ::  i, t
        double precision, intent(in) :: eps, sig
        integer, intent(in) :: n
        double precision, intent(in) :: dis(n,n)
        double precision :: d,l,p
        
v=0.D0
do i=1,n
 do t=i+1,n
  d=4.D0*eps
  l=(sig/dis(i,t))**12
  p=(sig/dis(i,t))**6
  v=v+d*(l-p)
!  write(*,*)"l=",l,"p=",p,"d=",d,"v=",v,"           i=",i,"t=",t     debug line for checking value of the component of potential energy
 end do
end do
V=v
end function V

! calculating the kinetic energy
double precision function T(Natoms, vel, mass)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: vel(Natoms,3)
        double precision, intent(in) :: mass(Natoms)
        double precision :: xvel, yvel, zvel
        integer :: i
t=0.D0
do i=1, Natoms
 xvel=vel(i,1)**2                       ! calculating the square of each component of the velocity
 yvel=vel(i,2)**2                       !        
 zvel=vel(i,3)**2                       !
 t=t+mass(i)*(xvel+yvel+zvel)           ! updating the total kinetic energy
enddo
T=0.5D0*t
end function T
