program main
implicit none
integer :: ind,i_stat
integer, parameter :: n=2
double precision :: eps, sig, dist, u, mass, V, T, dpot
double precision, allocatable :: dis(:,:), vel(:,:)
write(*,*)"Chose the function to test:"
write(*,*)"1-POTENTIAL ENERGY"
write(*,*)"2-KINETIC ENERGY"
write(*,*)"3-POTENTIAL ENERGY DERIVATE"
read(*,*)ind


if(ind == 1) then
 allocate(dis(n,n) , stat=i_stat)             ! allocation distance array
  if(i_stat /= 0) then
   print *, "MEMORY ALLOCATION FAILED"
   stop
 end if
 write(*,*)"Calculating LJ energy for two atoms"
 write(*,*)"Enter the sigma (Angstrom):"
 read(*,*)sig
 write(*,*)"Enter the epsilon (kj/mol):"
 read(*,*)eps
 write(*,*)"Enter the distance (Angstrom):"
 read(*,*)dis(1,2)
 dis(1,2)=dis(1,2)+10.D0
 dis(2,1)=dis(1,2)                              !conversion to nm
 sig=sig*10.D0                                  !conversion to nm
 u=V(eps,sig,n,dis)
 write(*,*)"The value of the potential energy for the given data is:",u,"kj/mol"
 deallocate(dis)
elseif(ind == 2) then
 allocate(vel(1,3) , stat=i_stat)             ! allocation vel  array
  if(i_stat /= 0) then
   print *, "MEMORY ALLOCATION FAILED"
   stop
 end if
 write(*,*)"Calculating kinetic energy for one atom"
 write(*,*)"Enter the three component of the speed:"
 read(*,*)vel(1,1),vel(1,2),vel(1,3)
 write(*,*)"Enter the mass of the atom:"
 read(*,*)mass
 u=T(1,vel,mass)
 write(*,*)"The value of the kinetic  energy for the given data is:",u
 deallocate(vel)
elseif(ind == 3) then
 write(*,*)"Enter the sigma (Angstrom):"
 read(*,*)sig
 write(*,*)"Enter the epsilon (kj/mol):"
 read(*,*)eps
 write(*,*)"Enter the distance (Angstrom):"
 read(*,*)dist
 dist=dist*10.D0                                !conversion to nm
 sig=sig*10.D0                                !conversion to nm
 u=dpot(eps, sig, dist)
 write(*,*)"The value of the potential energy derivate for the given data is:",u
else
 write(*,*)"OPTION NOT FOUND, TRY AGAIN"
 stop
endif



end program main
