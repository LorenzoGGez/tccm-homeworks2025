program mol_dyn
implicit none
integer :: Natoms,m,n,input_file,read_Natoms

! declaring array 
double precision, allocatable :: xyz(:,:)
! check allocatiosnn status
integer :: i_stat




! allocation 2D array
allocate(xyz(m,n) , stat=i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if



! deallocate the array
deallocate(xyz)



Natoms=read_Natoms(input_file)  







end program mol_dyn







































































! reading atom function
integer function read_Natoms(input_file) result(Natoms)
implicit none
integer, intent(in) :: input_file
read(input_file,*)Natoms
write(6,*)"ATOMS NUMBER =", Natoms
end function read_Natoms 


! distances calculations
subroutine compute_distances(Natoms, coord, dist)
implicit none
integer, intent(in) :: Natoms
double precision, intent(in) :: coord(Natoms,3)
double precision, intent(out) :: dist(Natoms,Natoms)
integer :: i,j
double precision :: xdist,ydist,zdist

do i=1, Natoms
 do j=i, Natoms
  xdist=coord(i,1)-coord(j,1)
  ydist=coord(i,2)-coord(j,2)
  zdist=coord(i,3)-coord(j,3)
  dist(i,j)=sqrt(xdist**2+ydist**2+zdist**2)          ! calculating the distance between two atom
  dist(j,i)=dist(i,j)                                 ! filling the symmetric elements
 enddo
enddo
end subroutine compute_distances





! reading atomic data
subroutine read_mol(input_file, n, c, m)
implicit none
integer, intent(in) :: input_file
integer :: i
integer, intent(in) :: n
doubleprecision, intent(out) :: c(n,3)
doubleprecision, intent(out) :: m(n)
do i=1,n
 read(*,*) c(i,1), c(i,2), c(i,3), m(i)
end do
end subroutine read_mol

double precision function v(eps, sig, n, dis)
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
 end do
end do
end function v


! calculating the kinetic energy
double precision function T(Natoms, vel, mass)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: vel(Natoms,3)
        double precision, intent(in) :: mass(Natoms)
        double precision :: xvel, yvel, zvel
        integer :: i
do i=1, Natoms
 xvel=vel(i,1)**2                       ! calculating the square of each component of the velocity
 yvel=vel(i,2)**2                       !        
 zvel=vel(i,3)**2                       !
 T=T+mass(i)*(xvel+yvel+zvel)           ! updating the total kinetic energy
enddo
T=0.5D0*T
end function T









































































































