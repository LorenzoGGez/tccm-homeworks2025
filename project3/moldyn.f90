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
