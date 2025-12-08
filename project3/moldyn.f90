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
double precision :: intent(in) :: coord(Natoms,3)
double precision :: intent(out) :: coord(Natoms,Natoms)








