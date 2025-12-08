program mol_dyn











end program mol_dyn


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
