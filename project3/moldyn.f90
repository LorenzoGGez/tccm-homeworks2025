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

double precision function v(eps, sig, n, dis)
        implicit none
        integer :: d, l, p, i, t
        double precision, intent(in) :: eps, sig
        integer, intent(in) :: n
        double precision, intent(in) :: dis(n,n)

do i=1,n
do t=i+1,n
d=4*eps
l=(sig/dis(i,t))**12
p=(sig/dis(i,t))**6
v=d*(l-p)
end do
end do
end function v











































































































