module sparse_mod
implicit none




contains



subroutine read_matrix(unit_id, ndim, R, C, V)
implicit none
integer :: ndim, i, notz, unit_id, i_stat
integer, allocatable :: R(:), C(:)
double precision, allocatable :: V(:)


! reading dimension of the matrix
read(unit_id,*)ndim
write(*,*)"Matrix dimension=",ndim

! R array allocation and read
allocate(R(ndim+1), stat=i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if
read(unit_id,*) R
write(*,*) R

! determine the lenght of C and V (the last value of R is the lenght)
notz=R(ndim+1)

! C array allocation and read
allocate(C(notz), stat=i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if
read(unit_id,*) C
write(*,*) C

! V array allocation and read
allocate(V(notz), stat=i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if
read(unit_id,*) V
write(*,*) V

end subroutine read_matrix

subroutine sparse_multiplication(R1, C1, V1, R2, C2, V2, C, n_mul, ndim)
implicit none
integer :: p1, p2, i, j
double precision :: y 
integer, allocatable, intent(in) :: R1(:), C1(:), R2(:), C2(:)
double precision, allocatable, intent(in) :: V1(:), V2(:)
integer, intent(in) :: ndim

integer, intent(out) :: n_mul
double precision, intent(out) :: C(:,:)

! inizializing the product matrix and the number of multiplication
C=0.D0
n_mul=0

! multiplication
do j=1, ndim        !rows of B
 do i=j, ndim    !rows of A  (exploiting the simmetry  of the matrices)
  y = 0.0d0    !accumulator for the new matrix element
  do p1=R1(i)+1,R1(i+1)   !number of non zero element per row of A
   do p2=R2(j)+1,R2(j+1)    !number of non zero element per row of B 
    if (C1(p1) == C2(p2)) then
     y=y+V1(p1)*V2(p2)
     n_mul=n_mul+1 
    endif
   enddo
  enddo
  C(i,j) = y
  C(j,i) = y
 enddo
enddo



end subroutine sparse_multiplication




end module sparse_mod
