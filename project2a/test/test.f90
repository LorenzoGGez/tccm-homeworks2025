program main
use sparse_mod
implicit none
integer :: i, j, ndim1, i_stat
double precision :: y
integer, allocatable :: R1(:), C1(:)
double precision, allocatable :: V1(:), mat1(:,:)
character(len=50) :: filename


! first matrix entering
write(*,*)"Insert the name of the matrix file"
read(*,*)filename
open(unit=10,file=trim(filename), status='old')
call read_matrix(10, ndim1, R1, C1, V1)
close(10)
write(*,*)



! FIRST MATRIX
! allocating the non-sparse format matrix
allocate(mat1(ndim1,ndim1),  stat= i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if

! turning spare matrix in non-sparse format
call desparsification(R1, C1, V1, mat1, ndim1)

write(*,*)"---------------------------------------------------------------------------------"
write(*,*)"Transforming the first matrix in the non-sparse format"
write(*,*)"---------------------------------------------------------------------------------"


do i=1, ndim1
 write(*,'(*(F12.6,1X))')mat1(i,:)                    
enddo
write(*,*)"---------------------------------------------------------------------------------"


deallocate(R1, C1, V1, mat1)
end program main


