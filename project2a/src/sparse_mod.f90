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

end module sparse_mod
