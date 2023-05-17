module output
	use constant
	use structure
	implicit none
	
contains
	subroutine outucd(times)
		implicit none

		integer, intent(in) :: times
		integer :: l
		
		open(unit=10, file='../output/z'//trim(adjustl(itoa(times)))//'.dat', status='unknown')
		do l = 1, num_levels
			write(10,*)  z(:,l)
		enddo
		close(10)
		open(unit=20, file='../output/fmgz'//trim(adjustl(itoa(times)))//'.dat', status='unknown')
		do l = 1, num_levels
			write(20,*)  fmgz(:,l)
		enddo
		close(20)
	end subroutine outucd

	function itoa(i)
    implicit none
    integer, intent(in) :: i
    character(len=32) :: itoa
    write(itoa, '(i0)') i
    itoa = adjustl(itoa)
  end function itoa

end module output