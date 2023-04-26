!********************************************
! ProlongationやInterpolationを計算
!********************************************

module transfer_mod
	use boundary_mod
	use calc_variables_mod
	implicit none
	
contains

	! I_k^{k-1}
	subroutine Prolongation(Auf,Azf,Auc,Azc,Nx)
		implicit none

		integer, intent(in) :: Nx !course grid number
		real(8), intent(in) :: Auf(0:2*Nx), Azf(1:2*Nx)
		real(8), intent(out) :: Auc(0:Nx), Azc(1:Nx)

		integer :: ic, iff

		do ic = 0, Nx
			iff = 2*ic
			Auc(ic) = Auf(iff)
			! Azc(ic,jc) = (Azf(iff,jff) + Azf(iff-1,jff) + Azf(iff,jff-1) + Azf(iff-1,jff-1))/4.d0
		end do

		call calc_Az(Auc,Nx,Azc)

	end subroutine Prolongation

	subroutine Prolongation_defect(df,dc,Nx)
		implicit none

		integer, intent(in) :: Nx !course grid number
		real(8), intent(in) :: df(1:2*Nx)
		real(8), intent(out) :: dc(1:Nx)

		integer :: ic, jc, iff, jff

		do ic = 1, Nx
			iff = 2*ic
			dc(ic) = (df(iff) + df(iff-1))*0.5d0
		end do

	end subroutine Prolongation_defect

	!I_k^{k+1}
	subroutine Interpolation_defect(dc,df,Nx)
		implicit none

		integer, intent(in) :: Nx !fine grid number
		real(8), intent(in)  :: dc(0:Nx/2+1)
		real(8), intent(out) :: df(0:Nx+1)

		integer :: ic, iff
		real(8) :: dc1, dc2, dc3

		do iff = 2, Nx, 2
			ic = iff/2
			dc1 = dc(ic)
			dc2 = dc(ic-1)
			dc3 = dc(ic+1)
			df(iff) = (3*dc1+dc3)*0.25
			df(iff-1) = (3*dc1+dc2)*0.25
		end do
		
	end subroutine Interpolation_defect
end module transfer_mod