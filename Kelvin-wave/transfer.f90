!********************************************
! ProlongationやInterpolationを計算
!********************************************

module transfer_mod
	use boundary_mod
	use calc_variables_mod
	implicit none
	
contains

	! I_k^{k-1}
	subroutine Prolongation(zf,gammaf,hf,zc,gammac,hc,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny !course grid number
		real(8), intent(in) :: zf(1:2*Nx,1:2*Ny), gammaf(1:2*Nx,1:2*Ny), hf(1:2*Nx,1:2*Ny)
		real(8), intent(out) :: zc(0:Nx+1,0:Ny+1), gammac(0:Nx+1,0:Ny+1), hc(0:Nx+1,0:Ny+1)

		integer :: ic, jc, iff, jff

		do jc = 1, Ny
			do ic = 1, Nx
				iff = 2*ic
				jff = 2*jc
				zc(ic,jc) = (zf(iff,jff) + zf(iff-1,jff) + zf(iff,jff-1) + zf(iff-1,jff-1))*0.25d0
				gammac(ic,jc) = (gammaf(iff,jff) + gammaf(iff-1,jff) + gammaf(iff,jff-1) + gammaf(iff-1,jff-1))*0.25d0
				hc(ic,jc) = (hf(iff,jff) + hf(iff-1,jff) + hf(iff,jff-1) + hf(iff-1,jff-1))*0.25d0
			end do
		end do

		call boundary_z(zc,Nx,Ny)
		call boundary_z(gammac,Nx,Ny) !暫定
		call boundary_z(hc,Nx,Ny) !暫定
		! zのboundaryは外側で実装する
		! gammaのboundaryは外側で実装する
		! hのboundaryは外側で実装する

	end subroutine Prolongation

	subroutine Prolongation_defect(df,dc,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny !course grid number
		real(8), intent(in) :: df(1:2*Nx,1:2*Ny)
		real(8), intent(out) :: dc(1:Nx,1:Ny)

		integer :: ic, jc, iff, jff

		do jc = 1, Ny
			do ic = 1, Nx
				iff = 2*ic
				jff = 2*jc
				dc(ic,jc) = (df(iff,jff) + df(iff,jff-1) + df(iff-1,jff) + df(iff-1,jff-1))/4.d0
			end do
		end do

	end subroutine Prolongation_defect

	!I_k^{k+1}
	subroutine Interpolation_defect(dc,df,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny !fine grid number
		real(8), intent(in)  :: dc(0:Nx/2+1,0:Ny/2+1)
		real(8), intent(out) :: df(0:Nx+1,0:Ny+1)

		integer :: ic, jc, iff, jff
		real(8) :: dc1, dc2, dc3, dc4, dc5

		do jff = 2, Ny, 2
			do iff = 2, Nx, 2
				ic = iff/2
				jc = jff/2
				dc1 = dc(ic,jc)
				dc2 = dc(ic-1,jc-1)
				dc3 = dc(ic-1,jc+1)
				dc4 = dc(ic+1,jc-1)
				dc5 = dc(ic+1,jc+1)
				df(iff  ,jff  ) = (3*dc1+dc5)*0.25
				df(iff  ,jff-1) = (3*dc1+dc3)*0.25
				df(iff-1,jff  ) = (3*dc1+dc4)*0.25
				df(iff-1,jff-1) = (3*dc1+dc2)*0.25
			end do
		end do
		
	end subroutine Interpolation_defect
end module transfer_mod