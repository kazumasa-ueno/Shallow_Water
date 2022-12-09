module transfer_mod
	use boundary_mod
	use calc_variables_mod
	implicit none
	
contains

	! I_k^{k-1}
	subroutine Prolongation(Auf,Avf,Azf,Auc,Avc,Azc,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny !course grid number
		real(8), intent(in) :: Auf(0:2*Nx,1:2*Ny), Avf(1:2*Nx,0:2*Ny), Azf(1:2*Nx,1:2*Ny)
		real(8), intent(out) :: Auc(0:Nx,1:Ny), Avc(1:Nx,0:Ny), Azc(1:Nx,1:Ny)

		integer :: ic, jc, iff, jff

		do jc = 1, Ny
			do ic = 1, Nx
				iff = 2*ic
				jff = 2*jc
				Auc(ic,jc) = (Auf(iff,jff) + Auf(iff,jff-1))/2.d0
				Avc(ic,jc) = (Avf(iff,jff) + Avf(iff-1,jff))/2.d0
				! Azc(ic,jc) = (Azf(iff,jff) + Azf(iff-1,jff) + Azf(iff,jff-1) + Azf(iff-1,jff-1))/4.d0
			end do
			Auc(0,jc) = (Auf(0,jff) + Auf(0,jff-1))/2.d0
		end do
		do ic = 1, Nx
			Avc(ic,0) = (Avf(iff,0) + Avf(iff-1,0))/2.d0
		end do
		call calc_Az(Auc,Avc,Nx,Ny,Azc)

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