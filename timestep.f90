module timestep
	implicit none
	
contains

	subroutine timesteps(nx, ny, dx, dy, dt, g, A, Height, f, u, v, h)

		integer, intent(in) :: nx, ny
		real(8), intent(in) :: dy, dt, g, A, Height
		real(8), intent(in) :: dx(:), f(:)
		real(8), intent(inout) :: u(:,:,:), v(:,:,:), h(:,:,:)
		real(8) :: term1, term2, term3
		integer :: i, j

!$acc kernels present(u,v,h)
!$acc loop independent
		do j = 2, ny-1
!$acc loop independent
			do i = 2, nx-1
				!calculating u
				term1 = f(j) * (v(i-1,j,2)+v(i-1,j+1,2)+v(i,j,2)+v(i,j+1,2))/4 !!後で重み付き平均に変更
				term2 = -g * (h(i,j,2)-h(i-1,j,2))/dx(j)
				term3 = A * ((u(i+1,j,1)-2*u(i,j,1)+u(i-1,j,1))/(dx(j))**2 + (u(i,j+1,1)-2*u(i,j,1)+u(i,j-1,1))/(dy)**2)
				u(i,j,3) = u(i,j,1) + 2*dt*(term1+term2+term3)

				!calculating v
				term1 = -f(j) * (u(i,j-1,2)+u(i+1,j-1,2)+u(i,j,2)+u(i+1,j,2))/4 !!後で重み付き平均に変更
				term2 = -g * (h(i,j,2)-h(i,j-1,2))/dy
				term3 = A * ((v(i+1,j,1)-2*v(i,j,1)+v(i-1,j,1))/(dx(j))**2 + (v(i,j+1,1)-2*v(i,j,1)+v(i,j-1,1))/(dy)**2)
				v(i,j,3) = v(i,j,1) + 2*dt*(term1+term2+term3)

				!calculating h
				h(i,j,3) = h(i,j,1) + 2*dt*(-Height*((u(i+1,j,2)-u(i,j,2))/dx(j) + (v(i,j+1,2)-v(i,j,2))/dy))
			end do
		end do
!$acc end kernels

	end subroutine timesteps

end module timestep