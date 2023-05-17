module structure
	use constant
	implicit none
	
	real(8), dimension(Nx,num_levels) :: u, z, h, XForce, fmgz
  real(8), dimension(Nx,num_levels) :: Au, Az, b,residual, rhs !係数

end module structure