!********************************************
!	中央部にアノマリーを置いたときのロスビー波や重力波の伝播
! 浅水波方程式をマルチグリッド法で解くプログラム
! Semi-implicit Semi-Lagtange法を使用して離散化
!********************************************

program main
	use boundary_mod
	use calc_variables_mod
	use transfer_mod
	use,intrinsic :: iso_fortran_env
	implicit none

	integer(int32) :: time_begin_c,time_end_c, CountPerSec, CountMax !時間測定用
	
	integer, parameter :: l = 2								!グリッドの深さ
	integer, parameter :: Nx = 64
	integer, parameter :: ntmax = 5000				!時間ステップ
	integer, parameter :: nu1 = 2, nu2 = 1		!マルチグリッドサイクル内のsmooth回数
	real(8), parameter :: g = 9.81d0 					!重力定数
	real(8), parameter :: Cz = 80.d0					!Chezy 摩擦係数
	real(8), parameter :: pi = 4*atan(1.d0)		!円周率
	real(8), parameter :: f0 = 4*pi/86400			!コリオリパラメータf0
	real(8), parameter :: X = 1.d8						!領域サイズ
	! real(8), parameter :: dt = 60.d0*4				!時間間隔
	real(8), parameter :: dt = 90.d0					!時間間隔
	real(8), parameter :: dtau = dt/10.d0			!移流計算用小時間間隔

	real(8) :: f 	!コリオリパラメータ
	real(8) :: h(0:Nx+1)	!基準面からの水深
	real(8) :: u(0:Nx), v(0:Nx+1), z(0:Nx+1), gamma(0:Nx+1), u_b(0:Nx), v_b(0:Nx+1) !u_b, v_bは移流計算実行用の一時格納配列
	real(8) :: Au(0:Nx), Az(1:Nx), b(1:Nx) !係数
	real(8) :: Fu(Nx), Fv(Nx)
	real(8) :: dx	!格子間隔
	real(8) :: Res, difference	!Resは残差のl2ノルム、differenceは前の時間との残差	
	integer :: times, cyc !時間ループ用と収束までの繰り返し用
	integer :: ios !ファイル開く用
	integer :: i, j !空間ループ用
	real(8) :: Fu_tmp, Fv_tmp

	!for debug
	real(8) :: Prev(0:Nx+1), tmp(Nx+2) !前の値を格納しておくための配列

	! open(unit=10, file="./output/u.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/u.txt"
	! open(unit=11, file="./output/v.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/v.txt"
	open(unit=12, file="./output/z.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z.txt"
	open(unit=20, file="./output/z1.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	open(unit=21, file="./output/z2.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z2.txt"
	open(unit=22, file="./output/z3.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z3.txt"
	open(unit=23, file="./output/z4.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	open(unit=24, file="./output/z5.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z2.txt"
	open(unit=25, file="./output/z6.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z3.txt"
	open(unit=26, file="./output/z7.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	open(unit=27, file="./output/z8.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	! open(unit=30, file="./output/res5.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/res.txt"
	

	!!コリオリパラメータfと格子間隔dxの計算
	f = f0
	dx = X/Nx

	!u,v,z,gamma,hの初期化
	call initialize(u,v,z,gamma,h,dx,Nx)

	!メインのループ
	do times = 0, ntmax-1
		!時間計測スタート
		call system_clock(time_begin_c, CountPerSec, CountMax)
		
		!係数計算
		call calc_Au(z,gamma,h,g,dt,dx,Nx,Au)
		call calc_Az(Au,Nx,Az)
		call calc_b(u,v,z,gamma,h,f,dt,dx,dtau,Nx,b)
		
		! do j = 1, Ny
		! 	Fv(1,j) = calc_Fv(1,j,u,v,dt,dx,dy,dtau,Nx,Ny)
		! 	Fv(2,j) = calc_Fv(Nx,j,u,v,dt,dx,dy,dtau,Nx,Ny)
		! end do
		! do i = 1, Nx
		! 	Fu(i,1) = calc_Fu(i,1,u,v,dt,dx,dy,dtau,Nx,Ny)
		! 	Fu(i,2) = calc_Fu(i,Ny,u,v,dt,dx,dy,dtau,Nx,Ny)
		! end do
		
		difference = 100 !大きな値にセット
		Res = 100 !同上
		cyc = 0
		
		!収束するまで繰り返し
		! do while(Res>1.e-17)
		! 	cyc = cyc + 1
		do cyc = 1, 200
			Prev(:) = z(:)
			
			!zの計算
			! call MGCYC(l,l,z,Au,Av,Az,b,nu1,nu2,Y,Nx,Ny,Nx/2,Ny/2,Res)
			! call MGCYC(l,l,z,Au,Az,b,nu1,nu2,Nx,Nx/2,Res,z(1:Nx),cyc,times)
			call smooth(z,Au,Az,b,Nx)
			Fu_tmp = calc_Fu(0,u,dt,dx,dtau,Nx)
			Fv_tmp = calc_Fv(0,u,v,dt,dx,dtau,Nx)
			call boundary(z,Nx,dt,dx,gamma,f,g,Fu_tmp,Fv_tmp,u)
			
			tmp(:) = reshape(Prev(:) - z(:),(/(Nx+2)/))
			difference = dot_product(tmp,tmp)
			
			call calc_res(z,Au,Az,b,Nx,Res)
			
			! if(times==4 .and. cyc<51) then
			! 	write(30,*) Res
			! end if
			! write(*,*) 'cyc = ', cyc, Res, difference
			
		end do
		! call boundary(z,Nx,Ny)
		
		write(*,*) 'times = ', times, sum(z)
		u_b(:) = u(:)
		v_b(:) = v(:)
		call calc_u(u,v_b,z,f,gamma,dt,dx,dtau,g,times,Nx)
		call calc_v(u_b,v,f,gamma,dt,dx,dtau,Nx)
		call calc_gamma(u,v,z,h,gamma,g,Cz,Nx)
		
		!時間計測終わり
		call system_clock(time_end_c)
		! print *,time_begin_c,time_end_c, CountPerSec,CountMax
		! write(*,*) 'nt = ', times, real(time_end_c - time_begin_c)/CountPerSec,"sec"
		
		!格子中心での値を記録
		if(mod(times,10)==0) then
			! do j = 1, Ny
			! 	do i = 1, Nx
			! 		write(10,*) z_frac(u(i-1:i,j))
			! 		write(11,*) z_frac(v(i,j-1:j))
			! 	end do
			! end do
			write(12,*) Az(1:Nx)
		endif
	end do

	stop
contains

	subroutine initialize(u,v,z,gamma,h,dx,Nx)
		implicit none

		integer, intent(in) :: Nx
		real(8), intent(in) :: dx
		real(8), intent(out) :: u(0:Nx), v(0:Nx+1), z(0:Nx+1), gamma(0:Nx+1), h(0:Nx+1)

		integer :: i, j

		u(:) = 0.5d0
		v(:) = 0.d0
		z(:) = 0.d0
		h(:) = 120.d0
		z(:) = 0.d0
		! do j = 1, Ny
		! 	do i = 1, Nx
				! if(i>Nx-5 .and. i<Nx-1 .and. j>Ny/2-2 .and. j<Ny/2+2) then
				! 	z(i,j) = 5.d0
				! end if
				! z(i,j) = 10*exp(-((i*dx-6.d6)**2+(j*dy-3.d6)**2)/2.d0/16.d4**2) !!Gaussian
				! z(i,j) = 10*exp(-((i-Nx/2)**2+(j-Ny/2)**2)/2.d0/2.d0**2) !!Gaussian
				! z(i,j) = 1*exp(-((i-Nx/4)**2+(j-Ny/2)**2)/2.d0/16.d0**2) + 1*exp(-((i-Nx/4*3)**2+(j-Ny/2)**2)/2.d0/16.d0**2)
				! h(i,j) = 1.d3 - 990.d0*(Nx-i)/Nx
		! 	end do
		! end do
		! h(:,:) = 1.d3 ![m]
		gamma(:) = 0.d0
		call boundary_u(u,Nx,1)
		
	end subroutine initialize

	recursive subroutine MGCYC(l,k,z,Au,Az,b,nu1,nu2,Nx,Nxc,Res,origin_zf,cyc,times)
		implicit none

		integer, intent(in) :: l, k, nu1, nu2, Nx, Nxc
		real(8), intent(in) :: Au(0:Nx), Az(1:Nx), b(1:Nx)
		real(8), intent(inout) :: z(0:Nx+1)
		real(8), intent(out) :: Res
		real(8), intent(in) :: origin_zf(1:Nx)
		integer, intent(in) :: cyc, times

		integer :: nt, i, ntmax
		real(8) :: df(1:Nx), dc(1:Nxc), wf(0:Nx+1), wc(0:Nxc+1) !defects and errors etc.
		real(8) :: Auc(0:Nxc), Azc(1:Nxc)
		real(8) :: origin_zc(1:Nxc)


		!Presmoothing
		do nt = 1, nu1
			call smooth(z,Au,Az,b,Nx)
		end do

		!Coarse grid correction
		!Compute the defect
		do i = 1, Nx
			df(i) = b(i) + Au(i-1)*z(i-1) + Au(i)*z(i+1) - Az(i)*z(i)
		end do

		!Restrice the defect
		dc(:) = 0.d0
		call Prolongation(Au,Az,Auc,Azc,Nxc)
		call Prolongation_defect(df,dc,Nxc)
		call Prolongation_defect(origin_zf,origin_zc,Nxc)

		!Compute an approximate solution v of the defect equation on k-1
		wc(:) = 0.d0
		if(k==1) then
			do nt = 1, 1
				call smooth(wc,Auc,Azc,dc,Nxc)
			end do
		else
			call MGCYC(l,k-1,wc,Auc,Azc,dc,nu1,nu2,Nxc,Nxc/2,Res,origin_zc,cyc,times)
		end if

		!Interpolate the correction
		wf(:) = 0.d0
		call Interpolation_defect(wc,wf,Nx)

		!Compute the corrected approximation on k
		z(:) = z(:) + wf(:)

		!Postsmoothing
		do nt = 1, nu2
			call smooth(z,Au,Az,b,Nx)
		end do

		if(mod(times,10)==0 .and. cyc==20) then
		! if(cyc==350) then
			select case(k)
			case(1)
				write(20,*) origin_zf(:) + z(1:Nx)
			case(2)
				write(21,*) origin_zf(:) + z(1:Nx)
			case(3)
				write(22,*) origin_zf(:) + z(1:Nx)
			case(4)
				write(23,*) origin_zf(:) + z(1:Nx)
			case(5)
				write(24,*) origin_zf(:) + z(1:Nx)
			case(6)
				write(25,*) origin_zf(:) + z(1:Nx)
			case(7)
				write(26,*) origin_zf(:) + z(1:Nx)
			case(8)
				write(27,*) origin_zf(:) + z(1:Nx)
			case default
			end select
		end if

	end subroutine MGCYC

	subroutine smooth(z,Au,Az,b,Nx)
		implicit none
		
		integer, intent(in) :: Nx
		real(8), intent(in) :: Au(0:Nx), Az(1:Nx), b(1:Nx)
		real(8), intent(inout) :: z(0:Nx+1)

		integer :: i

		do i = 1, Nx
		! do i = Nx, 1, -1
			z(i) = (b(i)+Au(i-1)*z(i-1)+Au(i)*z(i+1))/Az(i)
		end do

	end subroutine smooth

end program main
