!********************************************
!	赤道ケルビン波の伝播
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
	
	integer, parameter :: l = 6								!グリッドの深さ
	! integer, parameter :: Nx = 160, Ny = 80		!グリッド数
	integer, parameter :: Nx = 512, Ny = 128
	! integer, parameter :: Nx = 16, Ny = 8
	integer, parameter :: ntmax = 4000				!時間ステップ
	integer, parameter :: nu1 = 2, nu2 = 1		!マルチグリッドサイクル内のsmooth回数
	real(8), parameter :: g = 9.81d0 	!修正重力定数
	! real(8), parameter :: g = 9.81d0 	!重力定数
	real(8), parameter :: Cz = 80.d0					!Chezy 摩擦係数
	real(8), parameter :: pi = 4*atan(1.d0)		!円周率
	real(8), parameter :: f0 = 4*pi/86400			!コリオリパラメータf0
	real(8), parameter :: X = 30000d3, Y = 7500d3	!領域サイズ
	! real(8), parameter :: dt = 60.d0*4				!時間間隔
	real(8), parameter :: dt = 216.d0				!時間間隔
	real(8), parameter :: dtau = dt/10.d0			!移流計算用小時間間隔

	real(8) :: f(0:Ny+1) 	!コリオリパラメータ
	real(8) :: h(0:Nx+1,0:Ny+1) !基準面からの水深
	real(8) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), u_b(0:Nx,0:Ny+1), v_b(0:Nx+1,0:Ny) !u_b, v_bは移流計算実行用の一時格納配列
	! real(8) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny) !係数
	real(8) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny), q(1:Nx,1:Ny) !係数
	real(8) :: dx, dy	!格子間隔
	real(8) :: Res, difference	!Resは残差のl2ノルム、differenceは前の時間との残差	
	integer :: times, cyc !時間ループ用と収束までの繰り返し用
	integer :: ios !ファイル開く用
	integer :: i, j !空間ループ用

	!for debug
	real(8) :: Prev(0:Nx+1,0:Ny+1), tmp((Nx+2)*(Ny+2)) !前の値を格納しておくための配列

	! open(unit=10, file="./output/u.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/u.txt"
	! open(unit=11, file="./output/v.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/v.txt"
	open(unit=12, file="./output/z.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z.txt"
	! open(unit=20, file="./output/z1.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	! open(unit=21, file="./output/z2.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/z2.txt"
	! open(unit=22, file="./output/z3.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/z3.txt"
	! open(unit=23, file="./output/z4.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	! open(unit=24, file="./output/z5.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/z2.txt"
	! open(unit=25, file="./output/z6.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/z3.txt"
	! open(unit=26, file="./output/z7.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	! open(unit=27, file="./output/z8.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	! open(unit=30, file="./output/res5.txt", iostat=ios, status="replace", action="write")
	! if ( ios /= 0 ) stop "Error opening file ./output/res.txt"
	
	

	!!コリオリfと格子間隔dx,dyの計算
	call calc_f(f,Ny,f0,Y)
	call when_l(l,X,Y,Nx,Ny,dx,dy)

	!u,v,z,gamma,hの初期化
	call initialize(u,v,z,gamma,h,dx,dy,Nx,Ny)

	!メインのループ
	do times = 0, ntmax-1
		!時間計測スタート
		call system_clock(time_begin_c, CountPerSec, CountMax)

		!係数計算
		call calc_Au(z,gamma,h,g,dt,dx,Nx,Ny,Au)
		call calc_Av(z,gamma,h,g,dt,dy,Nx,Ny,Av)
		call calc_Az(Au,Av,Nx,Ny,Az)
		call calc_b(u,v,z,gamma,h,f,dt,dx,dy,dtau,Nx,Ny,times,b,q)
		write(*,*) Au(256,128), Av(256,128), Az(256,128), b(256,128), q(256,128)


		difference = 100 !大きな値にセット
		Res = 100 !同上
		cyc = 0

		!収束するまで繰り返し
		! do while(Res>1.e-18 .and. cyc<300)
		! 	cyc = cyc + 1
		do cyc = 1, 20
			Prev(:,:) = z(:,:)

			!zの計算
			call MGCYC(l,l,z,Au,Av,Az,b,nu1,nu2,Y,Nx,Ny,Nx/2,Ny/2,Res,&
			& u(1:Nx,1:Ny),v(1:Nx,1:Ny),z(1:Nx,1:Ny),gamma(1:Nx,1:Ny),h(1:Nx,1:Ny),cyc,times,dx,dy)
			! call smooth(z,Au,Av,Az,b,Nx,Ny)
			call boundary_z(z,Nx,Ny)


			tmp(:) = reshape(Prev(:,:) - z(:,:),(/(Nx+2)*(Ny+2)/))
			difference = dot_product(tmp,tmp)

			call calc_res(z,Au,Av,Az,b,Nx,Ny,Res)
			
			! if(times==4 .and. cyc<51) then
			! 	write(30,*) Res
			! end if
			write(*,*) 'cyc = ', cyc, Res, difference

		end do
		! call boundary(z,Nx,Ny)

		write(*,*) 'times = ', times, sum(z)
		u_b(:,:) = u(:,:)
		v_b(:,:) = v(:,:)
		call calc_u(u,v_b,z,f,gamma,dt,dx,dy,dtau,g,Nx,Ny)
		call calc_v(u_b,v,z,f,gamma,dt,dx,dy,dtau,g,Nx,Ny)
		call calc_gamma(u,v,z,h,gamma,g,Cz,Nx,Ny)

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
			! write(12,*) z(1:Nx,1:Ny)*0.0006d0
			write(12,*) z(1:Nx,1:Ny)
		endif
	end do

	stop
contains

	subroutine initialize(u,v,z,gamma,h,dx,dy,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: dx, dy
		real(8), intent(out) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1),h(0:Nx+1,0:Ny+1)

		integer :: i, j

		u(:,:) = 0.d0
		v(:,:) = 0.d0
		h(:,:) = 120.d0
		z(:,:) = 0.d0
		do j = 1, Ny
			do i = 1, Nx
		! 		! if(i>Nx-5 .and. i<Nx-1 .and. j>Ny/2-2 .and. j<Ny/2+2) then
		! 		! 	z(i,j) = 5.d0
		! 		! end if
				! z(i,j) = 10*exp(-((i*dx-6.d6)**2+(j*dy-3.d6)**2)/2.d0/16.d4**2) !!Gaussian
				! z(i,j) = 10*exp(-((i-Nx/2)**2+(j-Ny/2)**2)/2.d0/2.d0**2) !!Gaussian
				z(i,j) = 1*exp(-((i-Nx/4)**2+(j-Ny/2)**2)/2.d0/8.d0**2) + 1*exp(-((i-Nx/4*3)**2+(j-Ny/2)**2)/2.d0/8.d0**2) !!Gaussian
		! 		! h(i,j) = 1.d3 - 990.d0*(Nx-i)/Nx
			end do
		end do
		! h(:,:) = 1.d3 ![m]
		gamma(:,:) = 0.d0
		
	end subroutine initialize

	recursive subroutine MGCYC(l,k,z,Au,Av,Az,b,nu1,nu2,Y,Nx,Ny,Nxc,Nyc,Res, &
	& uf,vf,zf,gammaf,hf,cyc,times,dx,dy)
		implicit none

		integer, intent(in) :: l, k, nu1, nu2, Nx, Ny, Nxc, Nyc
		real(8), intent(in) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny), Y, dx, dy
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)
		real(8), intent(out) :: Res
		real(8), intent(in) :: uf(1:Nx,1:Ny), vf(1:Nx,1:Ny), zf(1:Nx,1:Ny), gammaf(1:Nx,1:Ny), hf(1:Nx,1:Ny)
		integer, intent(in) :: cyc, times

		integer :: nt, i, j, ntmax
		real(8) :: df(1:Nx,1:Ny), dc(1:Nxc,1:Nyc), wf(0:Nx+1,0:Ny+1), wc(0:Nxc+1,0:Nyc+1) !defects and errors etc.
		real(8) :: Auc(0:Nxc,1:Nyc), Avc(1:Nxc,0:Nyc), Azc(1:Nxc,1:Nyc), bc(1:Nxc,1:Nyc)
		real(8) :: uc(1:Nxc,1:Nyc), vc(1:Nxc,1:Nyc), zc(0:Nxc+1,0:Nyc+1),&
		& gammac(0:Nxc+1,0:Nyc+1), hc(0:Nxc+1,0:Nyc+1), cor(0:Nxc+1,0:Nyc+1)


		!Presmoothing
		do nt = 1, nu1
			call smooth(z,Au,Av,Az,b,Nx,Ny)
		end do

		!Coarse grid correction
		!Compute the defect
		do j = 1, Ny
			do i = 1, Nx
				df(i,j) = b(i,j) + Au(i-1,j)*z(i-1,j) + Au(i,j)*z(i+1,j) + Av(i,j-1)*z(i,j-1) + Av(i,j)*z(i,j+1) - Az(i,j)*z(i,j)
			end do
		end do

		!Restrice the defect
		dc(:,:) = 0.d0
		call Prolongation(zf,gammaf,hf,zc,gammac,hc,Nxc,Nyc)
		call Prolongation_defect(df,dc,Nxc,Nyc)
		call calc_Au(zc,gammac,hc,g,dt,dx*2,Nxc,Nyc,Auc)
		call calc_Av(zc,gammac,hc,g,dt,dy*2,Nxc,Nyc,Avc)
		call calc_Az(Auc,Avc,Nxc,Nyc,Azc)

		!compute the right-hand side
		do j = 1, Nyc
			do i = 1, Nxc
				bc(i,j) = dc(i,j) + Auc(i-1,j)*zc(i-1,j) + Auc(i,j)*zc(i+1,j) + Avc(i,j-1)*zc(i,j-1) + Avc(i,j)*zc(i,j+1) - Azc(i,j)*zc(i,j)
			end do
		end do

		!Compute an approximate solution v of the defect equation on k-1
		wc(:,:) = 0.d0
		if(k==1) then
			do nt = 1, 1
				! call smooth(wc,Auc,Avc,Azc,bc,Nxc,Nyc)
				call smooth(wc,Auc,Avc,Azc,dc,Nxc,Nyc)
			end do
		else
			! call MGCYC(l,k-1,wc,Auc,Avc,Azc,bc,nu1,nu2,Y,Nxc,Nyc,Nxc/2,Nyc/2,Res,uc,vc,zc,gammac,hc,cyc,times)
			call MGCYC(l,k-1,wc,Auc,Avc,Azc,dc,nu1,nu2,Y,Nxc,Nyc,Nxc/2,Nyc/2,Res,uc,vc,zc,gammac,hc,cyc,times,dx*2,dy*2)
		end if

		!Compute the correction
		! cor = wc - zc

		!Interpolate the correction
		wf(:,:) = 0.d0
		! call Interpolation_defect(cor,wf,Nx,Ny)
		call Interpolation_defect(wc,wf,Nx,Ny)

		!Compute the corrected approximation on k
		z(:,:) = z(:,:) + wf(:,:)

		!Postsmoothing
		do nt = 1, nu2
			call smooth(z,Au,Av,Az,b,Nx,Ny)
		end do

		! if(mod(times,10)==0 .and. cyc==20) then
		! ! if(cyc==350) then
		! 	select case(k)
		! 	case(1)
		! 		write(20,*) zf(:,:) + z(1:Nx,1:Ny)
		! 	case(2)
		! 		write(21,*) zf(:,:) + z(1:Nx,1:Ny)
		! 	case(3)
		! 		write(22,*) zf(:,:) + z(1:Nx,1:Ny)
		! 	case(4)
		! 		write(23,*) zf(:,:) + z(1:Nx,1:Ny)
		! 	case(5)
		! 		write(24,*) zf(:,:) + z(1:Nx,1:Ny)
		! 	case(6)
		! 		write(25,*) zf(:,:) + z(1:Nx,1:Ny)
		! 	case(7)
		! 		write(26,*) zf(:,:) + z(1:Nx,1:Ny)
		! 	case(8)
		! 		write(27,*) zf(:,:) + z(1:Nx,1:Ny)
		! 	case default
		! 	end select
		! end if

	end subroutine MGCYC

	subroutine smooth(z,Au,Av,Az,b,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny)
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

		integer :: i, j

		! do j = 1, Ny
		do j = Ny, 1, -1
			do i = 1, Nx
			! do i = Nx, 1, -1
				z(i,j) = (b(i,j)+Au(i-1,j)*z(i-1,j)+Au(i,j)*z(i+1,j)+Av(i,j-1)*z(i,j-1)+Av(i,j)*z(i,j+1))/Az(i,j)
			end do
		end do

	end subroutine smooth

end program main
