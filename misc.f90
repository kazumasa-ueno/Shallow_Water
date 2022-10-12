module misc
  implicit none
  real(8) :: t_s

contains

  subroutine update(u,v,h)
    real(8),intent(inout) :: u(:,:,:),v(:,:,:),h(:,:,:)

    u(:,:,-1) = u(:,:,0)
		u(:,:,0) = u(:,:,1)
		v(:,:,-1) = v(:,:,0)
		v(:,:,0) = v(:,:,1)
		h(:,:,-1) = h(:,:,0)
		h(:,:,0) = h(:,:,1)
  end subroutine update

  subroutine start_timer()
    real(8) :: omp_get_wtime
    t_s = omp_get_wtime()
  end subroutine start_timer

  double precision function get_elapsed_time()
    real(8) :: omp_get_wtime
    get_elapsed_time = omp_get_wtime() - t_s
  end function get_elapsed_time

end module misc
