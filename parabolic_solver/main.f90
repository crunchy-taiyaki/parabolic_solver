program parabolic_equation
    use precision
    use input_functions
    use parabolic_equation_solver
    implicit none
    integer(mp) :: N, M !size of x and t vectors, respectively
    real(mp) :: sigma
    real(mp), allocatable :: x(:), t(:)
    real(mp), allocatable :: solution(:,:), residuals(:,:)
    integer(mp) :: i,k ! counters
    N = 10
    M = 50
    sigma = 1.0_mp
    allocate(x(0:N))
    allocate(t(0:M))
    allocate(solution(0:N,0:M))
    allocate(residuals(0:N,0:M))
    
    !test beta
    !call coord_grid(x,t,N,M)
    !do k=0,M
    !    write(*,*) beta(t(k)), t(k)**3 + 4.0_mp
    !enddo
    
    !test L
    !call coord_grid(x,t,N,M)
    !do k=0,M
    !do i=0,N
    !    write(*,*)k,i, L(x(i),t(k)) - (9.0_mp*x(i)**2 + 6.0_mp*x(i) - u_ref(x(i),t(k)))
    !enddo
    !enddo
    
    !test f(x,t)
    !call coord_grid(x,t,N,M)
    !do k=0,M
    !do i=0,N
    !    write(*,*)k,i, f(x(i),t(k)) - (3.0_mp*t(k)**2 - (9.0_mp*x(i)**2 + 6.0_mp*x(i) - u_ref(x(i),t(k))))
    !enddo
    !enddo
    
    call parabolic_solver(N,M,phi,f,Lh_x_u,alpha,alpha1,alpha2,beta,beta1,beta2,sigma,x,t,solution)
    do k=0,M
        do i=0,N
            residuals(i,k) = solution(i,k) - u_ref(x(i),t(k))
            write(*,*)k,i,residuals(i,k) !x(i),solution(i,k),
        enddo
    enddo
    end program