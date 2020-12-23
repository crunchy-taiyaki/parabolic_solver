program parabolic_equation
    use precision
    use input_functions
    use parabolic_equation_solver
    implicit none
    integer(mp) :: N, M !size of x and t vectors, respectively
    real(mp) :: sigma
    real(mp), allocatable :: x(:), t(:)
    real(mp), allocatable :: solution(:,:), reference_solution(:,:), residuals(:,:)
    integer(mp) :: i,k ! counters
    N = 100
    M = 100
    sigma = 0.5_mp
    allocate(x(0:N))
    allocate(t(0:M))
    allocate(solution(0:N,0:M))
    allocate(reference_solution(0:N,0:M))
    allocate(residuals(0:N,0:M))
    
    call parabolic_solver(N,M,phi,f,Lh_x_u,alpha,alpha1,alpha2,beta,beta1,beta2,sigma,x,t,solution)
 
    do k=0,M
        do i=0,N
            reference_solution(i,k) = u_ref(x(i),t(k))
            residuals(i,k) = solution(i,k) - reference_solution(i,k)
 !           write(*,*) i, solution(i,k)- u_ref(x(i),t(k))
        enddo
    enddo
    write(*,*) maxval(abs(residuals))
    end program