program parabolic_equation
    use precision
    use input_functions
    use parabolic_equation_solver
    implicit none
    integer(mp) :: N, M, N2, M2 !size of x and t vectors, respectively
    real(mp) :: sigma
    real(mp), allocatable :: x(:), t(:), x2(:), t2(:)
    real(mp), allocatable :: solution(:,:),solution2(:,:), reference_solution(:,:),reference_solution2(:,:), residuals(:,:),residuals2(:,:)
    integer(mp) :: i,k ! counters
    N = 20
    M = 20
    N2 = N/2
    M2 = M
    sigma = 0.5_mp
    allocate(x(0:N),x2(0:N2))
    allocate(t(0:M),t2(0:M2))
    allocate(solution(0:N,0:M),solution2(0:N2,0:M2))
    allocate(reference_solution(0:N,0:M),reference_solution2(0:N2,0:M2))
    allocate(residuals(0:N,0:M),residuals2(0:N2,0:M2))
    
    call parabolic_solver(N,M,phi,f,Lh_x_u,alpha,alpha1,alpha2,beta,beta1,beta2,sigma,x,t,solution)
    call parabolic_solver(N2,M2,phi,f,Lh_x_u,alpha,alpha1,alpha2,beta,beta1,beta2,sigma,x2,t2,solution2)
    
    write(*,*)'U(x,t)'
    write(*,'("t \ x     ",6("    &    ",f10.3))') x(0), x(1), x(2), x(3), x(4), x(5)
    do k=0,M
        write(*,'(f10.3,6("    &    ",f10.5))') t(k), solution(0,k),&
                                                            &solution(1,k),&
                                                            &solution(2,k),&
                                                            &solution(3,k),&
                                                            &solution(4,k),&
                                                            &solution(5,k)
    enddo
    write(*,*)'U*(x,t)'
    write(*,'("t \ x     ",6("    &    ",f10.3))') x(0), x(1), x(2), x(3), x(4), x(5)
    do k=0,M
        write(*,'(f10.3,6("    &    ",f10.5))') t(k), u_ref(x(0),t(k)),&
                                                            &u_ref(x(1),t(k)),&
                                                            &u_ref(x(2),t(k)),&
                                                            &u_ref(x(3),t(k)),&
                                                            &u_ref(x(4),t(k)),&
                                                            &u_ref(x(5),t(k))
    enddo
 
    do k=0,M
        do i=0,N
            reference_solution(i,k) = u_ref(x(i),t(k))
            residuals(i,k) = solution(i,k) - reference_solution(i,k)
        enddo
    enddo
    
    do k=0,M
        do i=0,N2
            residuals2(i,k) = solution(i*2,k) - solution2(i,k)
        enddo
    enddo
    write(*,'("h, tau",2(f10.3))') 1.0/N, 0.1/M
    write(*,'(e10.5)') maxval(abs(residuals)), maxval(abs(residuals2))
    end program