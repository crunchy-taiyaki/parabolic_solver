module parabolic_equation_solver
    use precision
    use input_functions, only : ABC_coeff, u_ref
    implicit none
    contains
    
    subroutine coord_grid(x,t,N,M,h,tau)
    integer(mp), intent(in) :: N, M
    real(mp), intent(out) :: x(0:N), t(0:M)
    real(mp), intent(out) :: h,tau
    integer(mp) :: i,k
    h = 1.0_mp/N
    tau = 0.1_mp/M
    forall(i=0:N) x(i) = i*h
    forall(k=0:M) t(k) = k*tau    
    end subroutine coord_grid
    
    subroutine tridiagonal_matrix_algorithm(A,B,C,G,y)
    real(mp), intent(in) :: A(0:),B(0:),C(0:),G(0:)
    real(mp), intent(out) :: y(0:)
    real(mp),allocatable :: t(:),s(:)
    integer(mp) :: N
    integer(mp) :: i !counter
    N = size(A)-1
    allocate(t(0:N),s(0:N))
    s(0) = C(0)/B(0)
    t(0) = -G(0)/B(0)
    do i=1,N
        s(i) = C(i)/(B(i)-A(i)*s(i-1))
        t(i) = (A(i)*t(i-1)-G(i))/(B(i) - A(i)*s(i-1))
    enddo
    y(N) = t(N)
    do i=N-1,0,-1
        y(i) = s(i)*y(i+1) + t(i)
    enddo
    deallocate(t,s)
    end subroutine tridiagonal_matrix_algorithm
    
    subroutine parabolic_solver(N,M,phi,f,Lh_x_u,alpha,alpha1,alpha2,beta,beta1,beta2,sigma,x,t,u)
    interface
        function phi(x)
        integer, parameter :: mp=8      
        real(mp), intent(in) :: x
        real(mp) :: phi
        end function phi
        
        function f(x,t,dt,dx)
        integer, parameter :: mp=8
        real(mp), intent(in) :: x,t,dt,dx
        real(mp) :: f
        end function f
        
        function Lh_x_u(x,t,h,uk_prev,uk,uk_next)
        integer, parameter :: mp=8
        real(mp), intent(in) :: x,t,h
        real(mp), intent(in) :: uk_prev,uk,uk_next
        real(mp) :: Lh_x_u
        end function Lh_x_u
        
        function alpha(t,dx)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t,dx
        real(mp) :: alpha
        end function alpha
        
        function alpha1(t)
        integer, parameter :: mp=8
        real(mp) :: alpha1
        real(mp), intent(in) :: t
        end function alpha1
        
        function alpha2(t)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t
        real(mp) :: alpha2
        end function alpha2
        
        function beta(t,dx)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t,dx
        real(mp) :: beta
        end function beta
        
        function beta1(t)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t
        real(mp) :: beta1
        end function beta1
        
        function beta2(t)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t
        real(mp) :: beta2
        end function beta2        
    end interface
    
    integer(mp), intent(in) :: N, M
    real(mp), intent(in) :: sigma
    real(mp),intent(out) :: x(0:N), t(0:M)
    real(mp), intent(out) :: u(0:N,0:M)
    real(mp) :: h, tau
    real(mp) :: t_underlined(1:M)
    real(mp), allocatable :: A(:,:),B(:,:),C(:,:),G(:,:)
    integer(mp) :: i, k !counters
    allocate(A(0:N,0:M),B(0:N,0:M),C(0:N,0:M),G(0:N,0:M))
    call coord_grid(x,t,N,M,h,tau)
    do i=0,N !k=0
        u(i,0) = phi(x(i))
    enddo   
    if (sigma.eq.0.0_mp) then
        t_underlined = t(0:M-1)
    elseif (sigma.eq.0.5_mp) then
        t_underlined = t(1:M) - tau/2.0_mp
    elseif (sigma.eq.1.0_mp) then
         t_underlined = t(1:M)
    endif
    call ABC_coeff(x,t,h,tau,sigma,A,B,C)
        do k=1,M
            do i=1,N-1        
                G(i,k) = -u(i,k-1)/tau - (1.0_mp-sigma)*Lh_x_u(x(i),t(k-1),h,u(i-1,k-1),u(i,k-1),u(i+1,k-1)) - f(x(i),t_underlined(k),tau,h)
            enddo
            G(0,k) = alpha(t(k),h)
            G(N,k) = beta(t(k),h)
            call tridiagonal_matrix_algorithm(A(:,k),B(:,k),C(:,k),G(:,k),u(:,k))
        enddo
    deallocate(A,B,C,G)        
    end subroutine parabolic_solver
    
    end module parabolic_equation_solver