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
    real(mp),intent(out) :: y(0:)
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
    do i=N-1,0,-1
        y(i) = s(i)*y(i+1) + t(i)
    enddo
    deallocate(t,s)
    end subroutine tridiagonal_matrix_algorithm
    
    subroutine parabolic_solver(N,M,phi,f,Lh_x_u,alpha,alpha1,alpha2,beta,beta1,beta2,sigma,x,t,u)
    interface
        real function phi(x)
        integer, parameter :: mp=8      
        real(mp), intent(in) :: x
        end function phi
        
        real function f(x,t,dt)
        integer, parameter :: mp=8
        real(mp), intent(in) :: x,t,dt
        end function f
        
        real function Lh_x_u(x,t,h,uk_prev,uk,uk_next)
        integer, parameter :: mp=8
        real(mp), intent(in) :: x,t,h
        real(mp), intent(in) :: uk_prev,uk,uk_next
        end function Lh_x_u
        
        real function alpha(t,dx)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t,dx
        end function alpha
        
        real function alpha1(t)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t
        end function alpha1
        
        real function alpha2(t)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t
        end function alpha2
        
        real function beta(t,dx)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t,dx
        end function beta
        
        real function beta1(t)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t
        end function beta1
        
        real function beta2(t)
        integer, parameter :: mp=8
        real(mp), intent(in) :: t
        end function beta2        
    end interface
    
    integer(mp), intent(in) :: N, M
    real(mp), intent(in) :: sigma
    real(mp),intent(out) :: x(0:N), t(0:M)
    real(mp), intent(out) :: u(0:N,0:M)
    real(mp) :: h, tau
    real(mp) :: t_underlined(0:M)
    real(mp), allocatable :: A(:,:),B(:,:),C(:,:),G(:,:)
    integer(mp) :: i, k !counters
    allocate(A(0:N,0:M),B(0:N,0:M),C(0:N,0:M),G(0:N,0:M))
    call coord_grid(x,t,N,M,h,tau)
    do i=0,N !k=0
        u(i,0) = phi(x(i))
    enddo
    if (sigma.eq.0.0_mp) then
        do k=1,M
            do i=1,N-1
                u(i,k) = u(i,k-1) + tau*(Lh_x_u(x(i),t(k-1),h,u(i-1,k-1),u(i,k-1),u(i+1,k-1)) + f(x(i),t(k-1),tau))
                write(*,*)k,i, u(i,k)-u_ref(x(i),t(k))
            enddo
        u(0,k) = (alpha(t(k),h) + 2.0_mp*alpha2(t(k))*u(1,k)/h - alpha2(t(k))*u(2,k)/(2.0_mp*h))/&
                 & (alpha1(t(k)) + 3.0_mp*alpha2(t(k))/(2.0_mp*h))
        write(*,*)k,0, u(0,k)-u_ref(x(0),t(k))
        u(N,k) = (beta(t(k),h) + 2.0_mp*beta2(t(k))*u(N-1,k)/h - beta2(t(k))*u(N-2,k)/(2.0_mp*h))/&
                 & (beta1(t(k)) + 3.0_mp*beta2(t(k))/(2.0_mp*h))
        write(*,*)k,N, u(N,k)-u_ref(x(N),t(k))
        enddo
    elseif ((sigma == 0.5_mp).or.(sigma == 1.0_mp))  then
        if (sigma.eq.1.0_mp) then
        t_underlined = t
        else
        t_underlined = t - tau/2.0_mp
        endif
    call ABC_coeff(x,t,h,tau,sigma,A,B,C)
        do k=1,M
            do i=1,N-1        
                G(i,k) = -u(i,k-1)/tau - (1-sigma)*Lh_x_u(x(i),t(k-1),h,u(i-1,k-1),u(i,k-1),u(i+1,k-1)) - f(x(i),t_underlined(k),tau)
            enddo
            G(0,k) = alpha(t(k),h)
            G(N,k) = beta(t(k),h)
            call tridiagonal_matrix_algorithm(A(:,k),B(:,k),C(:,k),G(:,k),u(:,k))
        enddo
    deallocate(A,B,C,G)
    else
        write(*,*) sigma, 'this sigma is undefined, please chose one of these: 0.0, 0.5, 1.0'
    endif
        
    end subroutine parabolic_solver
    
    end module parabolic_equation_solver