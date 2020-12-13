module input_functions
    use precision
    implicit none
    contains
    
    real function p(x)
    real(mp), intent(in) :: x
    p = x+1.0_mp
    end function p
    
    real function b(x,t)
    real(mp), intent(in) :: x,t
    b = 0.0_mp
    end function b
    
    real function c(x,t)
    real(mp), intent(in) :: x,t
    c = -1.0_mp
    end function c
    
    real function alpha1(t)
    real(mp), intent(in) :: t
    alpha1 = 1.0_mp
    end function alpha1
    
    real function alpha2(t)
    real(mp), intent(in) :: t
    alpha2 = 0.0_mp
    end function alpha2
    
    real function beta1(t)
    real(mp), intent(in) :: t
    beta1 = 1.0_mp
    end function beta1
    
    real function beta2(t)
    real(mp), intent(in) :: t
    beta2 = 1.0_mp
    end function beta2 
    
    real function u_ref(x,t)
    real(mp), intent(in) :: x,t
    u_ref = x**3 + t**3
    end function u_ref
    
    real function alpha(t,dx)
    real(mp), intent(in) :: t, dx
    real(mp) :: x0
    x0 = 0.0_mp
    !dx = 1.0d-3
    alpha = alpha1(t)*u_ref(x0,t) + alpha2(t)*(u_ref(x0+dx,t)-u_ref(x0-dx,t))/(2.0_mp*dx)
    end function alpha
    
    real function beta(t,dx)
    real(mp), intent(in) :: t,dx
    real(mp) :: x0
    x0 = 1.0_mp
    !dx = 1.0d-3
    beta = beta1(t)*u_ref(x0,t) + beta2(t)*(u_ref(x0+dx,t)-u_ref(x0-dx,t))/(2.0_mp*dx) 
    end function beta
    
    real function L(x,t)
    real(mp), intent(in) :: x,t
    real(mp) :: dx
    dx = 1.0d-2
    L = p(x+dx/2.0_mp)*(u_ref(x+dx,t)-u_ref(x,t))/(dx*dx) - p(x-dx/2.0_mp)*(u_ref(x,t)-u_ref(x-dx,t))/(dx*dx) + &
        & b(x,t)*(u_ref(x+dx,t)-u_ref(x-dx,t))/(2.0_mp*dx) + c(x,t)*u_ref(x,t)
    end function L
    
    real function Lh_x_u(x,t,h,uk_prev,uk,uk_next)
    real(mp), intent(in) :: x,t,h
    real(mp), intent(in) :: uk_prev,uk,uk_next
    Lh_x_u = p(x+h/2.0_mp)*(uk_next-uk)/(h*h) - p(x-h/2.0_mp)*(uk-uk_prev)/(h*h) + &
             &  b(x,t)*(uk_next-uk_prev)/(2.0_mp*h) + c(x,t)*uk
    end function Lh_x_u
    
    subroutine ABC_coeff(x,t,h,tau,sigma,A_values,B_values,C_values)
    real(mp), intent(in) :: x(0:),t(0:),h,tau,sigma
    integer(mp) :: N, M
    real(mp), intent(out) :: A_values(0:size(x)-1,0:size(t)-1), B_values(0:size(x)-1,0:size(t)-1), C_values(0:size(x)-1,0:size(t)-1)
    integer :: i, k !counters
    N = size(x)-1
    M = size(t)-1
    do k=0,M
        do i=1,N-1
            A_values(i,k) = -sigma*(p(x(i)-h/2.0_mp)/(h*h) + b(x(i),t(k))/(2.0_mp*h))
            B_values(i,k) = sigma*(p(x(i)+h/2.0_mp)/(h*h) + p(x(i)-h/2.0_mp)/(h*h) - c(x(i),t(k))) + 1/tau
            C_values(i,k) = sigma*(p(x(i)+h/2.0_mp)/(h*h) + b(x(i),t(k))/(2.0_mp*h))
        enddo
        A_values(0,k) = 0.0_mp
        A_values(N,k) = -(beta2(t(k))/h)
        B_values(0,k) = -(alpha1(t(k)) + alpha2(t(k))/h)
        B_values(N,k) = -(beta1(t(k))+beta2(t(k)/h))
        C_values(0,k) = -alpha2(t(k))/h
        C_values(N,k) = 0.0_mp


    enddo
    end subroutine ABC_coeff
    
    real function phi(x)
    real(mp), intent(in) :: x
    phi = u_ref(x,0.0_mp)
    end function phi
    
    real function f(x,t,dt)
    real(mp), intent(in) :: x,t,dt
    f = (u_ref(x,t+dt)-u_ref(x,t-dt))/(2.0_mp*dt) - L(x,t)
    end function f 
    
end module input_functions