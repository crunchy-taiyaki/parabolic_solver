module input_functions
    use precision
    implicit none
    contains
    
    function p(x)
    real(mp), intent(in) :: x
    real(mp) :: p
    p = x+1.0_mp
    end function p
    
    function b(x,t)
    real(mp), intent(in) :: x,t
    real(mp) :: b
    b = 0.0_mp
    end function b
    
    function c(x,t)
    real(mp), intent(in) :: x,t
    real(mp) :: c
    c = -1.0_mp
    end function c
    
    function alpha1(t)
    real(mp), intent(in) :: t
    real(mp) :: alpha1
    alpha1 = 1.0_mp
    end function alpha1
    
    function alpha2(t)
    real(mp), intent(in) :: t
    real(mp) :: alpha2
    alpha2 = 0.0_mp
    end function alpha2
    
    function beta1(t)
    real(mp), intent(in) :: t
    real(mp) :: beta1
    beta1 = 1.0_mp
    end function beta1
    
    function beta2(t)
    real(mp), intent(in) :: t
    real(mp) :: beta2
    beta2 = 1.0_mp
    end function beta2 
    
    function u_ref(x,t)
    real(mp), intent(in) :: x,t
    real(mp) :: u_ref
    u_ref = x**3 + t**3
    end function u_ref
    
    function alpha(t,h)
    real(mp), intent(in) :: t, h
    real(mp) :: alpha
    real(mp) :: x0
    x0 = 0.0_mp
    alpha = alpha1(t)*u_ref(x0,t) + alpha2(t)*(u_ref(x0+h,t)-u_ref(x0,t))/h
    end function alpha
    
    function beta(t,h)
    real(mp), intent(in) :: t,h
    real(mp) :: beta
    real(mp) :: x1
    x1 = 1.0_mp
    beta = beta1(t)*u_ref(x1,t) + beta2(t)*(u_ref(x1,t)-u_ref(x1-h,t))/h
    end function beta
    
    function L(x,t,h)
    real(mp), intent(in) :: x,t,h
    real(mp) :: L

    L = p(x+h/2.0_mp)*(u_ref(x+h,t)-u_ref(x,t))/h**2 - p(x-h/2.0_mp)*(u_ref(x,t)-u_ref(x-h,t))/h**2 + &
        & b(x,t)*(u_ref(x+h,t)-u_ref(x-h,t))/2.0_mp/h + c(x,t)*u_ref(x,t)
    end function L
    
    function Lh_x_u(x,t,h,uk_prev,uk,uk_next)
    real(mp), intent(in) :: x,t,h
    real(mp), intent(in) :: uk_prev,uk,uk_next
    real(mp) Lh_x_u
    Lh_x_u = p(x+h/2.0_mp)*(uk_next-uk)/h**2 - p(x-h/2.0_mp)*(uk-uk_prev)/h**2 + &
             &  b(x,t)*(uk_next-uk_prev)/2.0_mp/h + c(x,t)*uk
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
            A_values(i,k) = sigma*( p(x(i)-h/2.0_mp)/h**2 - b(x(i),t(k))/2.0_mp/h )
            B_values(i,k) = sigma*( p(x(i)+h/2.0_mp)/h**2 + p(x(i)-h/2.0_mp)/h**2 - c(x(i),t(k)) ) + 1.0_mp/tau
            C_values(i,k) = sigma*( p(x(i)+h/2.0_mp)/h**2 + b(x(i),t(k))/2.0_mp/h )
        enddo
        A_values(0,k) = 0.0_mp 
        A_values(N,k) = -beta2(t(k)) / h 
        B_values(0,k) = -( alpha1(t(k)) + alpha2(t(k)) / h  )
        B_values(N,k) = -( beta1(t(k)) + beta2(t(k)) / h  )
        C_values(0,k) = -alpha2(t(k)) / h
        C_values(N,k) = 0.0_mp
    enddo
    end subroutine ABC_coeff
    
    function phi(x)
    real(mp), intent(in) :: x
    real(mp) :: phi
    phi = u_ref(x,0.0_mp)
    end function phi
    
    function f(x,t,tau,h)
    real(mp), intent(in) :: x,t,tau,h
    real(mp) :: f
    f = (u_ref(x,t+tau)-u_ref(x,t-tau))/2.0_mp/tau - L(x,t,h)
    end function f 
    
end module input_functions