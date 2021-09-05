program main
    implicit none
    real*8,parameter :: dp=0.01,dt=0.1,Pc=5.0
    integer n,i,Nmax
    parameter(Nmax=1000)
    complex B(Nmax,Nmax)
    n=Pc/dp
    call OperatorB(dp,dt,5000.0*dt,n,B)
    do i=1,n+1
        print *, B(i,i)
    end do
    pause
end program
    
subroutine OperatorB(dp,dt,t,n,B)
    !Êä³öB=exp(-i*dt*HI(t+dt/2))
    implicit none
    integer i,Nmax,n
    real*8 E0,omg,tau,pi
    parameter (Nmax=1000)
    parameter (E0=5.337e-4,omg=2.85e-2,tau=1102,pi=4.0*datan(1.0d0))
    complex HI(Nmax,Nmax),B(Nmax,Nmax)
    real*8 dp,dt,t
    real*8 At
    complex,parameter :: cj=(0.,1.)
    At=E0*sin((t+dt/2)*pi/tau)**2*cos(omg*(t+dt/2))/omg !At=A(t+dt/2)
    HI=0
    B=0
    do i=1,n+1
        HI(i,i)=i*dp*At+At**2/2 !HI(i,i)=pi*A(t)+A(t)^2/2
    end do
    do i=1,n+1
        B(i,i)=exp(-cj*dt*HI(i,i))
    end do
end subroutine