program main
    implicit none
    real*8,parameter:: dp=0.01,dt=0.1,Pc=5.0
    integer n
    n=Pc/dp
    call TimeProp(dt,dp,n)
    pause
end program
    
subroutine TimeProp(dt,dp,n)
    implicit none
    integer,parameter :: Nmax=1000
    real*8,parameter:: E0=5.337e-4,omg=2.85e-2,tau=1102,pi=4.0*datan(1.0d0)
    real*8 dt,dp,t,At
    complex B_prime(Nmax,Nmax)
    !Psi����ʱ������,A-���A;B-���B;D-���p+A(t)
    complex,allocatable :: dv(:),Psi(:),Psi1(:),A(:,:),B(:,:),D(:,:)
    integer n,m,i,j !m-ʱ����ĸ���,n-p�����
    character PsiIn*12,NameA*15
    m=floor(tau/dt)
    allocate(dv(0:m))
    allocate(Psi(1:n+1))
    allocate(Psi1(1:n+1))
    allocate(A(1:n+1,1:n+1))
    allocate(B(1:n+1,1:n+1))
    allocate(D(1:n+1,1:n+1))
    
    !��ȡ��ʼ������
    print "(a)","��ָ����ʼ���������ļ�����"
    read(*,*) PsiIn
    open(10,file=PsiIn)
    read(10,*) (Psi(i),i=1,n+1)
    close(10)
    
    !��ȡ���A=exp(=i*H0*0.5dt)
    !print "(a)" "��ָ�����A�������ļ�����"
    !NameA="OperatorA.txt"
    NameA="A_cmplx.txt"
    open(11,file=NameA)
    do i=1,n+1
    read(11,*) (A(i,:))
    end do
    close(11)
    !���D
    do i=1,n+1
        D(i,i)=i*dp
    end do
    B=0
    !��ʱ�ݻ�
    do j=0,m
       t=j*dt 
       At=E0*sin(t*pi/tau)**2*cos(omg*t)/omg !At=A(t)
       D=D+At
       dv(j)=dot_product(Psi,matmul(D,Psi))*dp
       call OperatorB(dp,dt,t,n,B_prime)
       do i=1,n+1
           B(i,i)=B_prime(i,i)
       end do
       !Split-Operator�ݻ�����Psi
       Psi1=matmul(A,Psi)
       Psi=matmul(B,Psi1)
       Psi=matmul(A,Psi)
    end do
    !open(12,file="Dipole.txt")
    open(12,file="Dipole_real.txt")
    do i=0,m
    !print "(f)", dv(i)
    !write(12,*) (dv(i))
        write(12,*) (real(dv(i)))
    end do
end subroutine
    
subroutine OperatorB(dp,dt,t,n,B)
    !���B=exp(-i*dt*HI(t+dt/2))
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