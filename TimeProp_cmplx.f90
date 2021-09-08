program main
    implicit none
    real*8,parameter:: dp=0.01,dt=0.1,Pc=6.0
    integer n
    n=Pc/dp
    call TimeProp(dt,dp,n)
    pause
end program
    
subroutine TimeProp(dt,dp,n)
    implicit none
    integer,parameter :: Nmax=1000
    real*8,parameter:: E0=5.337e-2,omg=5.7e-2,tau=1102,pi=4.0*datan(1.0d0)
    real*8 dt,dp,t,At
    double complex B_prime(Nmax,Nmax)
    !Psi即含时波函数,A-算符A;B-算符B;D-算符p+A(t)
    double complex,allocatable :: dv(:),Psi(:),Psi1(:),A(:,:),B(:,:),D(:,:),D2(:,:)
    integer n,m,i,j !m-时间轴的格数,n-p轴格数
    character PsiIn*12,NameA*15
    m=floor(tau/dt)
    allocate(dv(0:m))
    allocate(Psi(1:n+1))
    allocate(Psi1(1:n+1))
    allocate(A(1:n+1,1:n+1))
    allocate(B(1:n+1,1:n+1))
    allocate(D(1:n+1,1:n+1))
    allocate(D2(1:n+1,1:n+1))
    
    !读取初始波函数
 !   print "(a)","请指定初始波函数的文件名："
 !   read(*,*) PsiIn
    open(10,file='phi1.dat')
    do i=1,n+1
    read(10,*) Psi(i)
    enddo
    close(10)
    
    !读取算符A=exp(=i*H0*0.5dt)
    !print "(a)" "请指定算符A的输入文件名："
    !NameA="OperatorA.txt"
    NameA="OA.txt"
    open(11,file='OA.txt')
    do i=1,n+1
    do j=1,n+1
    read(11,*) A(i,j)
    end do
    enddo
    close(11)
    !算符D
    do i=1,n+1
        D(i,i)=(i-n/2-1)*dp    !i*dp
    end do
    B=0
    !含时演化
    do j=0,m
       t=j*dt 
       At=E0*sin(t*pi/tau)**2*cos(omg*t)/omg !At=A(t)
       do i=1,n+1
       D2(i,i)=D(i,i)+At
       enddo

       dv(j)=dot_product(Psi,matmul(D2,Psi))  !*dp
       call OperatorB(dp,dt,t,n,B_prime)
       do i=1,n+1
           B(i,i)=B_prime(i,i)
       end do
       !Split-Operator演化更新Psi
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
    !输出B=exp(-i*dt*HI(t+dt/2))
    implicit none
    integer i,Nmax,n
    real*8 E0,omg,tau,pi
    parameter (Nmax=1000)
    parameter (E0=5.337e-2,omg=5.7e-2,tau=1102,pi=4.0*datan(1.0d0))
    double complex HI(Nmax,Nmax),B(Nmax,Nmax)
    real*8 dp,dt,t
    real*8 At
    double complex,parameter :: cj=(0.,1.)
    At=E0*sin((t+dt/2)*pi/tau)**2*cos(omg*(t+dt/2))/omg !At=A(t+dt/2)
    HI=0
    B=0
    do i=1,n+1
        HI(i,i)=(i-n/2-1)*dp*At+At**2/2 !HI(i,i)=pi*A(t)+A(t)^2/2
    end do
    do i=1,n+1
        B(i,i)=exp(-cj*dt*HI(i,i))
    end do
end subroutine
