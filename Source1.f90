program main
implicit none

real*8,parameter :: dp=0.01,Pc=6.0
call GetH0(dp,Pc)

pause
end program

!求出算符H0并写入文件"H0_3.txt"
subroutine GetH0(dp,Pc)
implicit none

real*8 :: Vp,p,pi,dp,Pc
real*8,parameter :: h=0.01,a=-200.0,b=200.0
!,dp=0.01,Pc=5
!p=5
!call RSQ(h,a,b,p,Vp)
!print "(a,f10.5)", "The result is Vp=",Vp
integer :: n,i,j
real*8,allocatable :: H0(:,:)
n=Pc/dp
allocate(H0(-n/2:n/2,-n/2:n/2))
pi=4.d0*datan(1.d0)
do i=-n/2,n/2
    print *, i
    do j=-n/2,n/2
	    if (j==i) then
		    call RSQ(h,a,b,abs(i-j)*dp,Vp)
		    H0(i,j)=0.5*(i*dp)**2+(dp*Vp)/(2*pi)
		else if (j>i) then
		    call RSQ(h,a,b,abs(i-j)*dp,Vp)
		    H0(i,j)=(dp*Vp)/(2*pi)
		else
		    H0(i,j)=H0(j,i)
		end if
	end do
end do

open(unit=10,file='H0_3.txt')
do i=-n/2,n/2
    write(10,"(601f)") H0(i,:)
end do
close(10)

print *, "Success!"
end subroutine

!Simpson法求数值积分，精度h,积分限[a,b],V(p)=r
subroutine RSQ(h,a,b,p,r)
implicit none

real*8 :: h,a,b,p,r
integer :: n,i
complex,parameter :: cj=(0.,1.)
real*8,allocatable :: x(:)

n=(b-a)/(2*h)

allocate(x(0:2*n))
do i=0,2*n
x(i)=a+(b-a)/(2*n)*i
end do

r=0
do i=1,n
r=r-(2*exp(cj*p*x(2*i-1))/sqrt(1+(x(2*i-1))**2)+exp(cj*p*x(2*i))/sqrt(1+(x(2*i))**2))
!r=r-(2.0d0/sqrt(1+(x(2*i-1))**2)+1.0d0/sqrt(1+(x(2*i))**2))
!r=r+(2*sin(x(2*i-1))+sin(x(2*i)))
end do
r=r*2
r=h/3*(-exp(cj*p*a)/sqrt(1+(a)**2)+exp(cj*p*b)/sqrt(1+(b)**2)+r)
!r=h/3*(-1.0d0/sqrt(1+a**2)+1.0d0/sqrt(1+b**2)+r)
!r=h/3*(sin(a)-sin(b)+r)
end subroutine


