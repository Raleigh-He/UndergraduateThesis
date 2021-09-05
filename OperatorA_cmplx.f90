program main
    implicit none
    real*8 dp,Pc,dt
    dp=0.01
    Pc=5.0
    dt=0.1
    call OperatorA(dp,Pc,dt)
    pause
end program

subroutine OperatorA(dp,Pc,dt)
implicit none

complex cj
parameter(cj=(0.,1.))!������λ
real*8 dp,Pc,dt
integer i,j,k
integer n
character InName*12,OutName*15
real*8,allocatable :: EigVal(:),EigVec(:,:)
complex,allocatable :: A(:,:)
n=Pc/dp!��n=5.0/0.01=500
allocate(A(1:n+1,1:n+1))
allocate(EigVal(1:n+1))
allocate(EigVec(1:n+1,1:n+1))
write(*,*) '��ָ�������ļ�����'
read(*,*) InName
!��ȡH0������ֵ��������
open(10,file=InName)
read(10,*) (EigVal(i),i=1,n+1)
do i=1,n+1
    read(10,*) (EigVec(i,j),j=1,n+1)
end do
close(10)
!�������A=exp(-i*H0*dt)
do i=1,n+1
    do j=1,i
        A(i,j)=0
        do k=1,n+1
            A(i,j)=A(i,j)+exp(-cj*0.5*dt*EigVal(k))*EigVec(i,k)*EigVec(j,k)
        end do
    end do
    do j=1,i
        A(j,i)=A(i,j)
    end do
end do
write(*,*) '��ָ������ļ�����'
read(*,*) OutName
open(11,file=Outname)
do i=1,n+1
WRITE(11,*) A(i,:)
end do
close(11)
!print "(2e10.5)",A(1,2),A(2,1)
!write(*,*) '���A�ľ�����ʽ��'
!do i=1,n+1
!WRITE(*,*)(A(I,J),J=1,N+1)
!print "(/)"
!end do
    !do j=1,i
     !   EigVec(j,i)=EigVec(i,j)
    !end do

 !WRITE(*,*) 'A������ֵΪ:'
 !       WRITE(*,*) (EigVal(I),I=1,N+1)
 !       WRITE(*,*) 'A����������Ϊ:'
        !WRITE(*,*) '  X1      X2     X3     ...:'
        !DO I=1,N+1
        !WRITE(*,*)(EigVec(I,J),J=1,N+1)
 !       WRITE(*,*)(EigVec(2,J),J=1,N+1)
        !ENDDO

end subroutine