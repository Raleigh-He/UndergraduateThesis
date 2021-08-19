program main
    implicit none
    real*8,parameter:: dt=0.1,tau=1102
    integer m
    m=tau/dt
    call HHG(dt,m)
end program
    
subroutine HHG(dt,m)
    implicit none
    complex cj,sum
    parameter(cj=(0.,1.))
    real*8 dt
    real*8,parameter:: omg=2.85e-2
    real*8,parameter:: domg=0.1*omg,OmgMax=35*omg
    integer m,l,i,j !m=tau/dt,时间轴的格数
    real*8,allocatable :: dv(:),Fomg(:)
    l=OmgMax/domg !l=350,Omega轴的格数
    allocate(dv(0:m),Fomg(0:l))
    !allocate(Fomg(0:l))
    
    !读取dv
    open(10,file="Dipole.txt")
    do i=0,m
    read(10,*) (dv(i))
    end do
    close(10)
    
    do i=0,l
        sum=0
        do j=0,m
            sum=sum+dv(j)*exp(-cj*j*dt*i*domg)*dt
        end do
        Fomg(i)=sum*conjg(sum)
    end do
    
    !输出F(Omg)
    open(11,file="HHGOut.txt")
    !write(11,*) (Fomg(i),i=0,l)
    write(11,"(351e)") (Fomg)
    close(11)
    print *, "Finished!"
    pause
end subroutine