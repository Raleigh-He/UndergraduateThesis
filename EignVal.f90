        PROGRAM TEZHENG_Jacobi
! 用Jacobi方法求正交矩阵的特征值与特征向量
!   就是用平面旋转矩阵U不断对矩阵A作正交相似变换把A化为对角矩阵,
!   从而求出A的特征值与特征向量

        IMPLICIT NONE
        INTEGER NMAX,N,I,J,P,Q,change,k
        PARAMETER(NMAX=1000)
        REAL*8 A(NMAX,NMAX),AMAX,TEMP,ZEMP,COO,SII,CO,SI,APP,AQQ,APQ,API,AQI
        REAL*8 R(NMAX,NMAX),RIP,RIQ,temp2,temp3
        REAL*8 tempArr(NMAX)
        CHARACTER NAME*12,NAMEO*12,CHR*1
! 从文件中读入实对称矩阵A
        WRITE(*,*)'输入实对称矩阵维数n(n<1000):'
        READ(*,*) N
        WRITE(*,*)'输入矩阵文件:'
        READ(*,*) NAME
        OPEN(6,FILE=NAME)
        DO I=1,N
         READ(6,*) (A(I,J),J=1,I)
         DO J=1,I
         A(J,I)=A(I,J)
         ENDDO
        ENDDO
        CLOSE(6)
! R矩阵存放正交变换矩阵U,在这先初始化,即单位矩阵
        DO I=1,N
        DO J=1,N
        R(I,J)=0
        ENDDO
        R(I,I)=1
        ENDDO
! 在矩阵A的非主对角线元素中,找出按模最大的元素Apq
100        AMAX=dABS(A(2,1))
        P=2
        Q=1
        DO I=2,N
         DO J=1,I-1
         IF(dABS(A(I,J)).GT.AMAX) THEN
         AMAX=dABS(A(I,J))
         P=I
         Q=J
         ENDIF
         ENDDO
        ENDDO


! 当非主对角线元素化为0,即小于给定精度时,输出特征值与特征向量
        IF(AMAX.LE.1.0E-14) THEN
            do i=1,N-1
                change=0
                do j=1,N-1
                    if(A(j,j).GT.A(j+1,j+1)) then
                        change=1
                        temp2=A(j,j)
                        A(j,j)=A(j+1,j+1)
                        A(j+1,j+1)=temp2
                        do k=1,N
                           tempArr(k)=R(k,j)
                           R(k,j)=R(k,j+1)
                           R(k,j+1)=tempArr(k) 
                        end do
                    end if
                end do
                if(change.EQ.0) then
                    exit
                end if
            end do
        !WRITE(*,*) 'A的特征值为:'
        WRITE(*,*) (A(I,I),I=1,N)
        !WRITE(*,*) 'A的特征向量为:'
        !WRITE(*,*) '  X1      X2     X3     ...:'
!        DO I=1,N
!        WRITE(*,*)(R(I,J),J=1,N)
!        ENDDO

        WRITE(*,*) '是否将结果存入文件(Y/N)?'
        READ(*,*) CHR
        IF(CHR.EQ.'Y'.OR.CHR.EQ.'y') THEN
!        WRITE(*,*) '输入文件名(小于12字符):'
!        READ(*,*) NAMEO
        OPEN(8,FILE='EigE.txt')
        open(9,file='EigV.txt')
!        WRITE(8,*) 'A的特征值为:'
        WRITE(8,*) (A(I,I),I=1,N)
!        WRITE(8,*) 'A的特征向量为:'
!        WRITE(8,*) '  X1      X2     X3     ...:'
        DO I=1,N
        do j=1,N
        WRITE(9,*)R(I,J)
        ENDDO
        enddo
        CLOSE(8)
        open(2,file='phi1.dat')
        temp3=0.d0
        do i=1,N
        write(2,*)R(i,1)
        temp3=temp3+R(i,1)*R(i,2)
        enddo
        write(*,*)temp3
        ENDIF

        STOP
        ENDIF
!开始准备计算平面旋转矩阵U
        TEMP=2*A(P,Q)/(A(P,P)-A(Q,Q)+1.0e-30)
        ZEMP=(A(P,P)-A(Q,Q))/(2*A(P,Q))
       
        IF(ABS(TEMP).LT.1.0) THEN
        COO=(1+TEMP**2)**(-0.5)
        SII=TEMP*(1+TEMP**2)**(-0.5)
        ELSE
        COO=ABS(ZEMP)*(1+ZEMP**2)**(-0.5)
        SII=SIGN(1.0,ZEMP)*(1+ZEMP**2)**(-0.5)
        ENDIF
       
        CO=SQRT(0.5*(1+COO))
        SI=SII/(2*CO)
! 计算平面旋转矩阵U
        DO I=1,N
        RIP=R(I,P)*CO+R(I,Q)*SI
        RIQ=-R(I,P)*SI+R(I,Q)*CO
        R(I,P)=RIP
        R(I,Q)=RIQ
        ENDDO
! 对A进行变换       
        APP=A(P,P)*CO**2+A(Q,Q)*SI**2+2*A(P,Q)*CO*SI
        AQQ=A(P,P)*SI**2+A(Q,Q)*CO**2-2*A(P,Q)*CO*SI
        APQ=0.5*(A(Q,Q)-A(P,P))*SII+A(P,Q)*COO
        A(P,P)=APP
        A(Q,Q)=AQQ
        A(P,Q)=APQ
        A(Q,P)=A(P,Q)
       
        DO I=1,N
        IF(I.EQ.P.OR.I.EQ.Q) THEN
        ELSE
        API=A(P,I)*CO+A(Q,I)*SI
        AQI=-A(P,I)*SI+A(Q,I)*CO
        A(P,I)=API
        A(Q,I)=AQI
        A(I,P)=A(P,I)
        A(I,Q)=A(Q,I)
        ENDIF
        ENDDO
       
      GOTO 100  
        END
