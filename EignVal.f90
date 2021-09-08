        PROGRAM TEZHENG_Jacobi
! ��Jacobi�������������������ֵ����������
!   ������ƽ����ת����U���϶Ծ���A���������Ʊ任��A��Ϊ�ԽǾ���,
!   �Ӷ����A������ֵ����������

        IMPLICIT NONE
        INTEGER NMAX,N,I,J,P,Q,change,k
        PARAMETER(NMAX=1000)
        REAL*8 A(NMAX,NMAX),AMAX,TEMP,ZEMP,COO,SII,CO,SI,APP,AQQ,APQ,API,AQI
        REAL*8 R(NMAX,NMAX),RIP,RIQ,temp2,temp3
        REAL*8 tempArr(NMAX)
        CHARACTER NAME*12,NAMEO*12,CHR*1
! ���ļ��ж���ʵ�Գƾ���A
        WRITE(*,*)'����ʵ�Գƾ���ά��n(n<1000):'
        READ(*,*) N
        WRITE(*,*)'��������ļ�:'
        READ(*,*) NAME
        OPEN(6,FILE=NAME)
        DO I=1,N
         READ(6,*) (A(I,J),J=1,I)
         DO J=1,I
         A(J,I)=A(I,J)
         ENDDO
        ENDDO
        CLOSE(6)
! R�����������任����U,�����ȳ�ʼ��,����λ����
        DO I=1,N
        DO J=1,N
        R(I,J)=0
        ENDDO
        R(I,I)=1
        ENDDO
! �ھ���A�ķ����Խ���Ԫ����,�ҳ���ģ����Ԫ��Apq
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


! �������Խ���Ԫ�ػ�Ϊ0,��С�ڸ�������ʱ,�������ֵ����������
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
        !WRITE(*,*) 'A������ֵΪ:'
        WRITE(*,*) (A(I,I),I=1,N)
        !WRITE(*,*) 'A����������Ϊ:'
        !WRITE(*,*) '  X1      X2     X3     ...:'
!        DO I=1,N
!        WRITE(*,*)(R(I,J),J=1,N)
!        ENDDO

        WRITE(*,*) '�Ƿ񽫽�������ļ�(Y/N)?'
        READ(*,*) CHR
        IF(CHR.EQ.'Y'.OR.CHR.EQ.'y') THEN
!        WRITE(*,*) '�����ļ���(С��12�ַ�):'
!        READ(*,*) NAMEO
        OPEN(8,FILE='EigE.txt')
        open(9,file='EigV.txt')
!        WRITE(8,*) 'A������ֵΪ:'
        WRITE(8,*) (A(I,I),I=1,N)
!        WRITE(8,*) 'A����������Ϊ:'
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
!��ʼ׼������ƽ����ת����U
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
! ����ƽ����ת����U
        DO I=1,N
        RIP=R(I,P)*CO+R(I,Q)*SI
        RIQ=-R(I,P)*SI+R(I,Q)*CO
        R(I,P)=RIP
        R(I,Q)=RIQ
        ENDDO
! ��A���б任       
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
