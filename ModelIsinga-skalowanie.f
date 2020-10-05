      program Skalowanie
       implicit none
       real y1, x1, T, Tc, m10, m20, m40, m100
       integer L, i
       real Tempcritical
       open(10,file='av-magnL=10.txt')
       open(40,file='av-magnL=40.txt')
       open(100,file='av-magnL=100.txt')
       open(11,file='skaL=10.txt')
       open(41,file='skaL=40.txt')
       open(101,file='skaL=100.txt')
       Tc = Tempcritical(50)
       write(*,*) 'Tc = ',Tc
       pause
       write(*,*)'--------------------------------------------'
       DO i=1, 100
         read(10,*) T, m10
         write(*,*) T, m10
         L = 10
         y1 = log(m10) + 0.125*log(L+0.0)
         x1 = log(abs(1-T/Tc)*(L+0.0))
         write(11,*) x1, y1
       enddo
       close(11)
       close(10)
       write(*,*)'--------------------------------------------'
       DO i=1, 100
         read(40,*) T, m40
         write(*,*) T, m40
         L = 40
         y1 = log(m40) + 0.125*log(L+0.0)
         x1 = log(abs(1-T/Tc)*(L+0.0))
         write(41,*)x1, y1
       enddo
       close(40)
       close(41)
       write(*,*)'--------------------------------------------'
       DO i=1, 100
         read(100,*) T, m100
         write(*,*) T, m100
         L = 100
         y1 = log(m100) + 0.125*log(L+0.0)
         x1 = log(abs(1-T/Tc)*(L+0.0))
         write(101,*) x1, y1
       enddo
       close(100)
       close(101)
      end
      
      Function Tempcritical(k)
       integer i, k, error
       real U100(k), U10(k), U30(k)
       real Temp(k)
       real min1, min2, Tempcritical
       real a, b, T1, T2
       open(21,file='kumBindL=10.txt')
       open(22,file='kumBindL=30.txt')
       open(23,file='kumBindL=100.txt')
       i = 1
       DO i=1, k
         read(21,*) Temp(i), U10(i)
         read(22,*) x, U30(i)
         read(23,*) x, U100(i)
         write(*,*) Temp(i), U10(i), U30(i), U100(i)
       enddo
       close(21)
       close(22)
       close(23)
       min1 = abs(1 - U100(20)/U30(20))
       min2 = abs(1 - U30(20)/U10(20))
       do i=21,k
          a = abs(1 - U100(i)/U30(i))
          b = abs(1 - U30(i)/U10(i))
          if(a < min1) then
             min1 = a
             T1 = Temp(i)
          endif
          if(b < min2) then
             min2 = b
             T2 = Temp(i)
          endif
       enddo
       write(*,*) T1, T2
       Tempcritical = (T1 + T2)/2.0
       return
       end
      

