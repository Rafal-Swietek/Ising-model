      program Sth
       implicit none
       real T, a, b, c, d
       integer i
       open(5,file='poj.cieplnaL=5.txt')
       open(10,file='poj.cieplnaL=10.txt')
       open(30,file='poj.cieplnaL=30.txt')
       open(100,file='poj.cieplnaL=100.txt')
       open(15,file='poj.cieplna-L=5.txt')
       open(110,file='poj.cieplna-L=10.txt')
       open(130,file='poj.cieplna-L=30.txt')
       open(1100,file='poj.cieplna-L=100.txt')
       write(*,*)'--------------------------------------------'
       DO i=1, 1000
         read(5,*) T, a
         read(10,*) T, b
         read(30,*) T, c
         read(100,*) T, d
         a = a*25
         b = b*100
         c = c*9000
         d = d*10000
         write(*,*) T, ' ', a , ' ', b , ' ' , c, ' ' , d
         write(15,*) T, a
         write(110,*) T, b
         write(130,*) T, c
         write(1100,*) T, d
       enddo
       close(5)
       close(10)
       close(30)
       close(100)
      end
