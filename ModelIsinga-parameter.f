      program ModelIsinga
       implicit none
       integer d, i, j, k
       integer L, mcs, a
       parameter(L=100, mcs=230000)
       integer S(L,L), ni(L), pi(L)
       real Hamiltonian, ran1, magnetisation, average
       real susceptibility, kumulantBindera, thermalcapacity
       real R, dE
       real w, v, T, m, mag, X, U, C, H
       real y1, x1, Tc
       real magn(mcs/100), Energy(mcs/100)
*       parameter(kB = 1.38e-23)
       d=-1
       T=1.8
       open(11,file='av-magnL=100.txt')
       open(12,file='PodatnoscL=100.txt')
       open(13,file='kumBindL=100.txt')
       open(14,file='poj.cieplnaL=100.txt')
       open(15,file='HamiltonianL=100.txt')
       call neighbour(L,ni,pi)
       call cleanarray(mcs/100,magn,Energy)
       call randomset(L,S)
       
       write(*,*) 'T','  ','X(T)','  ','U(T)','   ','C(T)','    ','E(T)'
7      continue
        a = 1
        do k=1, mcs
          do i=1, L
             do j=1, L
              dE=2*S(i,j)*(S(ni(i),j)+S(pi(i),j)+S(i,ni(j))+S(i,pi(j)))
                v = exp(-dE/T)
                w = min(1.0,v)
                R = ran1(d)
                if(R<=w) S(i,j) = -S(i,j)
             enddo
          enddo
          if(k>=30000.and.mod(k,100)==0) then
*                magn(a) = abs(magnetisation(L,S))
                Energy(a) = Hamiltonian(L,S,ni,pi)
                a = a + 1
          endif
        enddo
        a = a - 1
         m = average(a,magn,1)
         X = susceptibility(L,a,S,T,magn)
         U = kumulantBindera(L,a,S,magn)
         C = thermalcapacity(L,a,S,T,Energy)
         H = average(a,Energy,1)/L**2
         write(11,*) T, m
         write(12,*) T, X
         if(T>=2.0.and.T<=2.5) write(13,*) T, U
         write(14,*) T, C
         write(15,*) T, H
         write(*,*) T, X, U, C, H
         T = T + 0.01
       if(T<=3.0) goto 7
       
       write(*,*) 'DONE!'
       PAUSE
      end
      
      subroutine randomset(k,A)
       integer k, A(k,k)
       real ran1, p
       integer d,i,j
       d = -1
       i = 1
       do i=1, k
         do j=1, k
            p = ran1(d)
            if(p<0.5) then
                A(i,j) = -1
            else
                A(i,j) = 1
            endif
         enddo
       enddo
       end
       
      subroutine neighbour(k,ni,pi)
       integer k, ni(k), pi(k), i
       do i=1, k
         ni(i) = i+1
         pi(i) = i-1
       enddo
       ni(k)=1
       pi(1)=k
       end
       
       subroutine cleanarray(n,p,q)
       integer n
       real p(n), q(n)
       do i=1, n
          p(i) = 0
          q(i) = 0
        enddo
       end
      
      function average(k,X,n)
       integer k, i, n
       real X(k)
       real average
       average = 0
       do i=1, k
        average = average + X(i)**n
       enddo
       average = average/(k+0.0)
       return
       end
       
      function Hamiltonian(k,A,ni,pi)
       integer k, i, j
       integer A(k,k), ni(k), pi(k)
       real H, Hamiltonian
       H=0
       do i=1, k
         do j=1,k
           H = H - A(i,j)*(A(ni(i),j)+A(pi(i),j)+A(i,ni(j))+A(i,pi(j)))
         enddo
       enddo
       Hamiltonian = H/2
       return
       end

      function magnetisation(k,A)
       integer k, i, j, A(k,k), n
       real m, magnetisation
       n = k**2
       m = 0
       do i=1, k
         do j=1, k
            m = m + A(i,j)
         enddo
       enddo
       magnetisation = m/(n+0.0)
       return
       end
      
      function susceptibility(k,m,A,T,tab)
       integer k, i, j, A(k,k), n, m
       real T, tab(m), susceptibility
       real average
       n = k**2
       susceptibility = n/T*(average(m,tab,2) - average(m,tab,1)**2)
       return
       end
       
      function thermalcapacity(k,m,A,T,Etab)
       integer k, A(k,k), m
       real T, Etab(m)
       real thermalcapacity, average
       thermalcapacity=(average(m,Etab,2)-average(m,Etab,1)**2)/(T**2)
       thermalcapacity = thermalcapacity/k**2
       return
       end
       
      function kumulantBindera(k,m,A,mag)
       integer k, A(k,k), n
       real mag(m), kumulantBindera, average
       n = k**2
       kumulantBindera = 1-average(m,mag,4)/(3*average(m,mag,2)**2)
       return
       end
       
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM

      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
