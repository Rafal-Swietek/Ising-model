      program ModelIsinga
       implicit none
       integer d, i, j, k, a
       integer L, mcs, N
       parameter(L=100, mcs=20000)
       integer S(L,L), ni(L), pi(L)
       real ran1, magnetisation
       real R, dE
       real w, v, T, m, X, U, C
       d=-1
       a=1
       T=3.3
*       open(11,file='magn.txt')
       call neighbour(L,ni,pi)
       call randomset(L,S)
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
*           m = magnetisation(L,S)
*           if(k>=30000) write(11,*) k, m
           write (*,*) k
           call showmatrix(L,S)
       enddo
       call showmatrix(L,S)
       write(*,*) 'Done!'
       read(*,*)
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
       
      function magnetisation(k,A)
       integer k, i, j, A(k,k), n
       real m, magnetisation
       n = k**2
       M = 0
       do i=1, k
         do j=1, k
            M = M + A(i,j)
         enddo
       enddo
       magnetisation = m/(n+0.0)
       return
       end
       
      subroutine showmatrix(n, B)
       integer n, B(n,n)
       !open(88,file='matrixL=200,T=3.3.txt')
       open(88,file='data.txt')
       do j=1, n
          write(88,*)(B(i,j), i=1, n)
       enddo
*       close(88)
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
