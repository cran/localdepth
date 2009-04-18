      SUBROUTINE lddiams(X, nt, nc, nrx, diam)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension diam(nt), x(nrx, nc), xsimplex(nc+1,nc)
      dimension isimplex(nc+1)

      external lddiam

C Inizializza il vettore degli indici degli spigoli dei simplessi
      do 10 i=1,(nc+1)
        isimplex(i) = i
 10   continue
      nsimp = 0
      i = nc+1
      do 20 while (isimplex(1).le.(nrx-nc))
        icont = 1
        do 30 while (i.lt.(nc+1).and.icont.eq.1)
          if (isimplex(i).lt.(nrx-nc+i)) then
            i = i + 1
          else 
            icont = 0
          endif
 30     continue
        if (isimplex(i).le.(nrx-nc-1+i)) then

CC Evaluate the diameter of a given simplex
CC          write(*,*) isimplex
          nsimp = nsimp+1
          do 50 ii=1,(nc+1)
            is = isimplex(ii)
            do 60 jj=1,nc
              xsimplex(ii,jj) = x(is,jj)
 60         continue
 50       continue
          call lddiam(xsimplex, nc, diam(nsimp))
CC End of Evaluate the diamter of a given simplex
          isimplex(i) = isimplex(i)+1
        else
          isimplex(i-1) = isimplex(i-1)+1
          j = i
          do 40 while (j.le.(nc+1))
            isimplex(j) = isimplex(j-1)+1
            j = j+1
 40       continue
          i = i-1
        endif      
 20   continue
CC      write(*,*) nsimp
      return
      end

      SUBROUTINE lddiamsa(X, nt, nc, nrx, diam)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension diam(nt), x(nrx, nc), xsimplex(nc+1,nc)
      dimension isimplex(nrx)

      external lddiam
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      do 20 nsimp=1,nt
        do 10 i=1,nrx
          isimplex(i)=i
 10     continue
        do 50 ii=1,(nc+1)
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          do 60 jj=1,nc
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue
        call lddiam(xsimplex, nc, diam(nsimp))
 20   continue
CC      write(*,*) nsimp

      call rndend()
      return
      end

      SUBROUTINE lddiam(X, nc, diam)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)
      parameter(ddue=2.0d00)

      dimension x(nc+1, nc)

      diam = dzero
      do 10 i=1,nc
        do 20 j=(i+1),(nc+1)
          temp = dzero
          do 30 k=1,nc
            temp = temp+(x(i,k)-x(j,k))**ddue 
 30       continue
          temp = dsqrt(temp)
          if (diam.le.temp) then
            diam = temp
          endif
 20     continue
 10   continue

      return
      end
