      SUBROUTINE ldareashs(X, nt, nc, nrx, darea)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension darea(nt), x(nrx, nc), xsimplex(nc,nc)
      dimension isimplex(nc)

      external ldareahs

C Inizializza il vettore degli indici degli spigoli dei simplessi
      do 10 i=1,nc
        isimplex(i) = i
 10   continue
      nsimp = 0
      i = nc
      do 20 while (isimplex(1).le.(nrx-nc+1))
        icont = 1
        do 30 while (i.lt.nc.and.icont.eq.1)
          if (isimplex(i).lt.(nrx-nc+1+i)) then
            i = i + 1
          else 
            icont = 0
          endif
 30     continue
        if (isimplex(i).le.(nrx-nc+i)) then

CC Evaluate the volume of a given simplex
CC          write(*,*) isimplex
          nsimp = nsimp+1
          do 50 ii=1,nc
            is = isimplex(ii)
            do 60 jj=1,nc
              xsimplex(ii,jj) = x(is,jj)
 60         continue
 50       continue
          call ldareahs(xsimplex, nc, darea(nsimp))
CC End of Evaluate the volume of a given simplex
          isimplex(i) = isimplex(i)+1
        else
          isimplex(i-1) = isimplex(i-1)+1
          j = i
          do 40 while (j.le.nc)
            isimplex(j) = isimplex(j-1)+1
            j = j+1
 40       continue
          i = i-1
        endif      
 20   continue
      return
      end

      SUBROUTINE areashsa(X, nt, nc, nrx, diam)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension diam(nt), x(nrx, nc), xsimplex(nc,nc)
      dimension isimplex(nrx)

      external ldareahs
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      do 20 nsimp=1,nt
        do 10 i=1,nrx
          isimplex(i)=i
 10     continue
        do 50 ii=1,nc
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          do 60 jj=1,nc
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue
        call ldareahs(xsimplex, nc, diam(nsimp))
 20   continue
      call rndend()
      return
      end

      SUBROUTINE ldareahs(X, nc, darea)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)

      dimension x(nc, nc)
      dimension darchi(nc)
      do 5 i=1,nc
        darchi(i) = dzero
 5    continue
      dperi = dzero
      ii = 0
      do 10 i=1,(nc-1)
        do 20 j=(i+1),(nc)
          temp = dzero
          do 30 k=1,nc
            temp = temp+x(i,k)*x(j,k) 
 30       continue
          ii = ii + 1
          darchi(ii) = dacos(temp)
          if (darchi(ii).lt.dzero) then
            darchi(ii) = dzero
          endif
          dperi= dperi + darchi(ii)
 20     continue
 10   continue
      dtemp = dtan(dperi/4.0d00)
     &  *dtan((dperi/2.0d00-darchi(1))/2.0d00)
     &  *dtan((dperi/2.0d00-darchi(2))/2.0d00)
     &  *dtan((dperi/2.0d00-darchi(3))/2.0d00)
      if (dtemp.lt.dzero) then
        dtemp = dzero
      endif
      darea = 4.0d00*datan(dsqrt(dtemp))
      return
      end
