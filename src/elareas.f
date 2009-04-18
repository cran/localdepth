      SUBROUTINE elareas(X, nt, nc, nrx, dduse, darea)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension darea(nt), x(nrx, nc), xsimplex(nc+1,nc)
      dimension isimplex(nc+1)

      external ldarea

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

CC Evaluate the volume of a given simplex
CC          write(*,*) isimplex
          nsimp = nsimp+1
          do 50 ii=1,(nc+1)
            is = isimplex(ii)
            do 60 jj=1,nc
              xsimplex(ii,jj) = x(is,jj)
 60         continue
 50       continue
          call elarea(xsimplex, nc, dduse, darea(nsimp))
CC End of Evaluate the volume of a given simplex
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

      SUBROUTINE elareasa(X, nt, nc, nrx, dduse, darea)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      dimension darea(nt), x(nrx, nc), xsimplex(nc+1,nc)
      dimension isimplex(nrx)

      external ldarea
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
        call elarea(xsimplex, nc, dduse, darea(nsimp))
 20   continue
CC      write(*,*) nsimp

      call rndend()
      return
      end

      SUBROUTINE elarea(X, nc, dduse, darea)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

      dimension x(nc+1, nc), xx(nc+1, nc+1)
      dimension ipvt(nc+1), ztemp(nc+1), dwork(nc+1), ddeth(2)

      external dgeco
      external dgedi
      external dgamma

      dc = nc

      do 10 i=1,(nc+1)
        do 20 j=2,(nc+1)
          xx(i,j) = x(i,j-1)
 20     continue
 10   continue
      do 30 k=1,(nc+1)
        xx(k,1) = duno
 30   continue

      call dgeco(xx,(nc+1),(nc+1),ipvt,rcond,ztemp)
      call dgedi(xx,(nc+1),(nc+1),ipvt,ddeth,dwork,10)	

      darea = (dduse**(dc/ddue))*dabs(ddeth(1)*(10.0d00**ddeth(2)))
     &        *(dpi**(dc/ddue))/dgamma(dc/ddue+duno)

      return
      end
