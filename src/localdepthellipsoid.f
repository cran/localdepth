C############################################################
C
C	Functions for the exact ellipsoid (local) depth
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: September, 14, 2008
C	Version: 0.1-1
C
C	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
C
C############################################################

      SUBROUTINE lde(X, Y, dtau, nc, nt, nrx, nry, nuse, dtol,
     & dduse, depth, depthlocal)
CC dduse e' la dimensione dell'elissoide al quadrato

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), depthlocal(nry), y(nry,nc), x(nrx, nc)
      dimension isimplex(nc+1), xsimplex(nc+1, nc), depthsim(nry)

      external ldei
      external elarea
      external eldiam

      do 5 i=1,nry
        depth(i) = dzero 
        depthlocal(i) = dzero 
 5    continue

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

CC Evaluate the depth for a given simplex
CC          write(*,*) isimplex
          nsimp = nsimp+1
          do 50 ii=1,(nc+1)
            is = isimplex(ii)
            do 60 jj=1,nc
              xsimplex(ii,jj) = x(is,jj)
 60         continue
 50       continue
CC calculate the dimension of the simplex
          if (nuse.eq.0) then
            call eldiam(xsimplex, nc, dduse, ddim)   
          else
            call elarea(xsimplex, nc, dduse, ddim)   
          endif

CC calculate the depth
          call ldei(xsimplex, y, nc, nry, dtol, 
     &      dduse, depthsim)
          do 70 kk=1,nry
            depth(kk) = depth(kk)+depthsim(kk)
 70       continue
          if (ddim.le.dtau) then
            do 80 kk=1,nry
              depthlocal(kk) = depthlocal(kk)+depthsim(kk) 
 80         continue
          endif
CC End of Evaluate the depth for a given simplex

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


      SUBROUTINE ldea(X, Y, dtau, nc, nt, nsamp, nrx, nry, nuse, dtol,
     & dduse, depth, depthlocal, dd, dld)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), depthlocal(nry), y(nry,nc), x(nrx, nc)
      dimension isimplex(nrx), xsimplex(nc+1, nc), depthsim(nry)

      external ldei
      external elarea
      external eldiam
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      do 5 i=1,nry
        depth(i) = dzero 
        depthlocal(i) = dzero 
 5    continue

      nsimp = 0
      ntot = 0
CCC          write(*,*) nt
CCC          write(*,*) nsamp

      do 20 while (nsimp.lt.nsamp.and.ntot.lt.nt)

        ntot = ntot+1
        do 10 i=1,nrx
          isimplex(i)=i
 10     continue

CC Evaluate the depth for a given simplex
CC          write(*,*) isimplex
        do 50 ii=1,(nc+1)
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          do 60 jj=1,nc
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue

CC calculate the dimension of the simplex
        if (nuse.eq.0) then
          call eldiam(xsimplex, nc, dduse, ddim)
        else
          call elarea(xsimplex, nc, dduse, ddim)   
        endif

CC calculate the depth
        call ldei(xsimplex, y, nc, nry, dtol,
     &    dduse, depthsim)
        do 70 kk=1,nry
          depth(kk) = depth(kk)+depthsim(kk)
 70     continue
        if (ddim.le.dtau) then
          do 80 kk=1,nry
            depthlocal(kk) = depthlocal(kk)+depthsim(kk) 
 80       continue
          nsimp = nsimp+1
        endif
CC End of Evaluate the depth for a given simplex

CCC      write(*,*) nsimp
CCC      write(*,*) ntot
 20   continue

      dld = nsimp
      dd = ntot

CCC      write(*,*) dld
CCC      write(*,*) dd

      call rndend()
      return
      end

CC calculate the depth for a single ellipsoid
      SUBROUTINE ldei(X, Y, nc, nry, dtol, dduse, depth)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

      dimension depth(nry), dmah(nry), y(nry,nc), x(nc+1, nc)
      dimension dcov(nc, nc), dmean(nc), dy(nc)
      dimension ipvt(nc), ztemp(nc), dwork(nc), ddeth(2)

      external ldcov
      external dgeco
      external dgedi
     
CC calcola media e covarianza dei vertici del simplesso
      call ldcov(x, nc, dmean, dcov)
CC      write(*,*) dmean
CC      write(*,*) dcov

CC calcola inversa della matrice di var/cov
      call dgeco(dcov,nc,nc,ipvt,rcond,ztemp)
      call dgedi(dcov,nc,nc,ipvt,ddeth,dwork,01)	
CC ora dcov contiene la matrice inversa, il determinante non e' calcolato
 
CC      write(*,*) dcov

CC calcola distanza di mahalanobis per le osservazioni 'y'
      do 100 i=1,nry
        dmah(i) = dzero
        depth(i) = dzero
        do 200 j=1,nc
          dy(j) = y(i,j) - dmean(j)
 200    continue

        do 400 l=1,nc
          dtemp = dzero
          do 500 m=1,nc
            dtemp = dtemp + dy(m)*dcov(m,l)
 500      continue
          dmah(i) = dmah(i) + dtemp*dy(l)
 400    continue
        if (dmah(i).ne.dmah(i)) then
          goto 101
        endif
        if (dmah(i).le.(dduse+dtol).and.dmah(i).ge.dzero) then
          depth(i) = duno
        endif
 101    continue
 100  continue
      return
      end

