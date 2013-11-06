C############################################################
C
C	Functions for the approximation of simplicial (local) depth
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: August, 26, 2008
C	Version: 0.2-1
C
C	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
C
C############################################################
C we have fix ldsai function so that it takes into account the fact
C that dmah may overflow or NaN. The last is check by dmah.ne.dmah


      SUBROUTINE ldsa(X, Y, dtau, nc, nt, nrx, nry, nuse, dtol,
     & depth, depthlocal, nnapprox, dapprox)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), depthlocal(nry), y(nry,nc), x(nrx, nc)
      dimension isimplex(nc+1), xsimplex(nc+1, nc), depthsim(nry)
      dimension nnapprox(nry), napprox(nry)

      external ldsai
      external ldarea
      external lddiam
      external diffvol

      dc = nc
      dapprox = dzero
      call diffvol(dc, dapprox)

      do 5 i=1,nry
        depth(i) = dzero 
        depthlocal(i) = dzero 
        nnapprox(i) = 0
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
            call lddiam(xsimplex, nc, ddim)   
          else
            call ldarea(xsimplex, nc, ddim)   
          endif

CC calculate the depth
          call ldsai(xsimplex, y, nc, nry, dtol, 
     &      napprox, dapprox, depthsim)
          do 70 kk=1,nry
            depth(kk) = depth(kk)+depthsim(kk)
            nnapprox(kk) = nnapprox(kk)+napprox(kk)
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


      SUBROUTINE ldsaa(X, Y, dtau, nc, nt, nsamp, nrx, nry, nuse, dtol,
     & depth, depthlocal, dd, dld, nnapprox, dapprox)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), depthlocal(nry), y(nry,nc), x(nrx, nc)
      dimension isimplex(nrx), xsimplex(nc+1, nc), depthsim(nry)
      dimension nnapprox(nry), napprox(nry)

      external ldsai
      external ldarea
      external lddiam
      external rndstart
      external rndend
      external rndunif
      external dgamma
      external diffvol

      call rndstart()

      dc = nc
      dapprox = dzero
      call diffvol(dc, dapprox)

      do 5 i=1,nry
        depth(i) = dzero 
        depthlocal(i) = dzero 
        nnapprox(i) = 0
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
          call lddiam(xsimplex, nc, ddim)
        else
          call ldarea(xsimplex, nc, ddim)   
        endif

CC calculate the depth
        call ldsai(xsimplex, y, nc, nry, dtol,
     &    napprox, dapprox, depthsim)
        do 70 kk=1,nry
          depth(kk) = depth(kk)+depthsim(kk)
          nnapprox(kk) = nnapprox(kk)+napprox(kk)
 70     continue
        if (ddim.le.dtau) then
          do 80 kk=1,nry
            depthlocal(kk) = depthlocal(kk)+depthsim(kk) 
 80       continue
          nsimp = nsimp+1
        endif
CC End of Evaluate the depth for a given simplex

CCC      write(*,*) 'Adesso sono quelle del simplesso'
CCC
CCC      call ldsai(xsimplex, xsimplex, nc, nc+1, 0,
CCC     &    napprox, dapprox, depthsim)


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

CC calculate the depth for a single simplex
      SUBROUTINE ldsai(X, Y, nc, nry, dtol, napprox, dapprox, depth)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

      dimension depth(nry), dmah(nry), y(nry,nc), x(nc+1, nc)
      dimension dcov(nc, nc), dmean(nc), dy(nc), napprox(nry)
      dimension ipvt(nc), ztemp(nc), dwork(nc), ddeth(2)

      external ldcov
      external dgeco
      external dgedi
      external dgamma

      dc = nc      

CC calcola media e covarianza dei vertici del simplesso
      call ldcov(x, nc, dmean, dcov)

CC calcola inversa della matrice di var/cov
      call dgeco(dcov,nc,nc,ipvt,rcond,ztemp)
      call dgedi(dcov,nc,nc,ipvt,ddeth,dwork,01)	
CC ora dcov contiene la matrice inversa, il determinante non e' calcolato
 
CC calcola distanza di mahalanobis per le osservazioni 'y'
      do 100 i=1,nry
        dmah(i) = dzero
        depth(i) = dzero
        napprox(i) = 0
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
        if (dmah(i).gt.(dc+dtol).or.dmah(i).lt.dzero) then
CC          write(*,*) 'Saltato 1'
          goto 101
        endif
        if (dmah(i).ne.dmah(i)) then
CC          write(*,*) 'Saltato 2'
          goto 101
        endif
        if (dmah(i).le.(duno/dc+dtol)) then
          depth(i) = duno
        elseif ((dmah(i).gt.(duno/dc-dtol)).and.
     &    (dmah(i).le.(dc+dtol))) then
          depth(i) = dapprox
          napprox(i) = 1
        endif
 101    continue
CC        write(*,*) dmah(i)
 100  continue
      
      return
      end

      SUBROUTINE ldcov (X, nc, dmean, dcov)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)
      parameter(ddue=2.0d00)

      dimension x(nc+1, nc), dcov(nc, nc), dmean(nc)

      dc = nc

CCC Calcolo medie
      do 1000 i=1,nc
        dmean(i) = dzero
        do 2000 j=1,(nc+1)
          dmean(i) = dmean(i)+x(j,i)
 2000   continue
        dmean(i) = dmean(i)/(dc+1)
 1000 continue

CCC Calcolo matrice var/cov

      do 3000 i=1,nc
        do 4000 j=i,nc
          dcov(i,j) = dzero
          do 5000 k=1,(nc+1)
            dcov(i,j) = dcov(i,j)+(x(k,i)-dmean(i))*(x(k,j)-dmean(j))
 5000     continue
          dcov(i,j) = dcov(i,j)/(dc+1)
          dcov(j,i) = dcov(i,j)
 4000   continue
 3000 continue

      return
      end


CCC corretta la formula il 08/08/2008
CCC aggiunto controllo nel caso di dimensione 1 la funzione riporta il valore 0

      SUBROUTINE diffvol(dc, dapprox)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

      external dgamma  

      if (dc.eq.1) then
        dapprox = dzero
      else 
        volbp = (dpi**(dc/ddue))/dgamma(dc/ddue+duno)
        dapprox = ((dc+duno)**((dc+duno)/ddue)/dgamma(dc+duno) - 
     &   dc**(-dc/ddue)*volbp)/((dc**(dc/ddue)-dc**(-dc/ddue))*volbp)
      endif
      return
      end





c$$$      SUBROUTINE dcalc(dc, dapprox, volbp, dnum, dden, dgam)
c$$$
c$$$      implicit double precision(a-h,o-z)
c$$$      implicit integer (n,i,j)
c$$$
c$$$      parameter(dzero=0.0d00)
c$$$      parameter(duno=1.0d00)
c$$$      parameter(ddue=2.0d00)
c$$$      parameter(dpi=3.141592654d00)
c$$$
c$$$      external dgamma  
c$$$
c$$$      volbp = (dpi**(dc/ddue))/dgamma(dc/ddue+duno)
c$$$      dnum = ((dc+duno)**((dc+duno)/ddue)/dgamma(dc+duno) - 
c$$$     &   dc**(-dc/ddue)*volbp)
c$$$      dden = ((dc**(dc/ddue)-dc**(-dc/ddue))*volbp)
c$$$      dapprox = dnum/dden
c$$$      dgam = dgamma(dc/ddue+duno)
c$$$      return
c$$$      end
