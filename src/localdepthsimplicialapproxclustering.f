C############################################################
C
C	Functions for the approximation of simplicial (local) depth
C       it report membership for each simplex/point and dimension
C       used in the construct of a clustering procedure entirely based on depth 
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: June, 01, 2010
C	Version: 0.1
C
C	Copyright (C) 2010 Claudio Agostinelli and Mario Romanazzi
C
C############################################################
C we have fix ldsai function so that it takes into account the fact
C that dmah may overflow or NaN. The last is check by dmah.ne.dmah


      SUBROUTINE ldsac(X, Y, dtau, nc, nt, nrx, nry, nuse, dtol,
     & depth, depthlocal, dapprox)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nt,nry+1), depthlocal(nt,nry+1)
      dimension y(nry,nc), x(nrx, nc)
      dimension isimplex(nc+1), xsimplex(nc+1, nc), depthsim(nry)
      dimension napprox(nry)
      external ldsai
      external ldarea
      external lddiam
      external diffvol

      dc = nc
      dapprox = dzero
      call diffvol(dc, dapprox)

      do 4 j=1,nt
        do 5 i=1,nry
          depth(j,i) = dzero 
          depthlocal(j,i) = dzero 
 5      continue
 4    continue
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
            depth(nsimp,kk) = depthsim(kk)
 70       continue
          depth(nsimp,nry+1) = ddim
          if (ddim.le.dtau) then
            do 80 kk=1,nry
              depthlocal(nsimp,kk) = depthsim(kk) 
 80         continue
          depthlocal(nsimp,nry+1) = ddim
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


      SUBROUTINE ldsaac(X, Y, dtau, nc, nt, nsamp, nrx, nry, nuse, dtol,
     & depth, depthlocal, dd, dld, dapprox)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nsamp,nry+1), depthlocal(nsamp,nry+1) 
      dimension y(nry,nc), x(nrx, nc)
      dimension isimplex(nrx), xsimplex(nc+1, nc), depthsim(nry)
      dimension napprox(nry)

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

      do 4 j=1,nsamp
        do 5 i=1,nry
          depth(j,i) = dzero
          depthlocal(j,i) = dzero
CC         nnapprox(j,i) = 0
 5      continue
 4    continue

      nsimp = 0
      ntot = 0
CCC          write(*,*) nt
CCC          write(*,*) nsamp

      do 20 while (nsimp.lt.nsamp)

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

        if (ntot.le.nsamp) then
CCC          write(*,*) 'ntot' 
CCC          write(*,*) ntot
          do 70 kk=1,nry
            depth(ntot,kk) = depthsim(kk)
CC            nnapprox(ntot,kk) = napprox(kk)
 70       continue
          depth(ntot,nry+1) = ddim 
        endif
        if (ddim.le.dtau) then
          nsimp = nsimp+1
CCC          write(*,*) 'nsimp' 
CCC          write(*,*) nsimp
          do 80 kk=1,nry
            depthlocal(nsimp,kk) = depthsim(kk) 
 80       continue
          depthlocal(nsimp,nry+1) = ddim
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
