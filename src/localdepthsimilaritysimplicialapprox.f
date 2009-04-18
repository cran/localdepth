C############################################################
C
C	Functions for the approximation of simplicial (local) depth
C       for the similarity matrix
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: August, 27, 2008
C	Version: 0.1-1
C
C	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
C
C############################################################

      SUBROUTINE ldssa(X, Y, dtau, nc, nt, nrx, nry, nuse, dtol, 
     & depth, depthlocal, nnapprox, dapprox)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry,nry), depthlocal(nry,nry)
      dimension y(nry,nc), x(nrx, nc)
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
        do 6 j=1,nry
          depth(i,j) = dzero 
          depthlocal(i,j) = dzero
 6      continue 
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
            do 75 kkk=kk,nry
              depth(kk,kkk) = depth(kk,kkk) +
     &          depthsim(kk)*depthsim(kkk)
              depth(kkk,kk) = depth(kk,kkk)
 75         continue
            nnapprox(kk) = nnapprox(kk)+napprox(kk)
 70       continue
          if (ddim.le.dtau) then
            do 80 kk=1,nry
              do 85 kkk=1,nry
                depthlocal(kk,kkk) = depthlocal(kk,kkk) +
     &            depthsim(kk)*depthsim(kkk) 
                depthlocal(kkk,kk) = depthlocal(kk,kkk)
 85           continue
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


      SUBROUTINE ldssaa(X, Y, dtau, nc, nt, nsamp, nrx, nry, nuse, dtol, 
     & depth, depthlocal, dd, dld, nnapprox, dapprox)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry,nry), depthlocal(nry,nry)
      dimension y(nry,nc), x(nrx, nc)
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
        do 6 j=1,nry
          depth(i,j) = dzero 
          depthlocal(i,j) = dzero
 6      continue 
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
          do 75 kkk=kk,nry
            depth(kk,kkk) 
     &        = depth(kk,kkk)+depthsim(kk)*depthsim(kkk)
            depth(kkk,kk) = depth(kk,kkk) 
 75       continue
          nnapprox(kk) = nnapprox(kk)+napprox(kk)
 70     continue
        if (ddim.le.dtau) then
          do 80 kk=1,nry
            do 85 kkk=1,nry
              depthlocal(kk,kkk) 
     &          = depthlocal(kk,kkk)+depthsim(kk)*depthsim(kkk) 
              depthlocal(kkk,kk) = depthlocal(kk,kkk) 
 85         continue
 80       continue
          nsimp = nsimp+1
        endif
CC End of Evaluate the depth for a given simplex

 20   continue
CC      write(*,*) nsimp

      dld = nsimp
      dd = ntot

      call rndend()
      return
      end

