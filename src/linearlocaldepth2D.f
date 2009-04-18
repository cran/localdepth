C############################################################
C
C	Functions for the linear (local) depth
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: November, 19, 2008
C	Version: 0.1
C
C	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
C
C############################################################

      SUBROUTINE lldal2D(X, Y, nrx, nry, dtau, nsamp, nmax,
     &  dsimp, dtot, depthlocal)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depthlocal(nry), y(nry,2), x(nrx, 2)
      dimension isimplex(nrx), xsimplex(3, 2)

      external lld2D
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      do 10 i=1,nry
        depthlocal(i) = dzero 
 10   continue
      do 20 i=1,nrx
        isimplex(i)=i
 20   continue
      nsimp = 0
      ntot = 0
CCC          write(*,*) nt
CCC          write(*,*) nsamp
      do 30 while (nsimp.lt.nsamp.and.ntot.lt.nmax)

        ntot = ntot+1

CC Extract randomly one simplex
        do 50 ii=1,3
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          isimplex(nrx-ii+1) = iis
          do 60 jj=1,2
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue

CC        write(*,*) isimplex

CC Evaluate the depth for a given simplex
        ntau=0
        call lld2D(xsimplex, y, nry, dtau, ntau, depthlocal)
        nsimp = nsimp+ntau
CC End of Evaluate the depth for a given simplex

CCC      write(*,*) nsimp
CCC      write(*,*) ntot
 30   continue

      dsimp = nsimp
      dtot = ntot

      call rndend()
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC Monte Carlo

      SUBROUTINE lldmc2D(X, Y, nrx, nry, dtau, nsamp, nmax,
     &  dsimp, dtot, depthlocal)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depthlocal(nry), y(nry,2), x(nrx, 2)
      dimension isimplex(nrx), xsimplex(3, 2)

      external lld2D
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      do 10 i=1,nry
        depthlocal(i) = dzero 
 10   continue
      do 20 i=1,nrx
        isimplex(i)=i
 20   continue
      nsimp = 0
      ntot = 0
CCC          write(*,*) nt
CCC          write(*,*) nsamp
      do 30 while (nsimp.lt.nsamp.and.ntot.lt.nmax)

        ntot = ntot+1

CC Extract randomly one simplex
        do 50 ii=1,3
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          isimplex(nrx-ii+1) = iis
          do 60 jj=1,2
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue

CC        write(*,*) isimplex

CC Evaluate the depth for a given simplex
        ntau=0
        call lld2D(xsimplex, y, nry, dtau, ntau, depthlocal)
        nsimp = nsimp+ntau
CC End of Evaluate the depth for a given simplex

CCC      write(*,*) nsimp
CCC      write(*,*) ntot
 30   continue

      dsimp = nsimp
      dtot = ntot

      call rndend()
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC Similarity Monte Carlo

      SUBROUTINE lldmcs2D(X, Y, nrx, nry, dtau, nsamp, nmax,
     &  dsimp, dtot, depthsim)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depthlocal(nry), y(nry,2), x(nrx, 2)
      dimension isimplex(nrx), xsimplex(3, 2)
      dimension depthsim(nry,nry)

      external lld2D
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      dry = nry

      do 10 i=1,nry
        depthlocal(i) = dzero 
        do 15 j=1,nry
          depthsim(i,j) = dzero 
 15     continue
 10   continue

      do 20 i=1,nrx
        isimplex(i)=i
 20   continue
      nsimp = 0
      ntot = 0
CCC          write(*,*) nt
CCC          write(*,*) nsamp
      do 30 while (nsimp.lt.nsamp.and.ntot.lt.nmax)

        ntot = ntot+1

CC Extract randomly one simplex
        do 50 ii=1,3
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          isimplex(nrx-ii+1) = iis
          do 60 jj=1,2
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue

CC        write(*,*) isimplex

CC Evaluate the depth for a given simplex
        ntau=0

        do 70 i=1,nry
          depthlocal(i) = dzero 
 70     continue

        call lld2D(xsimplex, y, nry, dtau, ntau, depthlocal)
        nsimp = nsimp+ntau
CC End of Evaluate the depth for a given simplex

CC        dpesi = dzero
CC        do 75 ki=1,nry
CC          dpesi = dpesi+depthlocal(ki)
CC 75     continue
CC        dpesi = dpesi/dry

        do 80 i=1,nry
          do 90 j=1,nry
            depthsim(i,j) = depthsim(i,j) + depthlocal(i)*depthlocal(j)
CCCC     &      *dpesi
 90       continue
 80     continue

CCC      write(*,*) nsimp
CCC      write(*,*) ntot
 30   continue

      dsimp = nsimp
      dtot = ntot

      call rndend()
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC Dimension Monte Carlo

      SUBROUTINE lldmcd2D(X, nrx, nsamp, dimtube)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      dimension x(nrx, 2), dimtube(nsamp*3), dlength(3)
      dimension isimplex(nrx), xsimplex(3, 2)

      external lld2Ddim
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      ipos = 0
      do 20 i=1,nrx
        isimplex(i)=i
 20   continue

CC      write(*,*) nsamp
      do 30 j=1,nsamp
CC        write(*,*) j
CC Extract randomly one simplex
        do 50 ii=1,3
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          isimplex(nrx-ii+1) = iis
          do 60 jj=1,2
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue

CC        write(*,*) isimplex

CC Evaluate the dimension for a given simplex
        call lld2Ddim(xsimplex, dlength)
CC End of Evaluate the depth for a given simplex

        do 70 k=1,3
          ipos = ipos+1
          dimtube(ipos) = dlength(k)
 70     continue
 30   continue


      call rndend()
      return
      end

CCCCCCCCCCCCCCCCC Only dimension of the tube

      SUBROUTINE lld2Ddim(X, dlength)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      dimension x(3, 2), dlength(3)

      juno = 1
      jdue = 2
      jtre = 3
      do 100 i=1,3
        da = (x(juno,2)+x(jdue,2))/ddue
        db = (x(juno,1)+x(jdue,1))/ddue
        dc = da*x(jtre,1)-db*x(jtre,2)
        da = x(jtre,2)-da
        db = db-x(jtre,1)

CC        write(*,*) da
CC        write(*,*) db
CC        write(*,*) dc

        dlength(i) = dabs(da*x(juno,1)+db*x(juno,2)+dc)
        jquattro = jtre 
        jtre = jdue
        jdue = juno
        juno = jquattro
 100  continue
 
CC      write(*,*) dlength

      return
      end

CC calculate the linear local depth for a single simplex
      SUBROUTINE lld2D(X, Y, nry, dtau, ntau, depthlocal)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      dimension x(3, 2), y(nry,2), depthlocal(nry)

      juno = 1
      jdue = 2
      jtre = 3
      do 100 i=1,3
        da = (x(juno,2)+x(jdue,2))/ddue
        db = (x(juno,1)+x(jdue,1))/ddue
        dc = da*x(jtre,1)-db*x(jtre,2)
        da = x(jtre,2)-da
        db = db-x(jtre,1)
        dlength = dabs(da*x(juno,1)+db*x(juno,2)+dc)
        if (dlength.le.dtau) then
          ntau = 1
          do 200 j=1,nry
            dy = dabs(da*y(j,1)+db*y(j,2)+dc)
            if (dy.le.dlength) then
              depthlocal(j) = depthlocal(j) + duno
            endif
 200      continue
        endif
        jquattro = jtre 
        jtre = jdue
        jdue = juno
        juno = jquattro
 100  continue
      return
      end




CCCCCCCCCCCCCCCCCCCCCCCCCCCC depth
CCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC Monte Carlo

      SUBROUTINE dlmc2D(X, Y, nrx, nry, nsamp, depth)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), y(nry,2), x(nrx, 2)
      dimension isimplex(nrx), xsimplex(3, 2)

      external lld2D
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      do 10 i=1,nry
        depth(i) = dzero 
 10   continue
      do 20 i=1,nrx
        isimplex(i)=i
 20   continue

CCC          write(*,*) nsamp
      do 30 k=1,nsimp

CC Extract randomly one simplex
        do 50 ii=1,3
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          isimplex(nrx-ii+1) = iis
          do 60 jj=1,2
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue

CC        write(*,*) isimplex

CC Evaluate the depth for a given simplex
        ntau=0
        call dl2D(xsimplex, y, nry, depth)
CC End of Evaluate the depth for a given simplex
 30   continue
      call rndend()
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC Similarity Monte Carlo

      SUBROUTINE dlmcs2D(X, Y, nrx, nry, nsamp, depthsim)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), y(nry,2), x(nrx, 2)
      dimension isimplex(nrx), xsimplex(3, 2)
      dimension depthsim(nry,nry)

      external lld2D
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      dry = nry

      do 10 i=1,nry
        depth(i) = dzero 
        do 15 j=1,nry
          depthsim(i,j) = dzero 
 15     continue
 10   continue

      do 20 i=1,nrx
        isimplex(i)=i
 20   continue

      do 30 k=1,nsamp

CC Extract randomly one simplex
        do 50 ii=1,3
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          isimplex(nrx-ii+1) = iis
          do 60 jj=1,2
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue

CC        write(*,*) isimplex

CC Evaluate the depth for a given simplex
        do 70 i=1,nry
          depth(i) = dzero 
 70     continue

        call dl2D(xsimplex, y, nry, depth)
CC End of Evaluate the depth for a given simplex

CC        dpesi = dzero
CC        do 75 ki=1,nry
CC          dpesi = dpesi+depth(ki)
CC 75     continue
CC        dpesi = dpesi/dry

        do 80 i=1,nry
          do 90 j=1,nry
            depthsim(i,j) = depthsim(i,j) + depth(i)*depth(j)
CCCCC*dpesi
 90       continue
 80     continue

 30   continue
      call rndend()
      return
      end

CC calculate the linear depth for a single simplex
      SUBROUTINE dl2D(X, Y, nry, depth)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      dimension x(3, 2), y(nry,2), depth(nry)

      juno = 1
      jdue = 2
      jtre = 3
      do 100 i=1,3
        da = (x(juno,2)+x(jdue,2))/ddue
        db = (x(juno,1)+x(jdue,1))/ddue
        dc = da*x(jtre,1)-db*x(jtre,2)
        da = x(jtre,2)-da
        db = db-x(jtre,1)
        dlength = dabs(da*x(juno,1)+db*x(juno,2)+dc)
        do 200 j=1,nry
          dy = dabs(da*y(j,1)+db*y(j,2)+dc)
          if (dy.le.dlength) then
           depth(j) = depth(j) + duno
          endif
 200    continue
        jquattro = jtre 
        jtre = jdue
        jdue = juno
        juno = jquattro
 100  continue
      return
      end


