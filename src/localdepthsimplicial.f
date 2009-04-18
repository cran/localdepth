C############################################################
C
C	Functions for the simplicial (local) depth
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: September, 25, 2008
C	Version: 0.1
C
C	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
C
C############################################################

      SUBROUTINE ldse(X, Y, dtau, nc, nt, nrx, nry, nuse, dtol,
     & depth, depthlocal)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), depthlocal(nry), y(nry,nc), x(nrx, nc)
      dimension isimplex(nc+1), xsimplex(nc+1, nc), depthsim(nry)

      external ldsei
      external ldarea
      external lddiam

      dc = nc

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
            call lddiam(xsimplex, nc, ddim)   
          else
            call ldarea(xsimplex, nc, ddim)   
          endif

CC calculate the depth
          call ldsei(xsimplex, y, nc, nry, dtol, 
     &      depthsim)
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


      SUBROUTINE ldsea(X, Y, dtau, nc, nt, nsamp, nrx, nry, nuse, dtol,
     & depth, depthlocal, dd, dld)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), depthlocal(nry), y(nry,nc), x(nrx, nc)
      dimension isimplex(nrx), xsimplex(nc+1, nc), depthsim(nry)

      external ldsei
      external ldarea
      external lddiam
      external rndstart
      external rndend
      external rndunif

      call rndstart()

      dc = nc

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
          call lddiam(xsimplex, nc, ddim)
        else
          call ldarea(xsimplex, nc, ddim)   
        endif

CC calculate the depth
        call ldsei(xsimplex, y, nc, nry, dtol,
     &    depthsim)
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
      SUBROUTINE ldsei(X, Y, nc, nry, dtol, depth)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

      dimension depth(nry), y(nry,nc), x(nc+1, nc), xuno(nc, nc)
      dimension xdue(nc, nry), xtre(nc, nry)
      dimension ipvt(nc), ztemp(nc), dwork(nc), ddeth(2)

      external dgeco
      external dgedi
      external dgemm

      dc = nc      

      do 10 i=1,nc
        do 20 j=1,nc
          xuno(j,i) = x(i,j) - x(nc+1,j) 
 20     continue
 10   continue

      do 30 i=1,nry
        depth(i) = dzero
        do 40 j=1,nc
          xdue(j,i) = y(i,j) - x(nc+1,j) 
 40     continue
 30   continue
      
CC calcola inversa della xuno
      call dgeco(xuno,nc,nc,ipvt,rcond,ztemp)
      call dgedi(xuno,nc,nc,ipvt,ddeth,dwork,01)	
CC ora xuno contiene la matrice inversa, il determinante non e' calcolato
      call dgemm('N','N',nc,nry,nc,duno,xuno,nc,
     & xdue,nc,dzero,xtre,nc)
CC ora xtre contiene i coefficienti

      do 50 i=1,nry
        negat=0
        dsum=dzero
        do 60 j=1,nc
          dsum=dsum+xtre(j,i)
          if (xtre(j,i).le.-dtol) then
            negat=1
          endif
 60     continue
          if (dsum.le.(duno+dtol).and.negat.eq.0) then
            depth(i) = duno
          endif
 50   continue
      return
      end

