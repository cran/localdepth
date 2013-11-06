C############################################################
C
C	Functions for the simplicial (local) depth in hyperspheres
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: September, 01, 2011
C	Version: 0.2
C
C	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
C
C############################################################

      SUBROUTINE ldsehs(X, Y, dtau, nc, nt, nrx, nry, nuse, dtol,
     & depth, depthlocal)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), depthlocal(nry), y(nry,nc), x(nrx, nc)
      dimension isimplex(nc), xsimplex(nc, nc), depthsim(nry)

      external ldseihs
      external ldareahs
      external lddiamhs

      do 5 i=1,nry
        depth(i) = dzero 
        depthlocal(i) = dzero 
 5    continue

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

CC Evaluate the depth for a given simplex
CC          write(*,*) isimplex
          nsimp = nsimp+1
          do 50 ii=1,nc
            is = isimplex(ii)
            do 60 jj=1,nc
              xsimplex(ii,jj) = x(is,jj)
 60         continue
 50       continue
CC calculate the dimension of the simplex
          if (nuse.eq.0) then
            call lddiamhs(xsimplex, nc, ddim)   
          else
            call ldareahs(xsimplex, nc, ddim)   
          endif
CC calculate the depth

CC          write(*,*) 'matrice xsimplex'
CC          do 75 ii=1,nc
CC             write(*,"(1X,100g15.5)") ( xsimplex(ii,jj), jj=1,nc )
CC 75       continue

          call ldseihs(xsimplex, y, nc, nry, dtol, 
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
          do 40 while (j.le.nc)
            isimplex(j) = isimplex(j-1)+1
            j = j+1
 40       continue
          i = i-1
        endif      
 20   continue
CC      write(*,*) nsimp
      return
      end

      SUBROUTINE ldseahs(X, Y, dtau, nc, nt, nsamp, nrx, nry, nuse,
     & dtol, depth, depthlocal, dd, dld)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)

      dimension depth(nry), depthlocal(nry), y(nry,nc), x(nrx, nc)
      dimension isimplex(nrx), xsimplex(nc, nc), depthsim(nry)

      external ldseihs
      external ldareahs
      external lddiamhs
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
        do 50 ii=1,nc
          is = (nrx-ii) * rndunif()+1
          iis = isimplex(is)
          isimplex(is) = isimplex(nrx-ii+1)
          do 60 jj=1,nc
            xsimplex(ii,jj) = x(iis,jj)
 60       continue
 50     continue

CC calculate the dimension of the simplex
        if (nuse.eq.0) then
          call lddiamhs(xsimplex, nc, ddim)
        else
          call ldareahs(xsimplex, nc, ddim)   
        endif

CC calculate the depth
        call ldseihs(xsimplex, y, nc, nry, dtol,
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
      SUBROUTINE ldseihs(X, Y, nc, nry, dtol, depth)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

      dimension depth(nry), y(nry,nc), x(nc, nc), xuno(nc, nc)
      dimension xunobis(nc, nc), xdue(nc), xtre(nc)
      dimension ipvt(nc), ztemp(nc), dwork(nc), ddeth(2)

      external dgeco
      external dgedi
      external dgemv

CC      write(*,*) 'matrice X'
CC      do 5 i=1,nc
CC        write(*,"(1X,100g15.5)") ( x(i,j), j=1,nc )
CC 5    continue
 
      do 10 j=1,nc
        do 20 i=1,(nc-1)
          xuno(j,i) = x(nc,j) - x(i,j)
CC          write(*,*) 'i=', i
CC          write(*,*) 'j=', j
CC          write(*,*) 'x(nc,j)=', x(nc,j)
CC          write(*,*) 'x(i,j)=', x(i,j)
CC          write(*,*) 'xuno(j,i)=', xuno(j,i)
 20     continue
          xdue(j) = x(nc,j)
 10   continue
      do 30 i=1,nry
CC Controlliamo che sia un vertice
        isver = 0
        do 35 j=1,nc
          isver1 = 0
          do 37 jj=1,nc
            if (dabs(y(i,jj)-x(j,jj)).le.dzero) then
              isver1 = isver1 + 1
            endif
 37       continue
          if (isver1.eq.nc) then
            isver = 1           
          endif
 35     continue
        if (isver.eq.1) then
          depth(i) = duno
CC altrimenti NON e' un vertice
        else
          depth(i) = dzero
          do 40 j=1,nc
          do 43 jj=1,(nc-1)
            xunobis(j,jj) = xuno(j,jj)
 43       continue
            xunobis(j,nc) = y(i,j)
 40       continue
CC       write(*,*) 'i=', i
CC        write(*,*) 'matrice A'
CC        do 45 ii=1,nc
CC          write(*,"(1X,100g15.5)") ( xunobis(ii,j), j=1,nc )
CC 45     continue
CC calcola inversa della xunobis
          call dgeco(xunobis,nc,nc,ipvt,rcond,ztemp)
          call dgedi(xunobis,nc,nc,ipvt,ddeth,dwork,11)
CC        write(*,*) 'matrice inversa di A'
CC        do 46 ii=1,nc
CC          write(*,"(1X,100g15.5)") ( xunobis(ii,j), j=1,nc )
CC 46     continue
CC ora xuno contiene la matrice inversa, il determinante non e' calcolato
          call dgemv('N',nc,nc,duno,xunobis,nc,
     &      xdue,1,dzero,xtre,1)
CC ora xtre contiene i coefficienti
CC          write(*,*) 'rcond=', rcond
CC          write(*,*) 'det=', ddeth(1) * 10.0**ddeth(2)
CC          write(*,*) 'coeff, beta, gamma, alpha'
CC          write(*,*) xtre
          no=0
          dsum=dzero
          do 50 j=1,(nc-1)
            dsum=dsum+xtre(j)
            if (xtre(j).lt.-dtol) then
              no=1
            endif
            if (xtre(j).gt.duno+dtol) then
              no=1
            endif
 50       continue
          if (xtre(nc).lt.-dtol) then
            no=1
          endif
          if (xtre(nc).gt.duno+dtol) then
            no=1
          endif
CC          write(*,*) 'somma =', dsum
          if (dsum.le.(duno+dtol).and.no.eq.0) then
            depth(i) = duno
          endif
        endif
 30   continue
CC      write(*,*) depth
      return
      end
