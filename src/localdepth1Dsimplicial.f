      SUBROUTINE LD1DS (XDATA,YDATA,NXSIZE,NYSIZE,NSIMP,TAU,NUSE,
     & ldepth, depth, diameter)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	ld1ds function
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: August, 14, 2007
C	Version: 0.2
C
C	Copyright (C) 2007 Claudio Agostinelli and Mario Romanazzi
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; version 2 of the License.
C
C   This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    NUSE=1 use diameter
C    NUSE=2 use area
C    NUSE=3 use spherical

      implicit double precision(a-h,l,o-z)
      implicit integer (n,i,j)

      dimension ydata(nysize),xdata(nxsize) 
      dimension diameter(nsimp)
      dimension depth(nysize), ldepth(nysize)

      do 10 ind2=1, nysize
        depth(ind2)=0.0d00
        ldepth(ind2)=0.0d00
 10   continue

      ind1=0
      do 20 i1=1, nxsize-1
        do 30 i2=i1+1, nxsize
          ind1=ind1+1
          x1=min(xdata(i1),xdata(i2))
          x2=max(xdata(i1),xdata(i2))
          diameter(ind1)=x2-x1
          do 40 ind2=1, nysize 
              if (ydata(ind2).ge.x1.and.ydata(ind2).le.x2) then
                 depth(ind2)=depth(ind2)+1.0d00
                 if (nuse.eq.1.or.nuse.eq.2) then
                    if (diameter(ind1).le.tau) then
                      ldepth(ind2)=ldepth(ind2)+1.0d00
                    endif
                 else
                   spherical = max(x2-ydata(ind2),ydata(ind2)-x1)
                   if (spherical.le.tau) then
                     ldepth(ind2)=ldepth(ind2)+1.0d00
                   endif
                 endif
              endif
 40       continue
 30     continue
 20   continue
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LD1DSALL (XDATA,YDATA,NXSIZE,NYSIZE,NSIMP,
     & depth, diameter)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	ld1dsall function
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: August, 10, 2007
C	Version: 0.1
C
C	Copyright (C) 2007 Claudio Agostinelli and Mario Romanazzi
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; version 2 of the License.
C
C   This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      dimension ydata(nysize),xdata(nxsize) 
      dimension diameter(nsimp)
      dimension depth(nsimp,nysize)

      ind1=0
      do 10 i1=1, nxsize-1
        do 20 i2=i1+1, nxsize
          ind1=ind1+1
          x1=min(xdata(i1),xdata(i2))
          x2=max(xdata(i1),xdata(i2))
          diameter(ind1)=x2-x1
          do 30 ind2=1, nysize
            depth(ind1,ind2)=0
            if (ydata(ind2).ge.x1.and.ydata(ind2).le.x2) then
               depth(ind1,ind2)=1
            endif
 30       continue
 20     continue
 10   continue
      return
      end
