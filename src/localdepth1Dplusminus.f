C############################################################
C
C	Functions for the halfspace (local) depth
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: November, 16, 2011
C	Version: 0.1
C
C	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
C
C############################################################
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

      SUBROUTINE LD1PM (XDATA,YDATA,NXSIZE,NYSIZE,TAU,
     & poldepth, neldepth, posdepth, negdepth)

      implicit double precision(a-h,l,o-z)
      implicit integer (n,i,j)

      dimension ydata(nysize),xdata(nxsize)
      dimension posdepth(nysize), poldepth(nysize)
      dimension negdepth(nysize), neldepth(nysize)

      do 10 ind1=1, nysize
        posdepth(ind1)=0.0d00
        poldepth(ind1)=0.0d00
        negdepth(ind1)=0.0d00
        neldepth(ind1)=0.0d00
 10   continue
      do 20 ind1=1, nysize
        depthneg=0.0d00
        ldepthneg=0.0d00
        depthpos=0.0d00
        ldepthpos=0.0d00
        do 30 ind2=1, nxsize
          if (xdata(ind2).eq.ydata(ind1)) then
            depthneg = depthneg + 1.0d00
            ldepthneg = ldepthneg + 1.0d00
            depthpos = depthpos + 1.0d00
            ldepthpos = ldepthpos + 1.0d00
          else if (xdata(ind2).lt.ydata(ind1)) then
            depthneg = depthneg + 1.0d00
            dist = ydata(ind1) - xdata(ind2)
            if (dist.le.tau) then
              ldepthneg = ldepthneg + 1.0d00
            endif
          else
            depthpos = depthpos + 1.0d00
            dist = xdata(ind2) - ydata(ind1)
            if (dist.le.tau) then
              ldepthpos = ldepthpos + 1.0d00
            endif
          endif
 30     continue
        negdepth(ind1) = depthneg
        posdepth(ind1) = depthpos
        neldepth(ind1) = ldepthneg
        poldepth(ind1) = ldepthpos 
 20   continue
      return
      end
