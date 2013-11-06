C############################################################
C
C	Functions for the Tukey circular (local) depth
C	Author: Claudio Agostinelli and Mario Romanazzi
C	E-mail: claudio@unive.it
C	Date: July, 19, 2011
C	Version: 0.2
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

      SUBROUTINE LDTC (XDATA,YDATA,NXSIZE,NYSIZE,TAU,DTOL,
     & depth, ldepth)

      implicit double precision(a-h,l,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(dmuno=-1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)
CCC      parameter(dduepi=6.283185307d00)

      dimension ydata(nysize), xdata(nxsize)
      dimension depth(nysize), ldepth(nysize)
      dimension sectors(nxsize,2), freq(nxsize,4)
      dimension lsectorp(nxsize,2), lsectorm(nxsize,2) 
      dimension lfreqp(nxsize,3), lfreqm(nxsize,3)
      dimension naddpip(nxsize), naddpim(nxsize)

      dduepi = ddue*dpi

      dxsize = nxsize
      dxmax = nxsize + 10
      do 10 j=1, nysize
        depth(j)=dxsize
        ldepth(j)=dxmax
        ydata(j)= dmod(ydata(j), dduepi)
 10   continue
      do 20 i=1, nxsize
        xdata(i) = dmod(xdata(i), dduepi)
        xtemp = dmod(xdata(i)+dpi, dduepi)
        lxtempp = xdata(i)+tau
        lxtempm = xdata(i)-tau
        if (lxtempp.ge.dduepi) then
          naddpip(i) = 1
        else
          naddpip(i) = 0
        endif
        if (lxtempm.lt.dzero) then
          naddpim(i) = 1
        else
          naddpim(i) = 0
        endif
        lxtempp = dmod(lxtempp+dduepi, dduepi)
        lxtempm = dmod(lxtempm+dduepi, dduepi)
        freq(i,1) = dzero
        freq(i,2) = dzero
        freq(i,3) = dzero
        freq(i,4) = dzero
        lfreqp(i,1) = dzero
        lfreqp(i,2) = dzero
        lfreqp(i,3) = dzero
        lfreqm(i,1) = dzero
        lfreqm(i,2) = dzero
        lfreqm(i,3) = dzero
        if (xdata(i) < xtemp) then
          sectors(i,1) = xdata(i)
          sectors(i,2) = xtemp
        else
          sectors(i,1) = xtemp
          sectors(i,2) = xdata(i)
        endif
        if (xdata(i) < lxtempp) then
          lsectorp(i,1) = xdata(i)
          lsectorp(i,2) = lxtempp
        else
          lsectorp(i,1) = lxtempp
          lsectorp(i,2) = xdata(i)
        endif
        if (xdata(i) < lxtempm) then
          lsectorm(i,1) = xdata(i)
          lsectorm(i,2) = lxtempm
        else
          lsectorm(i,1) = lxtempm
          lsectorm(i,2) = xdata(i)
        endif
        do 30 j=1,nxsize
          if (dabs(xdata(j)-sectors(i,1)).lt.dtol) then
             freq(i,3) = freq(i,3)+duno
          else if (dabs(xdata(j)-sectors(i,2)).lt.dtol) then
             freq(i,4) = freq(i,4)+duno
          else if (xdata(j).gt.sectors(i,1)
     &      .and.xdata(j).lt.sectors(i,2)) then
             freq(i,1) = freq(i,1)+duno
          else
             freq(i,2) = freq(i,2)+duno
          endif
          if (naddpip(i).eq.1) then
            daddpi = dpi
            dsegno = dmuno
          else
            daddpi = dzero
            dsegno = duno
          endif
          if (dabs(xdata(j)-lsectorp(i,1)).lt.dtol) then
             lfreqp(i,2) = lfreqp(i,2)+duno
          else if (dabs(xdata(j)-lsectorp(i,2)).lt.dtol) then
             lfreqp(i,3) = lfreqp(i,3)+duno
          else if (dsegno*dmod(xdata(j)+daddpi,dduepi)
     &      .gt.dsegno*dmod(lsectorp(i,1)+daddpi,dduepi)
     &      .and.dsegno*dmod(xdata(j)+daddpi,dduepi)
     &      .lt.dsegno*dmod(lsectorp(i,2)+daddpi,dduepi)) then
             lfreqp(i,1) = lfreqp(i,1)+duno
          endif
          if (naddpim(i).eq.1) then
            daddpi = dpi
            dsegno = dmuno
          else
            daddpi = dzero
            dsegno = duno
          endif
          if (dabs(xdata(j)-lsectorm(i,1)).lt.dtol) then
             lfreqm(i,2) = lfreqm(i,2)+duno
          else if (dabs(xdata(j)-lsectorm(i,2)).lt.dtol) then
             lfreqm(i,3) = lfreqm(i,3)+duno
          else if (dsegno*dmod(xdata(j)+daddpi,dduepi)
     &      .gt.dsegno*dmod(lsectorm(i,1)+daddpi,dduepi)
     &      .and.dsegno*dmod(xdata(j)+daddpi,dduepi)
     &      .lt.dsegno*dmod(lsectorm(i,2)+daddpi,dduepi)) then
             lfreqm(i,1) = lfreqm(i,1)+duno
          endif
 30     continue
C      write(*,*) i
C      write(*,*) sectors(i,1), sectors(i,2)
C      write(*,*) freq(i,1), freq(i,2), freq(i,3), freq(i,4)
C      if (naddpip(i).eq.1) then
C        write(*,*) 'segno meno'
C      else
C        write(*,*) 'segno piu'
C      endif
C      write(*,*) lsectorp(i,1), lsectorp(i,2)
C      write(*,*) lfreqp(i,1), lfreqp(i,2), lfreqp(i,3)
C      if (naddpim(i).eq.1) then
C        write(*,*) 'segno meno'
C      else
C        write(*,*) 'segno piu'
C      endif
C      write(*,*) lsectorm(i,1), lsectorm(i,2)
C      write(*,*) lfreqm(i,1), lfreqm(i,2), lfreqm(i,3)
 20   continue
      do 40 i=1, nxsize 
        do 50 j=1, nysize
          if (dabs(ydata(j)-sectors(i,1)).lt.dtol) then
            dfreq = freq(i,3)+min(freq(i,1),freq(i,2))
          else if (dabs(ydata(j)-sectors(i,2)).lt.dtol) then
            dfreq = freq(i,4)+min(freq(i,1),freq(i,2))
          else if (ydata(j).gt.sectors(i,1)
     &      .and.ydata(j).lt.sectors(i,2)) then
             dfreq = freq(i,1)+min(freq(i,3),freq(i,4))
          else
             dfreq = freq(i,2)+min(freq(i,3),freq(i,4))
          endif
          depth(j) = min(depth(j), dfreq)

          ldfreqp = dxmax
          if (naddpip(i).eq.1) then
            daddpi = dpi
            dsegno = dmuno
          else
            daddpi = dzero
            dsegno = duno
          endif
          if (dabs(ydata(j)-lsectorp(i,1)).lt.dtol) then
            ldfreqp = lfreqp(i,1)+lfreqp(i,2)
          else if (dabs(ydata(j)-lsectorp(i,2)).lt.dtol) then
            ldfreqp = lfreqp(i,1)+lfreqp(i,3)
          else if (dsegno*dmod(ydata(j)+daddpi,dduepi)
     &       .gt.dsegno*dmod(lsectorp(i,1)+daddpi,dduepi)
     &       .and.dsegno*dmod(ydata(j)+daddpi,dduepi)
     &       .lt.dsegno*dmod(lsectorp(i,2)+daddpi,dduepi)) then
            ldfreqp = lfreqp(i,1)+min(lfreqp(i,2),lfreqp(i,3))
          endif
          ldfreqm = dxmax
          if (naddpim(i).eq.1) then
            daddpi = dpi
            dsegno = dmuno
          else
            daddpi = dzero
            dsegno = duno
          endif
          if (dabs(ydata(j)-lsectorm(i,1)).lt.dtol) then
            ldfreqm = lfreqm(i,1)+lfreqm(i,2)
          else if (dabs(ydata(j)-lsectorm(i,2)).lt.dtol) then
            ldfreqm = lfreqm(i,1)+lfreqm(i,3)
          else if (dsegno*dmod(ydata(j)+daddpi,dduepi)
     &       .gt.dsegno*dmod(lsectorm(i,1)+daddpi,dduepi)
     &       .and.dsegno*dmod(ydata(j)+daddpi,dduepi)
     &      .lt.dsegno*dmod(lsectorm(i,2)+daddpi,dduepi)) then
            ldfreqm = lfreqm(i,1)+min(lfreqm(i,2),lfreqm(i,3))
          endif
          ldfreq = min(ldfreqp, ldfreqm)
          ldepth(j) = min(ldepth(j), ldfreq)
 50     continue
 40   continue
      do 60 j=1, nysize
        if (ldepth(j).gt.dxsize) then
          ldepth(j) = dzero
        endif
 60   continue
      return
      end

