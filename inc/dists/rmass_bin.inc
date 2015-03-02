C ffbar invariant mass
      if(m_rmass.eq.1)then
c rmass bin size.
        rmassw=(rmassmax-rmassmin)/ndiv_rmass
c generate bin in rmass.
        do i=1,ndiv_rmass
          xrmass(i)=rmassmin+rmassw*(i-1)+rmassw/2.d0
        end do
      end if
