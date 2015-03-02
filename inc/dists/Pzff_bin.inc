C ffbar system z-momentum
      if(m_Pzff.eq.1)then
c Pzff bin size.
        Pzffw=(Pzffmax-Pzffmin)/ndiv_Pzff
c generate bin in Pzff.
        do i=1,ndiv_Pzff
          xPzff(i)=Pzffmin+Pzffw*(i-1)+Pzffw/2.d0
        end do
      end if
