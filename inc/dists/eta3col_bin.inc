C anti-fermion lab frame pseudorapidity
      if(m_eta3col.eq.1)then
c eta3col bin size.
        eta3colw=(eta3colmax-eta3colmin)/ndiv_eta3col
c generate bin in eta3col.
        do i=1,ndiv_eta3col
          xeta3col(i)=eta3colmin+eta3colw*(i-1)+eta3colw/2.d0
          do it=1,20
            fxeta3col(i,it)=0d0
          end do
        end do
      end if
