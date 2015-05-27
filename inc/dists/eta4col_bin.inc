C fermion lab frame pseudorapidity
      if(m_eta4col.eq.1)then
c eta4col bin size.
        eta4colw=(eta4colmax-eta4colmin)/ndiv_eta4col
c generate bin in eta4col.
        do i=1,ndiv_eta4col
          xeta4col(i)=eta4colmin+eta4colw*(i-1)+eta4colw/2.d0
          do it=1,20
            fxeta4col(i,it)=0d0
          end do
        end do
      end if
