C fermion lab frame polar angle
      if(m_cost4col.eq.1)then
c cost4col bin size.
        cost4colw=(cost4colmax-cost4colmin)/ndiv_cost4col
c generate bin in cost4col.
        do i=1,ndiv_cost4col
          xcost4col(i)=cost4colmin+cost4colw*(i-1)+cost4colw/2.d0
        end do
      end if
