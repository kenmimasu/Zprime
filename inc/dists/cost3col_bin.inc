C anti-fermion lab frame polar angle
      if(m_cost3col.eq.1)then
c cost3col bin size.
        cost3colw=(cost3colmax-cost3colmin)/ndiv_cost3col
c generate bin in cost3col.
        do i=1,ndiv_cost3col
          xcost3col(i)=cost3colmin+cost3colw*(i-1)+cost3colw/2.d0
        end do
      end if
