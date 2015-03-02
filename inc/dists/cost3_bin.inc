C anti-fermion CM frame polar angle
      if(m_cost3.eq.1)then
c cost3 bin size.
        cost3w=(cost3max-cost3min)/ndiv_cost3
c generate bin in cost3.
        do i=1,ndiv_cost3
          xcost3(i)=cost3min+cost3w*(i-1)+cost3w/2.d0
        end do
      end if
