C fermion CM frame polar angle
      if(m_cost4.eq.1)then
c cost4 bin size.
        cost4w=(cost4max-cost4min)/ndiv_cost4
c generate bin in cost4.
        do i=1,ndiv_cost4
          xcost4(i)=cost4min+cost4w*(i-1)+cost4w/2.d0
          do it=1,20
            fxcost4(i,it)=0d0
          end do
        end do
      end if
