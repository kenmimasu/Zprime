C fermion lab frame polar angle
      if(m_cost4col.eq.1)then
c plot distribution in cost4col.
        do j=1,ndiv_cost4col
          do i=1,it
            fxcost4col(j,i)=fxcost4col(j,i)*avgi/cnorm(i)/cost4colw
          end do
        end do
        do j=1,ndiv_cost4col
          fxcost4coltot(j)=0.d0 
          do i=1,it
            fxcost4coltot(j)=fxcost4coltot(j)+fxcost4col(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "cos4"'
        write(jh,*)'TITLE: "',f,' Polar Angle Distribution"'
        write(jh,*)'AXES: "cos{/Symbol q} (GeV)"',
     &' "d{/Symbol s}/cos{/Symbol q} (pb/GeV)"'
        do i=1,ndiv_cost4col
          write(jh,964)xcost4col(i),fxcost4coltot(i)
        end do
        write(jh,*)'END'
      end if
