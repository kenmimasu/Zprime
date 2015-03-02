C fermion CM frame polar angle
      if(m_cost4.eq.1)then
c plot distribution in cost4.
        do j=1,ndiv_cost4
          do i=1,it
            fxcost4(j,i)=fxcost4(j,i)*avgi/cnorm(i)/cost4w
          end do
        end do
        do j=1,ndiv_cost4
          fxcost4tot(j)=0.d0 
          do i=1,it
            fxcost4tot(j)=fxcost4tot(j)+fxcost4(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "cos4cm"'
        write(jh,*)'TITLE: "',f,' (CM) Polar Angle Distribution"'
        write(jh,*)'AXES: "cos{/Symbol q} (GeV)"',
     &' "d{/Symbol s}/cos{/Symbol q} (pb/GeV)"'
        do i=1,ndiv_cost4
          write(jh,964)xcost4(i),fxcost4tot(i)
        end do
        write(jh,*)'END'
      end if
