C anti-fermion lab frame polar angle
      if(m_cost3col.eq.1)then
c plot distribution in cost3col.
        do j=1,ndiv_cost3col
          do i=1,it
            fxcost3col(j,i)=fxcost3col(j,i)*avgi/cnorm(i)/cost3colw
          end do
        end do
        do j=1,ndiv_cost3col
          fxcost3coltot(j)=0.d0 
          do i=1,it
            fxcost3coltot(j)=fxcost3coltot(j)+fxcost3col(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "cos3"'
        write(jh,*)'TITLE: "',f,'bar Polar Angle Distribution"'
        write(jh,*)'AXES: "cos{/Symbol q} (GeV)"',
     &' "d{/Symbol s}/cos{/Symbol q} (pb/GeV)"'
        do i=1,ndiv_cost3col
          write(jh,964)xcost3col(i),fxcost3coltot(i)
        end do
        write(jh,*)'END'
      end if
