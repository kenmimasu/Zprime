C anti-fermion CM frame polar angle
      if(m_cost3.eq.1)then
c plot distribution in cost3.
        do j=1,ndiv_cost3
          do i=1,it
            fxcost3(j,i)=fxcost3(j,i)*avgi/cnorm(i)/cost3w
          end do
        end do
        do j=1,ndiv_cost3
          fxcost3tot(j)=0.d0 
          do i=1,it
            fxcost3tot(j)=fxcost3tot(j)+fxcost3(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "cos3cm"'
        write(jh,*)'TITLE: "',f,'bar (CM) Polar Angle Distribution"'
        write(jh,*)'AXES: "cos{/Symbol q} (GeV)"',
     &' "d{/Symbol s}/cos{/Symbol q} (pb/GeV)"'
        do i=1,ndiv_cost3
          write(jh,964)xcost3(i),fxcost3tot(i)
        end do
        write(jh,*)'END'
      end if
