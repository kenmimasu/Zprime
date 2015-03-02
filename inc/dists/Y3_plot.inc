C anti-fermion CM frame rapidity
      if(m_Y3.eq.1)then
c plot distribution in Y3.
        do j=1,ndiv_Y3
          do i=1,it
            fxY3(j,i)=fxY3(j,i)*avgi/cnorm(i)/Y3w
          end do
        end do
        do j=1,ndiv_Y3
          fxY3tot(j)=0.d0 
          do i=1,it
            fxY3tot(j)=fxY3tot(j)+fxY3(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y3cm"'
        write(jh,*)'TITLE: "',f,'bar (CM) Rapidity Distribution"'
        write(jh,*)'AXES: "y (GeV)"',
     &' "d{/Symbol s}/y (pb/GeV)"'
        do i=1,ndiv_Y3
          write(jh,964)xY3(i),fxY3tot(i)
        end do
        write(jh,*)'END'
      end if
