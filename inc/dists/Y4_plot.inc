C fermion CM frame rapidity
      if(m_Y4.eq.1)then
c plot distribution in Y4.
        do j=1,ndiv_Y4
          do i=1,it
            fxY4(j,i)=fxY4(j,i)*avgi/cnorm(i)/Y4w
          end do
        end do
        do j=1,ndiv_Y4
          fxY4tot(j)=0.d0 
          do i=1,it
            fxY4tot(j)=fxY4tot(j)+fxY4(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y4cm"'
        write(jh,*)'TITLE: "',f,' (CM) Rapidity Distribution"'
        write(jh,*)'AXES: "y (GeV)"',
     &' "d{/Symbol s}/y (pb/GeV)"'
        do i=1,ndiv_Y4
          write(jh,964)xY4(i),fxY4tot(i)
        end do
        write(jh,*)'END'
      end if
