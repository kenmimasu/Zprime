C anti-fermion lab frame rapidity
      if(m_Y3col.eq.1)then
c plot distribution in Y3col.
        do j=1,ndiv_Y3col
          do i=1,it
            fxY3col(j,i)=fxY3col(j,i)*avgi/cnorm(i)/Y3colw
          end do
        end do
        do j=1,ndiv_Y3col
          fxY3coltot(j)=0.d0 
          do i=1,it
            fxY3coltot(j)=fxY3coltot(j)+fxY3col(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y3"'
        write(jh,*)'TITLE: "',f,'bar Rapidity Distribution"'
        write(jh,*)'AXES: "y (GeV)"',
     &' "d{/Symbol s}/y (pb/GeV)"'
        do i=1,ndiv_Y3col
          write(jh,964)xY3col(i),fxY3coltot(i)
        end do
        write(jh,*)'END'
      end if
