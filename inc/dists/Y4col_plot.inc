C fermion lab frame rapidity
      if(m_Y4col.eq.1)then
c plot distribution in Y4col.
        do j=1,ndiv_Y4col
          do i=1,it
            fxY4col(j,i)=fxY4col(j,i)*avgi/cnorm(i)/Y4colw
          end do
        end do
        do j=1,ndiv_Y4col
          fxY4coltot(j)=0.d0 
          do i=1,it
            fxY4coltot(j)=fxY4coltot(j)+fxY4col(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y4"'
        write(jh,*)'TITLE: "',f,' Rapidity Distribution"'
        write(jh,*)'AXES: "y (GeV)"',
     &' "d{/Symbol s}/y (pb/GeV)"'
        do i=1,ndiv_Y4col
          write(jh,964)xY4col(i),fxY4coltot(i)
        end do
        write(jh,*)'END'
      end if
