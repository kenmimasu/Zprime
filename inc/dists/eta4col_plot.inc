C fermion lab frame pseudorapidity
      if(m_eta4col.eq.1)then
c plot distribution in eta4col.
        do j=1,ndiv_eta4col
          do i=1,it
            fxeta4col(j,i)=fxeta4col(j,i)*avgi/cnorm(i)/eta4colw
          end do
        end do
        do j=1,ndiv_eta4col
          fxeta4coltot(j)=0.d0 
          do i=1,it
            fxeta4coltot(j)=fxeta4coltot(j)+fxeta4col(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "eta4"'
        write(jh,*)'TITLE: "',f,' Pseudorapidity Distribution"'
        write(jh,*)'AXES: "{/Symbol h} (GeV)"',
     &' "d{/Symbol s}/{/Symbol h} (pb/GeV)"'
        do i=1,ndiv_eta4col
          write(jh,964)xeta4col(i),fxeta4coltot(i)
        end do
        write(jh,*)'END'
      end if
