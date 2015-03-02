C anti-fermion lab frame pseudorapidity
      if(m_eta3col.eq.1)then
c plot distribution in eta3col.
        do j=1,ndiv_eta3col
          do i=1,it
            fxeta3col(j,i)=fxeta3col(j,i)*avgi/cnorm(i)/eta3colw
          end do
        end do
        do j=1,ndiv_eta3col
          fxeta3coltot(j)=0.d0 
          do i=1,it
            fxeta3coltot(j)=fxeta3coltot(j)+fxeta3col(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "eta3"'
        write(jh,*)'TITLE: "',f,'bar Pseudorapidity Distribution"'
        write(jh,*)'AXES: "{/Symbol h} (GeV)"',
     &' "d{/Symbol s}/{/Symbol h} (pb/GeV)"'
        do i=1,ndiv_eta3col
          write(jh,964)xeta3col(i),fxeta3coltot(i)
        end do
        write(jh,*)'END'
      end if
