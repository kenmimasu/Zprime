C anti-fermion CM frame pseudorapidity
      if(m_eta3.eq.1)then
c plot distribution in eta3.
        do j=1,ndiv_eta3
          do i=1,it
            fxeta3(j,i)=fxeta3(j,i)*avgi/cnorm(i)/eta3w
          end do
        end do
        do j=1,ndiv_eta3
          fxeta3tot(j)=0.d0 
          do i=1,it
            fxeta3tot(j)=fxeta3tot(j)+fxeta3(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "eta3cm"'
        write(jh,*)'TITLE: "',f,'bar (CM) Pseudorapidity Distribution"'
        write(jh,*)'AXES: "{/Symbol h} (GeV)"',
     &' "d{/Symbol s}/{/Symbol h} (pb/GeV)"'
        do i=1,ndiv_eta3
          write(jh,964)xeta3(i),fxeta3tot(i)
        end do
        write(jh,*)'END'
      end if
