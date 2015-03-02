C fermion CM frame pseudorapidity
      if(m_eta4.eq.1)then
c plot distribution in eta4.
        do j=1,ndiv_eta4
          do i=1,it
            fxeta4(j,i)=fxeta4(j,i)*avgi/cnorm(i)/eta4w
          end do
        end do
        do j=1,ndiv_eta4
          fxeta4tot(j)=0.d0 
          do i=1,it
            fxeta4tot(j)=fxeta4tot(j)+fxeta4(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "eta4cm"'
        write(jh,*)'TITLE: "',f,' (CM) Pseudorapidity Distribution"'
        write(jh,*)'AXES: "{/Symbol h} (GeV)"',
     &' "d{/Symbol s}/{/Symbol h} (pb/GeV)"'
        do i=1,ndiv_eta4
          write(jh,964)xeta4(i),fxeta4tot(i)
        end do
        write(jh,*)'END'
      end if
