C fermion energy
      if(m_Et.eq.1)then
c plot distribution in Et.
        do j=1,ndiv_Et
          do i=1,it
            fxEt(j,i)=fxEt(j,i)*avgi/cnorm(i)/Etw
          end do
        end do
        do j=1,ndiv_Et
          fxEttot(j)=0.d0 
          do i=1,it
            fxEttot(j)=fxEttot(j)+fxEt(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "E',f,'"'
        write(jh,*)'TITLE: "E_',f,' Distribution"'
        write(jh,*)'AXES: "E_',f,' (GeV)"',
     &' "d{/Symbol s}/E_',f,' (pb/GeV)"'
        do i=1,ndiv_Et
          write(jh,964)xEt(i),fxEttot(i)
        end do
        write(jh,*)'END'
      end if
