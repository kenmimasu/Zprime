C ffbar system z-momentum
      if(m_Pzff.eq.1)then
c plot distribution in Pzff.
        do j=1,ndiv_Pzff
          do i=1,it
            fxPzff(j,i)=fxPzff(j,i)*avgi/cnorm(i)/Pzffw
          end do
        end do
        do j=1,ndiv_Pzff
          fxPzfftot(j)=0.d0 
          do i=1,it
            fxPzfftot(j)=fxPzfftot(j)+fxPzff(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "pz',ff,'"'
        write(jh,*)'TITLE: "p^{z}_{',ff,'} Distribution"'
        write(jh,*)'AXES: "p^{z}_{',ff,'} (GeV)"',
     &' "d{/Symbol s}/p^{z}_{',ff,'} (pb/GeV)"'
        do i=1,ndiv_Pzff
          write(jh,964)xPzff(i),fxPzfftot(i)
        end do
        write(jh,*)'END'
      end if
