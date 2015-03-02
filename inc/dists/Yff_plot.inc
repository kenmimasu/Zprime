C ffbar system rapidity
      if(m_Yff.eq.1)then
c plot distribution in Yff.
        do j=1,ndiv_Yff
          do i=1,it
            fxYff(j,i)=fxYff(j,i)*avgi/cnorm(i)/Yffw
          end do
        end do
        do j=1,ndiv_Yff
          fxYfftot(j)=0.d0 
          do i=1,it
            fxYfftot(j)=fxYfftot(j)+fxYff(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "y',ff,'"'
        write(jh,*)'TITLE: "y_{',ff,'} Distribution"'
        write(jh,*)'AXES: "y_{',ff,'} (GeV)"',
     &' "d{/Symbol s}/y_{',ff,'} (pb/GeV)"'
        do i=1,ndiv_Yff
          write(jh,964)xYff(i),fxYfftot(i)
        end do
        write(jh,*)'END'
      end if
