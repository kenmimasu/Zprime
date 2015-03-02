C ffbar invariant mass
      if(m_rmass.eq.1)then
c plot distribution in rmass.
        do j=1,ndiv_rmass
          do i=1,it
            fxrmass(j,i)=fxrmass(j,i)*avgi/cnorm(i)/rmassw
          end do
        end do
        do j=1,ndiv_rmass
          fxrmasstot(j)=0.d0 
          do i=1,it
            fxrmasstot(j)=fxrmasstot(j)+fxrmass(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "M',ff,'"'
        write(jh,*)'TITLE: "',mff,' Distribution"'
        write(jh,*)'AXES: "',mff,' (GeV)"',
     &' "d{/Symbol s}/',mff,' (pb/GeV)"'
        do i=1,ndiv_rmass
          write(jh,964)xrmass(i),fxrmasstot(i)
        end do
        write(jh,*)'END'
      end if
