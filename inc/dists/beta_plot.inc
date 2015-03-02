C Lorentz boost
      if(m_beta.eq.1)then
c plot distribution in beta.
        do j=1,ndiv_beta
          do i=1,it
            fxbeta(j,i)=fxbeta(j,i)*avgi/cnorm(i)/betaw
          end do
        end do
        do j=1,ndiv_beta
          fxbetatot(j)=0.d0 
          do i=1,it
            fxbetatot(j)=fxbetatot(j)+fxbeta(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "beta"'
        write(jh,*)'TITLE: "{\Symbol b} Distribution"'
        write(jh,*)'AXES: "{\Symbol b} (GeV)"',
     &' "d{/Symbol s}/{\Symbol b} (pb/GeV)"'
        do i=1,ndiv_beta
          write(jh,964)xbeta(i),fxbetatot(i)
        end do
        write(jh,*)'END'
      end if
