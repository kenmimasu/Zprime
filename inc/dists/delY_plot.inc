C ffbar rapidity difference
      if(m_delY.eq.1)then
c plot distribution in delY.
        do j=1,ndiv_delY
          do i=1,it
            fxdelY(j,i)=fxdelY(j,i)*avgi/cnorm(i)/delYw
          end do
        end do
        do j=1,ndiv_delY
          fxdelYtot(j)=0.d0 
          do i=1,it
            fxdelYtot(j)=fxdelYtot(j)+fxdelY(j,i)
          end do
        end do
        write(jh,*)'DISTRIBUTION:'
        write(jh,*)'DSTRING: "dely"'
        write(jh,*)'TITLE: "{\Symbol D}y Distribution"'
        write(jh,*)'AXES: "{\Symbol D}y (GeV)"',
     &' "d{/Symbol s}/{\Symbol D}y (pb/GeV)"'
        do i=1,ndiv_delY
          write(jh,964)xdelY(i),fxdelYtot(i)
        end do
        write(jh,*)'END'
      end if
