C ffbar rapidity difference
      if(m_delY.eq.1)then
c delY bin size.
        delYw=(delYmax-delYmin)/ndiv_delY
c generate bin in delY.
        do i=1,ndiv_delY
          xdelY(i)=delYmin+delYw*(i-1)+delYw/2.d0
          do it=1,20
            fxdelY(i,it)=0d0
          end do
        end do
      end if
