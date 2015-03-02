C fermion lab frame rapidity
      if(m_Y4col.eq.1)then
c Y4col bin size.
        Y4colw=(Y4colmax-Y4colmin)/ndiv_Y4col
c generate bin in Y4col.
        do i=1,ndiv_Y4col
          xY4col(i)=Y4colmin+Y4colw*(i-1)+Y4colw/2.d0
        end do
      end if
