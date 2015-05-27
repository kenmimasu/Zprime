C ffbar system rapidity
      if(m_Yff.eq.1)then
c Yff bin size.
        Yffw=(Yffmax-Yffmin)/ndiv_Yff
c generate bin in Yff.
        do i=1,ndiv_Yff
          xYff(i)=Yffmin+Yffw*(i-1)+Yffw/2.d0
          do it=1,20
            fxYff(i,it)=0d0
          end do
        end do
      end if
