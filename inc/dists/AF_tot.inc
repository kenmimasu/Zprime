c Forward Asymmetry
      if(dabs(y4col).gt.yc)then
         asycross(it,+1,iAF)=asycross(it,+1,iAF)+hist
         asyerror(it,+1,iAF)=asyerror(it,+1,iAF)
     &                 +asycross(it,+1,iAF)**2
      end if
      if(dabs(y3col).gt.yc)then
         asycross(it,-1,iAF)=asycross(it,-1,iAF)+hist
         asyerror(it,-1,iAF)=asyerror(it,-1,iAF)
     &                 +asycross(it,-1,iAF)**2
      end if
