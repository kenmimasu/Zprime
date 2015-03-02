c Forward Asymmetry
      if(dabs(y4col).gt.yc)then
         asycross(it,+1,6)=asycross(it,+1,6)+hist
         asyerror(it,+1,6)=asyerror(it,+1,6)
     &                 +asycross(it,+1,6)**2
      end if
      if(dabs(y3col).gt.yc)then
         asycross(it,-1,6)=asycross(it,-1,6)+hist
         asyerror(it,-1,6)=asyerror(it,-1,6)
     &                 +asycross(it,-1,6)**2
      end if
