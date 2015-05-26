c Central Asymmetry
      if(dabs(y4col).le.yc)then
         asycross(it,+1,iAC)=asycross(it,+1,iAC)+hist
         asyerror(it,+1,iAC)=asyerror(it,+1,iAC)
     &                 +asycross(it,+1,iAC)**2
      end if
      if(dabs(y3col).le.yc)then
         asycross(it,-1,iAC)=asycross(it,-1,iAC)+hist
         asyerror(it,-1,iAC)=asyerror(it,-1,iAC)
     &                 +asycross(it,-1,iAC)**2
      end if
