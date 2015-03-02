c Central Asymmetry
      if(dabs(y4col).le.yc)then
         asycross(it,+1,5)=asycross(it,+1,5)+hist
         asyerror(it,+1,5)=asyerror(it,+1,5)
     &                 +asycross(it,+1,5)**2
      end if
      if(dabs(y3col).le.yc)then
         asycross(it,-1,5)=asycross(it,-1,5)+hist
         asyerror(it,-1,5)=asyerror(it,-1,5)
     &                 +asycross(it,-1,5)**2
      end if
