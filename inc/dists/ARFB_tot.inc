c ARFB (rapidity dependent AFB)
      if((dabs(y4col).eq.dabs(y3col)).or.(abs(yff).lt.yffcut))then
        continue
      else if((dabs(y4col).gt.dabs(y3col)))then
        asycross(it,+1,8)=asycross(it,+1,8)+hist
        asyerror(it,+1,8)=asyerror(it,+1,8)
     &                +asycross(it,+1,8)**2
      else 
        asycross(it,-1,8)=asycross(it,-1,8)+hist
        asyerror(it,-1,8)=asyerror(it,-1,8)
     &                +asycross(it,-1,8)**2
      end if
      
