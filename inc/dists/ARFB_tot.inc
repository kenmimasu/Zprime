c ARFB (rapidity dependent AFB)
      if((dabs(y4col).eq.dabs(y3col)).or.(abs(yff).lt.yffcut))then
        continue
      else if((dabs(y4col).gt.dabs(y3col)))then
        asycross(it,+1,iARFB)=asycross(it,+1,iARFB)+hist
        asyerror(it,+1,iARFB)=asyerror(it,+1,iARFB)
     &                +asycross(it,+1,iARFB)**2
      else 
        asycross(it,-1,iARFB)=asycross(it,-1,iARFB)+hist
        asyerror(it,-1,iARFB)=asyerror(it,-1,iARFB)
     &                +asycross(it,-1,iARFB)**2
      end if
      
