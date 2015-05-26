c AOFB (one sided AFB)
      if((dabs(y4col).eq.dabs(y3col)).or.(abs(pzff).lt.pzffcut))then
        continue
      else if((dabs(y4col).gt.dabs(y3col)))then
        asycross(it,+1,iAOFB)=asycross(it,+1,iAOFB)+hist
        asyerror(it,+1,iAOFB)=asyerror(it,+1,iAOFB)
     &                +asycross(it,+1,iAOFB)**2
      else 
        asycross(it,-1,iAOFB)=asycross(it,-1,iAOFB)+hist
        asyerror(it,-1,iAOFB)=asyerror(it,-1,iAOFB)
     &                +asycross(it,-1,iAOFB)**2
      end if
      
      
