c FB asymmetry. IN Collider FRAME!
      if(cost4col.eq.0.d0)then
        continue
      else if(cost4col.gt.0.d0)then
        asycross(it,+1,iAFB)=asycross(it,+1,iAFB)+hist
        asyerror(it,+1,iAFB)=asyerror(it,+1,iAFB)
     &                 +asycross(it,+1,iAFB)**2
      else if(cost4col.lt.0.d0)then
        asycross(it,-1,iAFB)=asycross(it,-1,iAFB)+hist
        asyerror(it,-1,iAFB)=asyerror(it,-1,iAFB)
     &                 +asycross(it,-1,iAFB)**2
      end if 
