c FB asymmetry. IN Collider FRAME!
      if(cost4col.eq.0.d0)then
        continue
      else if(cost4col.gt.0.d0)then
        asycross(it,+1,4)=asycross(it,+1,4)+hist
        asyerror(it,+1,4)=asyerror(it,+1,4)
     &                 +asycross(it,+1,4)**2
      else if(cost4col.lt.0.d0)then
        asycross(it,-1,4)=asycross(it,-1,4)+hist
        asyerror(it,-1,4)=asyerror(it,-1,4)
     &                 +asycross(it,-1,4)**2
      end if 
