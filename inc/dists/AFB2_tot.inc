c FB asymmetry. IN CM FRAME!
      if(cost4.eq.0.d0)then
        continue
      else if(cost4.gt.0.d0)then
        asycross(it,+1,7)=asycross(it,+1,7)+fxn1*wgt 
        asyerror(it,+1,7)=asyerror(it,+1,7)
     &                 +asycross(it,+1,7)**2
        asycross(it,-1,7)=asycross(it,-1,7)+fxn2*wgt
        asyerror(it,-1,7)=asyerror(it,-1,7)
     &                 +asycross(it,-1,7)**2
      else if(cost4.lt.0.d0)then
        asycross(it,-1,7)=asycross(it,-1,7)+fxn1*wgt
        asyerror(it,-1,7)=asyerror(it,-1,7)
     &                 +asycross(it,-1,7)**2
        asycross(it,+1,7)=asycross(it,+1,7)+fxn2*wgt
        asyerror(it,+1,7)=asyerror(it,+1,7)
     &                 +asycross(it,-1,7)**2
      end if
