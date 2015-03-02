c AFBSTAR for different yff cuts
      if(costst.gt.0.d0)then
        asycross(it,+1,10)=asycross(it,+1,10)+hist
        asyerror(it,+1,10)=asyerror(it,+1,10)
     &                +asycross(it,+1,10)**2
         do j=1,8
           if(dabs(yff).gt.yffs(j))then
               asycross(it,+1,10+j)=asycross(it,+1,10+j)+hist
               asyerror(it,+1,10+j)=asyerror(it,+1,10+j)
     &                +asycross(it,+1,10+j)**2
           end if
         enddo
       else if(costst.lt.0.d0)then
        asycross(it,-1,10)=asycross(it,-1,10)+hist
        asyerror(it,-1,10)=asyerror(it,-1,10)
     &                +asycross(it,-1,10)**2
         do j=1,8
           if(dabs(yff).gt.yffs(j))then
               asycross(it,-1,10+j)=asycross(it,-1,10+j)+hist
               asyerror(it,-1,10+j)=asyerror(it,-1,10+j)
     &                +asycross(it,-1,10+j)**2
           end if
         enddo
        end if 
