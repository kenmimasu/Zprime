c AFBSTAR for different yff cuts
      if(costst.gt.0.d0)then
        asycross(it,+1,iAFBM)=asycross(it,+1,iAFBM)+hist
        asyerror(it,+1,iAFBM)=asyerror(it,+1,iAFBM)
     &                +asycross(it,+1,iAFBM)**2
         do j=1,3
           if(dabs(yff).gt.yffs(j))then
               asycross(it,+1,iAFBM+j)=asycross(it,+1,iAFBM+j)+hist
               asyerror(it,+1,iAFBM+j)=asyerror(it,+1,iAFBM+j)
     &                +asycross(it,+1,iAFBM+j)**2
           end if
         enddo
       else if(costst.lt.0.d0)then
        asycross(it,-1,iAFBM)=asycross(it,-1,iAFBM)+hist
        asyerror(it,-1,iAFBM)=asyerror(it,-1,iAFBM)
     &                +asycross(it,-1,iAFBM)**2
         do j=1,3
           if(dabs(yff).gt.yffs(j))then
               asycross(it,-1,iAFBM+j)=asycross(it,-1,iAFBM+j)+hist
               asyerror(it,-1,iAFBM+j)=asyerror(it,-1,iAFBM+j)
     &                +asycross(it,-1,iAFBM+j)**2
           end if
         enddo
        end if 
