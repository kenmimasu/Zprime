c AFBSTAR for different yff cuts
      if(costst.gt.0.d0)then
        asycross(it,+1,iAFBSTAR)=asycross(it,+1,iAFBSTAR)+hist
        asyerror(it,+1,iAFBSTAR)=asyerror(it,+1,iAFBSTAR)
     &                +asycross(it,+1,iAFBSTAR)**2
         do j=1,8
           if(dabs(yff).gt.yffs(j))then
               asycross(it,+1,iAFBSTAR+j)=asycross(it,+1,iAFBSTAR+j)
     &                + hist
               asyerror(it,+1,iAFBSTAR+j)=asyerror(it,+1,iAFBSTAR+j)
     &                + asycross(it,+1,iAFBSTAR+j)**2
           end if
         enddo
       else if(costst.lt.0.d0)then
        asycross(it,-1,iAFBSTAR)=asycross(it,-1,iAFBSTAR)+hist
        asyerror(it,-1,iAFBSTAR)=asyerror(it,-1,iAFBSTAR)
     &                + asycross(it,-1,iAFBSTAR)**2
         do j=1,8
           if(dabs(yff).gt.yffs(j))then
               asycross(it,-1,iAFBSTAR+j)=asycross(it,-1,iAFBSTAR+j)
     &                + hist
               asyerror(it,-1,iAFBSTAR+j)=asyerror(it,-1,iAFBSTAR+j)
     &                + asycross(it,-1,iAFBSTAR+j)**2
           end if
         enddo
        end if 
