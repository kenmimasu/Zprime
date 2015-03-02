c Subroutine modifies t, b, W, H, masses from standard values if needed by model
c also toggles the jf switch to 5 or 6 for non universal couplings to heavy quarks
c finally calculates any deviations in Z couplings, from file read in.
      SUBROUTINE SMCORR(rmt,gamt,rmb,rmW,rmH,
     &           DgZu,DgZd,DgZe,DgZn,DgZt,DgZb,DgZta,DgZnt,iuniv)
      implicit none
C Arguments
      real*8 rmt,gamt,rmb,rmW,rmH
      real*8 DgZu(2),DgZd(2),DgZe(2),DgZn(2),
     &       DgZt(2),DgZb(2),DgZta(2),DgZnt(2) 
      integer i,iuniv
c do nothing to masses
c set all corrective factors to 1
      do i=1,2
            DgZu(i)=1.d0
            DgZd(i)=1.d0
            DgZe(i)=1.d0
            DgZn(i)=1.d0
            DgZt(i)=1.d0
            DgZb(i)=1.d0
            DgZta(i)=1.d0
            DgZnt(i)=1.d0
      end do
      RETURN
      END