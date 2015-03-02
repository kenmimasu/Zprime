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
C Locals
      real*8 Dgvu,Dgau,Dgvd,Dgad,Dgve,Dgae,
     &       Dgvn,Dgan,Dgvt,Dgat,Dgvb,Dgab
C Commons 
      include 'names.inc'
      include 'ewcoups.inc'
CC     
      iuniv=1
C       print*, 'reading Models/'//model(1:imodel)//'_SMcorr.txt'
c read in SMcorr file
      open(unit=43,file='Models/'//model(1:imodel)//'_SMcorr.txt'
     &                                    ,status='old')
c read in model file and set couplings
      read(43,*) rmt ! mass corrections
      read(43,*) rmb
      read(43,*) rmW  
      read(43,*) rmH
      read(43,*) Dgvu ! new vector and axial couplings of the Z
      read(43,*) Dgau
      read(43,*) Dgvd
      read(43,*) Dgad
      read(43,*) Dgve
      read(43,*) Dgae
      read(43,*) Dgvn
      read(43,*) Dgan
      read(43,*) Dgvt
      read(43,*) Dgat
      read(43,*) Dgvb
      read(43,*) Dgab
c work out corrective factors
      DgZu(1)=-(Dgvu+Dgau) ! new couplings
      DgZu(2)=-(Dgvu-Dgau)
      DgZd(1)=-(Dgvd+Dgad)
      DgZd(2)=-(Dgvd-Dgad)
      DgZe(1)=-(Dgve+Dgae)
      DgZe(2)=-(Dgve-Dgae)
      DgZn(1)=-(Dgvn+Dgan)
      DgZn(2)=-(Dgvn-Dgan)
      DgZta(1)=-(Dgve+Dgae)
      DgZta(2)=-(Dgve-Dgae)
      DgZnt(1)=-(Dgvn+Dgan)
      DgZnt(2)=-(Dgvn-Dgan)
      DgZt(1)=-(Dgvt+Dgat)
      DgZt(2)=-(Dgvt-Dgat)
      DgZb(1)=-(Dgvb+Dgab)
      DgZb(2)=-(Dgvb-Dgab)
      do i=1,2
            DgZu(i)=DgZu(i)/GZU(i) ! rescale by old couplings
            DgZd(i)=DgZd(i)/GZD(i)
            DgZe(i)=DgZe(i)/GZL(i)
            DgZn(i)=DgZn(i)/GZN(i)
            DgZt(i)=DgZt(i)/GZU(i)
            DgZb(i)=DgZb(i)/GZD(i)
      end do

      RETURN
      END