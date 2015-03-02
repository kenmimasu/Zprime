c zprime model params
      real *8 param(10),rmZZ(10),gamZZ(10),gp(10)
      integer nzp
      common/pars/param,nzp
      common/ZZ/rmZZ,gamZZ,gp
c Possible corrections to Z couplings from Z-Z' mixing
      REAL*8 DgZu(2),DgZd(2),DgZe(2),DgZn(2),
     &       DgZt(2),DgZb(2),DgZta(2),DgZnt(2) 
      common/ZCORR/DgZu,DgZd,DgZe,DgZn,DgZt,DgZb,DgZta,DgZnt
