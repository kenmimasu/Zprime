C PS integral params
      integer iseed
      common/rndm/iseed
      integer ierr
      common/errors/ierr
      integer npoints,nctpoints
      common/stat/npoints,nctpoints
      real*8 fac
      common/limfac/fac
C Collider information and acceptance cuts
      real*8 s
      common/par/s
      real*8 ecm_coll
      common/coll/ecm_coll
      integer icoll
      common/collider/icoll
      real*8 y,ptCut,dycut
      common/ptr/y,ptCut,dycut

C PDFs, factorisation scale, running alpha_s
      integer istructure
      common/PARTDIST/ISTRUCTURE
      real*8 QQ
      common/PDF/QQ
      real*8 rlambdaQCD4
      integer nloops
      common/ALFASTRONG/rlambdaQCD4,nloops
c Heavy neutrino parameters
      real*8 NeuMass
      common/neu/NeuMass
      integer iNeu
      common/neutype/iNeu
c final state switches, parameters
      real*8 rmf
      common/fermion/rmf
      integer jf,mf
      common/fints/jf,mf
      character*2 ff
      character*1 f
      common/fstr/ff,f
c contribution and QCD switches
      integer CONT,QCD,EW,BSM
      common/switch/CONT,QCD,EW,BSM
      include 'names.inc'