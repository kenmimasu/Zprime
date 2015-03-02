      SUBROUTINE ZZwidthMulti(nzp,mass,width,neumass,ineu)
c calculates Zp width contributions from light fermions
      implicit real*8 (a-h,o-z)
c arguments
      real*8 mass(10),width(10)
      real*8 NeuMass
      integer nzprime,ineu
c
      include 'names.inc'
ccccccccccccccccccccc
c initialise
      do n=1,10
            width(n)=0.d0
      enddo
c read in widths from Models/model_widths.txt
      open(unit=44,file='Models/'//model(1:imodel)//
     &                                    '_widths.txt',status='old')
      read(44,*) width(1),width(2),width(3),width(4),width(5)

      return
      end
