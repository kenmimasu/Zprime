      SUBROUTINE ZZwidthMulti(nzprime,mass,width,neumass,ineu)
c  reads in widths from file named Models/{$MODEL}_widths.txt
      implicit none 
c arguments
      real*8 mass(10),width(10)
      real*8 NeuMass
      integer nzprime,ineu,n
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
      read(44,*) width(1),width(2),width(3),width(4),width(5),
     &           width(6),width(7),width(8),width(9),width(10)  
      return
      end
