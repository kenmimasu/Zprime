C fermion energy
      if(m_Et.eq.1)then
c Et bin size.
        Etw=(Etmax-Etmin)/ndiv_Et
c generate bin in Et.
        do i=1,ndiv_Et
          xEt(i)=Etmin+Etw*(i-1)+Etw/2.d0
          do it=1,20
            fxEt(i,it)=0d0
          end do
        end do
      end if
