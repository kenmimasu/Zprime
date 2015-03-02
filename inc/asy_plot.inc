      if((m_sigp.eq.1).and.(m_sigm.eq.1))then
c plot distribution in asymmetries.
          do j =1,nasy
          ave(j)=1.d0
          enddo
 1399   FORMAT(F10.4,F15.4,2E15.5)
        do 1234 jasy=1,nasy
        snorm(jasy)=0.d0
        do j=1,ndiv_sigp
          do i=1,it
            fxsigp(jasy,j,i)=fxsigp(jasy,j,i)*avgi/cnorm(i)/sigpw
          end do
        end do
        do j=1,ndiv_sigp
          fxsigptot(jasy,j)=0.d0 
          do i=1,it
            fxsigptot(jasy,j)=fxsigptot(jasy,j)+fxsigp(jasy,j,i)
          end do
        end do
        do j=1,ndiv_sigm
          do i=1,it
            fxsigm(jasy,j,i)=fxsigm(jasy,j,i)*avgi/cnorm(i)/sigmw
          end do
        end do
        do j=1,ndiv_sigm
          fxsigmtot(jasy,j)=0.d0 
          do i=1,it
            fxsigmtot(jasy,j)=fxsigmtot(jasy,j)+fxsigm(jasy,j,i)
          end do
        end do
        write(jh,*)'ASYMMETRY:'
 1425   FORMAT('TITLE: "',A10,' Distribution"')
 1426   FORMAT('TITLE: "',A10,'(',F3.1,') Distribution"')
 1427   FORMAT('TITLE: "',A10,'(',F3.0,') Distribution"')
 1428   FORMAT('AXES: "',A6,' (GeV)" "',A10,'"')
        if(jasy.eq.1)then
              write(jh,*)'DSTRING: "ALL"'
              write(jh,1425)'A_{LL}'
              write(jh,1428)mff,'A_{LL}'
        else if(jasy.eq.2)then
              write(jh,*)'DSTRING: "AL"'
              write(jh,1425)'A_{L}'
              write(jh,1428)mff,'A_{L}'
        else if(jasy.eq.3)then
              write(jh,*)'DSTRING: "APV"'
              write(jh,1425)'A_{PV}'
              write(jh,1428)mff,'A_{PV}'
        else if(jasy.eq.4)then
              write(jh,*)'DSTRING: "AFB2"'
              write(jh,1425)'A_{FB}(CM)'
              write(jh,1428)mff,'A_{FB}'
        else if(jasy.eq.5)then
              write(jh,*)'DSTRING: "AC"'
              write(jh,1426)'A_{C}',yc
              write(jh,1428)mff,'A_{C}'
        else if(jasy.eq.6)then
              write(jh,*)'DSTRING: "AF"'
              write(jh,1426)'A_{F}',yc
              write(jh,1428)mff,'A_{F}'
        else if(jasy.eq.7)then
              write(jh,*)'DSTRING: "AFB"'
              write(jh,1425)'A_{FB}'
              write(jh,1428)mff,'A_{FB}'
        else if(jasy.eq.8)then
              write(jh,*)'DSTRING: "ARFB"'
              write(jh,1426)'A_{RFB}',yffcut
              write(jh,1428)mff,'A_{RFB}'
        else if(jasy.eq.9)then
              write(jh,*)'DSTRING: "AOFB"'
              write(jh,1427)'A_{OFB}',pff
              write(jh,1428)mff,'A_{OFB}'
        else if(jasy.eq.10)then
              write(jh,*)'DSTRING: "AFBSTAR"'
              write(jh,1425)'A^{*}_{FB}'
              write(jh,1428)mff,'A^{*}_{FB}'
        end if
        do j=1,8
          if(jasy.eq.(10+j))then
            write(intstr,'(i1)') j 
            write(jh,*)'DSTRING: "AFBSTAR'//intstr//'"'
            write(jh,1426)'A^{*}_{FB}',yffs(j)
            write(jh,1428)'A^{*}_{FB}'
          endif
        enddo
        ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
        do i=1,ndiv_sig
          if(fxsigptot(jasy,i)+fxsigmtot(jasy,i).eq.0.d0)then
            write(jh,1399)(xsigm(i)+xsigp(i))/2.d0,0.d0,
     &                fxsigptot(jasy,i),
     &                fxsigmtot(jasy,i)           
            snorm(jasy)=snorm(jasy)+0.d0
          else  
            write(jh,1399)(xsigm(i)+xsigp(i))/2.d0,
     &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
     &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))/ave(jasy),
     &                fxsigptot(jasy,i),
     &                fxsigmtot(jasy,i)
            snorm(jasy)=snorm(jasy)+
     &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
     &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))/ave(jasy)
     &               *fxrmasstot(i)*rmassw/avgi
          end if
        end do
        write(jh,"(A12,F10.4)")'INTEGRATED:',snorm(jasy)
        write(jh,*)'END'
        
 1234   continue
      end if