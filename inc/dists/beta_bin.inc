C Lorentz boost
      if(m_beta.eq.1)then
c beta bin size.
        betaw=(betamax-betamin)/ndiv_beta
c generate bin in beta.
        do i=1,ndiv_beta
          xbeta(i)=betamin+betaw*(i-1)+betaw/2.d0
          do it=1,20
            fxbeta(i,it)=0d0
          end do
        end do
      end if
