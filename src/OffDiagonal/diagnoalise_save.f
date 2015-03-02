      subroutine diag(p1,p2,m,w,gi,gf,
     &                m_eff,w_eff,gi_eff,gf_eff,c_eff)
c Routine implements diagonalisation of the matrix 2 propagator system, 
c outputting momentum dependent effective masses, widths and couplings 
c to feed into HELAS subroutines
c     IN
c     real : p1(0:3),p2(0:3)        initial state madgraph momenta
c     real : m(2),w(2,2)            masses and width matrix
c     real : gi(2,2)                initial state couplings g(L/R,i)
c     real : gf(2,2)                final state couplings g(L/R,i)
     
c     OUT
c     cmplx : m_eff(2),w_eff(2,2)   eff. masses and width matrix
c     cmplx : gi(2,2)               eff. gi g(L/R,i)
c     cmplx : gf(2,2)               eff. gf g(L/R,i)
c     cmplx : c_1,c_2               complex coeffs. for eff. couplings

c     NOTE: generalised for chiral couplings, for scalar or vector 
c     interactions, have g(2,i)=0

      implicit none
c Arguments
      real*8 p1(0:3),p2(0:3),m(2),w(2,2),gi(2,2),gf(2,2)
c Returns
      real*8 m_eff2
      complex*16 m_eff,w_eff
      complex*16 gi_eff(2,2),gf_eff(2,2),c_eff(2)

c Locals
      integer j,k
      real*8 mbar,mtilde,sbar,stilde,psq
      real*8 s(2,2)
      complex*16 X,Xsq,i,d(2,2),dinv(2,2)
c ----------
c Begin code
c ----------
c Definitions
      i = dcmplx(0d0, 1d0)
c Sigma(i,j) matrix
      s(1,1) = m(1)*w(1,1)
      s(2,2) = m(2)*w(2,2)
      s(1,2) = m(1)*w(1,2)
      s(2,1) = m(2)*w(2,1)
c
      mbar   = ( m(1)**2 + m(2)**2 )/2d0
      mtilde = ( m(1)**2 - m(2)**2 )/2d0
      sbar   = ( s(1,1) + s(2,2) )/2d0
      stilde = ( s(1,1) - s(2,2) )/2d0
      Xsq= (mtilde - i*stilde)**2 - s(1,2)*s(2,1)
      X=sqrt(Xsq)
c Amplitude ~ gi.P.gf -> gi.d.P(diag).dinv.gf;  P(diag)=dinv.P.d
c Rotation matrices for eigenbasis 
      if ((w(1,2).eq.0d0).or.(w(2,1).eq.0d0)) then
          d(1,1) =1d0
          d(1,2) =0d0
          d(2,1) =0d0
          d(2,2) =1d0
          dinv(1,1) =1d0
          dinv(1,2) =0d0
          dinv(2,1) =0d0
          dinv(2,2) =1d0
      else
          d(1,1) = ( stilde + i*( mtilde - X ) )/s(2,1)
          d(1,2) = ( stilde + i*( mtilde + X ) )/s(2,1)
          d(2,1) = dcmplx( 1d0, 0d0 )
          d(2,2) = dcmplx( 1d0, 0d0 )
          dinv(1,1) = i*s(2,1)/X/2d0
          dinv(1,2) = (mtilde + X - i*stilde)/X/2d0
          dinv(2,1) = -i*s(2,1)/X/2d0
          dinv(2,2) = (-mtilde + X + i*stilde)/X/2d0
      endif    
c rotate coupling basis like (gi).d and dinv.gf
      do j=1,2
          do k=1,2
              gi_eff(k,j)=gi(k,1)*d(1,j)+gi(k,2)*d(2,j)
              gf_eff(k,j)=dinv(j,1)*gf(k,1)+dinv(j,2)*gf(k,2)
          enddo
      enddo
C       do j=1,2
C           do k=1,2
C               if (abs(REAL(gi_eff(j,k))).lt.(1d-10)) then
C                   print*,'gi(',j,k,') = ', gi(j,k)
C                   print*,'gi_eff(',j,k,') = ', gi_eff(j,k)
C               endif
C               if (abs(REAL(gf_eff(j,k))).lt.(1d-10)) then
C                   print*,'gf(',j,k,') = ',gf(j,k)
C                   print*,'gf_eff(',j,k,') = ', gf_eff(j,k)
C               endif
C           enddo
C       enddo
c sqrt(s)
      psq = (P1(0)+P2(0))**2 - (P1(1)+P2(1))**2 
     &    - (P1(2)+P2(2))**2 - (P1(3)+P2(3))**2 
c Effective coupling, mass, width
      c_eff(1) = ( 1d0 - (mbar - i*sbar + X)/psq )
      c_eff(2) = ( 1d0 - (mbar - i*sbar - X)/psq )
      m_eff2 = 2d0*mbar - ( m(1)**2*m(2)**2 + s(1,2)*s(2,1)
     &                                     - s(1,1)*s(2,2) )/psq
      if (m_eff2.lt.0d0) then
          m_eff = dcmplx(0d0,sqrt(abs(m_eff2)))
      else
          m_eff = dcmplx(sqrt(abs(m_eff2)),0d0)
      endif
      w_eff = dcmplx(2d0*sbar - 
     &          ( m(1)**2*s(2,2) + m(2)**2*s(1,1) )/psq, 0d0)
      w_eff = w_eff/m_eff
C       print*, 'm_eff,w_eff = ',m_eff,w_eff
      
C       print*,'m',m
C       print*,'w',w
C       print*,'mbar',mbar
C       print*,'psq',psq
C       print*,'mtilde',mtilde
C       print*,'m_eff',m_eff

      END