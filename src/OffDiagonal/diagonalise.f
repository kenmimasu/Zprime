      subroutine diag(m,w,gi,gf,
     &                m_eff,w_eff,gi_eff,gf_eff)
c Routine implements diagonalisation of the matrix 2 propagator system, 
c outputting momentum dependent effective masses, widths and couplings 
c to feed into HELAS subroutines
c     IN
c     real : p1(0:3),p2(0:3)        initial state madgraph momenta
c     real : m(2),w(2,2)            masses and width matrix
c     real : gi(2,2)                initial state couplings g(L/R,i)
c     real : gf(2,2)                final state couplings g(L/R,i)
     
c     OUT
c     real  : m_eff(2),w_eff(2,2)   eff. masses and width matrix
c     cmplx : gi(2,2)               eff. gi g(L/R,i)
c     cmplx : gf(2,2)               eff. gf g(L/R,i)

c     NOTE: generalised for chiral couplings, for scalar or vector 
c     interactions, have g(2,i)=0

      implicit none
c Arguments
      real*8 m(2),w(2,2),gi(2,2),gf(2,2)
c Returns
      real*8 m_eff2(2),m_eff(2),w_eff(2)
      complex*16 eig(2),gi_eff(2,2),gf_eff(2,2)
c Locals
      integer j,k
      real*8 mbar,mtilde,sbar,stilde
      real*8 s(2,2)
      complex*16 X,Xsq,i,d(2,2),dinv(2,2)
      logical decoup,diagonal
c ----------
c Begin code
c ----------
      i = dcmplx(0d0, 1d0)
c check that Z's are sufficiently coupled to initial and final state and that OD effects are non-zero
      decoup=.false.    
      do j =1,2
        decoup = decoup.or.(((gi(1,j).eq.0d0).and.(gf(1,j).eq.0d0)
     &                  .and.(gi(2,j).eq.0d0).and.(gf(2,j).eq.0d0)))
      enddo
      diagonal = (w(1,2).eq.0d0).or.(w(2,1).eq.0d0)
      if (decoup.or.diagonal) then
          d(1,1) = 1d0
          d(1,2) = 0d0
          d(2,1) = 0d0
          d(2,2) = 1d0
          dinv(1,1) = 1d0
          dinv(1,2) = 0d0
          dinv(2,1) = 0d0
          dinv(2,2) = 1d0
          do j=1,2
              m_eff(j)=m(j)
              w_eff(j)=w(j,j)
          enddo
       else
c Definitions
        mbar   = ( m(1)**2 + m(2)**2 )/2d0
        mtilde = ( m(1)**2 - m(2)**2 )/2d0
c Sigma(i,j) matrix
        s(1,1) = m(1)*w(1,1)
        s(2,2) = m(2)*w(2,2)
        s(1,2) = sqrt(mbar)*w(1,2)
        s(2,1) = sqrt(mbar)*w(2,1)
c
        sbar   = ( s(1,1) + s(2,2) )/2d0
        stilde = ( s(1,1) - s(2,2) )/2d0
        Xsq= (mtilde - i*stilde)**2 - s(1,2)*s(2,1)
        X=sqrt(Xsq)
c Eigenvalues
        eig(1) = -mbar + i*sbar + X
        eig(2) = -mbar + i*sbar - X
        do j=1,2
            m_eff2(j) = -REAL(eig(j)) ! Effective mass                    
        enddo
        do j=1,2
            if (m_eff2(j).ge.0d0) then ! Check physical soln. 
                m_eff(j) = SQRT(m_eff2(j))
                w_eff(j) = IMAG(eig(j))/m_eff(j)
            else
                m_eff(j) = 0d0
                w_eff(j) = 0d0
            endif
        enddo
c Amplitude ~ gi.P.gf -> gi.d.P(diag).dinv.gf;  P(diag)=dinv.P.d
c Rotation matrices for eigenbasis 
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
      END