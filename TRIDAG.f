c     $Id$
c     $Source$

      SUBROUTINE TRIDAG(A,B,C,R,U,N) !!A,B,C,R,U are N dimensional vectors
c     Solves the matrix equation MU=R where M is a tridiagonal matrix with 
c     diagonal B, sub-diagonal A? and super-diagonal C
      INTEGER NMAX
      PARAMETER (NMAX=300)
      DOUBLE PRECISION A,B,C,R,U,GAM,BET
      INTEGER N, J
      DIMENSION GAM(NMAX),A(N),B(N),C(N),R(N),U(N)
      IF (B(1) .EQ. 0.) PAUSE 'TRIDAG: rewrite eqns.'
      BET=B(1)
      U(1)=R(1)/BET
      DO 11 J=2,N
         GAM(J) = C(J-1)/BET
         BET=B(J)-A(J)*GAM(J)
         IF(BET .EQ. 0.) PAUSE 'TRIDAG failed'
         U(J)=(R(J)-A(J)*U(J-1))/BET
 11      CONTINUE
      DO 12 J=N-1,1,-1
         U(J)=U(J)-GAM(J+1)*U(J+1)
 12      CONTINUE
      RETURN
      END


