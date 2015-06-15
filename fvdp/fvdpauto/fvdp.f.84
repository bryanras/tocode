C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
C   vdp :    Forced van der Pol oscillator
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
C     ---------- ---- 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(NDIM),PAR(*),F(NDIM)
C 
C     Define some things to ease the computation.
C
      SSR = U(1)*U(1) + U(2)*U(2) 
      SRT = DSQRT(SSR)
C
      F(1) = PAR(2)*U(1)*U(1)/SSR -  U(1)*U(3)/SRT - PAR(3)*U(2)
      F(2) = PAR(2)*U(1)*U(2)/SSR -  U(2)*U(3)/SRT + PAR(3)*U(1)
      F(3) = SRT +  PAR(1)*U(3)*(1.0-U(3)*U(3)/3.0) - 5.0
C 
      RETURN 
      END 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE STPNT(NDIM,U,PAR,T) 
C     ---------- ----- 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(NDIM),PAR(*) 
C     
C 
      PAR(1)=0.4 
      PAR(2)=0.0 
      PAR(3)=SQRT(0.84)
C 
      RETURN 
      END 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE BCND 
      RETURN 
      END 
C 
      SUBROUTINE ICND 
      RETURN 
      END 
C 
      SUBROUTINE FOPT 
      RETURN 
      END 
C 
      SUBROUTINE PVLS
      RETURN 
      END 
C---------------------------------------------------------------------- 
