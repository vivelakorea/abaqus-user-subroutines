C======================================================================
C========               YIELD FUNCTION                        =========
C======================================================================

      SUBROUTINE HILL48(STRESS, HILL_PARAM, dfds, ddfdss, EQSTRESS)
      
        REAL*8 STRESS(6), HILL_PARAM(6), 
     $          dfds(6), ddfdss(6,6), EQSTRESS, tmp66(6,6)
        PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)


        
        EQSTRESS = SQRT(
     $                   HILL_PARAM(1)*(STRESS(2) - STRESS(3))**TWO +
     $                   HILL_PARAM(2)*(STRESS(3) - STRESS(1))**TWO + 
     $                   HILL_PARAM(3)*(STRESS(1) - STRESS(2))**TWO +
     $                   TWO*HILL_PARAM(4)*STRESS(6)**TWO +
     $                   TWO*HILL_PARAM(5)*STRESS(5)**TWO +
     $                   TWO*HILL_PARAM(6)*STRESS(4)**TWO
     $                 )


      dfds(1:6) = 0.d0
      ddfdss(1:6,1:6) = 0.d0

      if(EQSTRESS.GT.1.d-6) then
      
        dfds(1) = ONE/EQSTRESS
     $   *(-HILL_PARAM(2)*(STRESS(3) - STRESS(1)) 
     $       + HILL_PARAM(3)*(STRESS(1) - STRESS(2)))
      
      
      
        dfds(2) = ONE/EQSTRESS
     $   *( HILL_PARAM(1)*(STRESS(2) - STRESS(3)) 
     $       - HILL_PARAM(3)*(STRESS(1) - STRESS(2)))

        dfds(3) = ONE/EQSTRESS
     $   *(-HILL_PARAM(1)*(STRESS(2) - STRESS(3)) 
     $       + HILL_PARAM(2)*(STRESS(3) - STRESS(1)))

        dfds(4) = ONE/EQSTRESS*(TWO*HILL_PARAM(6)*STRESS(4))
        dfds(5) = ONE/EQSTRESS*(TWO*HILL_PARAM(5)*STRESS(5))
        dfds(6) = ONE/EQSTRESS*(TWO*HILL_PARAM(4)*STRESS(6))
        
        do i = 1,6
          do j =1,6
            ddfdss(i,j) = ZERO
          enddo
        enddo

        ddfdss(1,1) =  ONE/EQSTRESS*(HILL_PARAM(2) + HILL_PARAM(3))
        ddfdss(2,2) =  ONE/EQSTRESS*(HILL_PARAM(1) + HILL_PARAM(3))
        ddfdss(3,3) =  ONE/EQSTRESS*(HILL_PARAM(1) + HILL_PARAM(2))
        ddfdss(1,2) =  ONE/EQSTRESS*(-HILL_PARAM(3))
        ddfdss(2,1) =  ONE/EQSTRESS*(-HILL_PARAM(3))
        ddfdss(1,3) =  ONE/EQSTRESS*(-HILL_PARAM(2))
        ddfdss(3,1) =  ONE/EQSTRESS*(-HILL_PARAM(2))
        ddfdss(2,3) =  ONE/EQSTRESS*(-HILL_PARAM(1))
        ddfdss(3,2) =  ONE/EQSTRESS*(-HILL_PARAM(1))
        ddfdss(4,4) =  ONE/EQSTRESS*(TWO*HILL_PARAM(6))
        ddfdss(5,5) =  ONE/EQSTRESS*(TWO*HILL_PARAM(5))
        ddfdss(6,6) =  ONE/EQSTRESS*(TWO*HILL_PARAM(4))

        call DYADIC(dfds, dfds, tmp66)

        ddfdss = ddfdss - ONE/(EQSTRESS**THREE)*tmp66

      endif
 
      RETURN
      END SUBROUTINE


C=======================================================================
C========       isotropic hardening function                  ==========
C=======================================================================
      SUBROUTINE VOCE(eqstrain, HARD_PARAM, ystress, hard_mod) 
        REAL*8 eqstrain, HARD_PARAM(3), ystress
        PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0)
        
        ystress = HARD_PARAM(1) + HARD_PARAM(2)
     $                           *(ONE - exp(-HARD_PARAM(3)*eqstrain))
        hard_mod = HARD_PARAM(2)*
     $                      HARD_PARAM(3)*exp(-HARD_PARAM(3)*eqstrain)

      RETURN
      END SUBROUTINE

C=======================================================================
C========               algorithmic moduli                    ==========
C=======================================================================

      SUBROUTINE ALGO_MOD(CE_INV, DLAM, DDFDSS, ALGMOD)
        REAL*8 CE_INV(6,6), DLAM, DDFDSS(6,6), ALGMOD(6,6), TMP(6,6)
        INTEGER ITMP(6,6)

        TMP(1:6,1:6) = CE_INV + DLAM*DDFDSS
        CALL MATINV(TMP, 6, 6, ITMP, ALGMOD)
       
      RETURN
      END SUBROUTINE

        

C=======================================================================
C========              matrix manipulations                   ==========
C=======================================================================

      SUBROUTINE DYADIC(A, B, C)
      
        REAL*8 A(6), B(6),
     $         C(6,6)
        
        do i = 1,6
          do j = 1,6
            C(i,j) = A(i)*B(j)
          enddo
        enddo

      RETURN
      END SUBROUTINE


C *******************************************************************
      SUBROUTINE M3INV(A,AINV)

C --   THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX [A]
C --   AND PLACES THE RESULT IN [AINV]. 
C --   IF DET(A) IS ZERO, THE CALCULATION
C --   IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C ----------------------------------------------------------------------
C   VARIABLES
C
      IMPLICIT REAL*8(A-H,O-Z)   
      REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)
C
C      A(3,3)   -- THE MATRIX WHOSE INVERSE IS DESIRED.
C      DET      -- THE COMPUTED DETERMINANT OF [A].
C      ACOFAC(3,3)   -- THE MATRIX OF COFACTORS OF A(I,J).
C               THE SIGNED MINOR (-1)**(I+J)*M_IJ
C               IS CALLED THE COFACTOR OF A(I,J).
C      AADJ(3,3)   -- THE ADJOINT OF [A]. IT IS THE MATRIX
C               OBTAINED BY REPLACING EACH ELEMENT OF
C               [A] BY ITS COFACTOR, AND THEN TAKING
C               TRANSPOSE OF THE RESULTING MATRIX.
C      AINV(3,3)   -- RETURNED AS INVERSE OF [A].
C               [AINV] = [AADJ]/DET.
C ----------------------------------------------------------------------

      CALL MDET(A,DET)
   
      IF ( DET .EQ. 0.0 ) THEN
        WRITE(91,10)
C       CALL XIT
      ENDIF
     
      CALL MCOFAC(A,ACOFAC)
      CALL MTRANS(ACOFAC,AADJ)
      DO  I = 1,3
        DO  J = 1,3
          AINV(I,J) = AADJ(I,J)/DET
        ENDDO
      ENDDO

C ----------------------------------------------------------------------
C   FORMAT
C

 10   FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +       10X,'PROGRAM TERMINATED')
C----------------------------------------------------------------------
      RETURN
      END

C **********************************************************************
      SUBROUTINE MCOFAC(A,ACOFAC)
C --
C --   THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C --   AND PLACES THE RESULT IN ACOFAC. 
C -----------------------------------------------------------------
C   VARIABLES
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8  A(3,3), ACOFAC(3,3)
C ------------------------------------------------------------------
C   COMPUTATION
C
      ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
      ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
      ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
      ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
      ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
      ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
      ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
   
      RETURN
      END

C ****************************************************************

C ******************************************************************
      SUBROUTINE MTRANS(A,ATRANS)
C --
C --   THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C --   MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C -----------------------------------------------------------------
C   VARIABLES
C
         IMPLICIT REAL*8(A-H,O-Z)
         REAL*8  A(3,3), ATRANS(3,3)
C ------------------------------------------------------------------
C   COMPUTATION
C
      CALL ONEM(ATRANS)
      DO I = 1, 3
        DO J = 1, 3
          ATRANS(I,J) = A(J,I)
        ENDDO
      ENDDO
      
      RETURN
      END
C ****************************************************************
C ******************************************************************
      SUBROUTINE MDET(A,DET)
C --
C --   THIS SUBROUTINE CALCULATES THE DETERMINANT
C --   OF A 3 BY 3 MATRIX [A].
C ------------------------------------------------------------------
C   VARIABLES
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(3,3)
C ------------------------------------------------------------------
C   COMPUTATION
C
      DET =     A(1,1)*A(2,2)*A(3,3) 
     +   + A(1,2)*A(2,3)*A(3,1)
     +   + A(1,3)*A(2,1)*A(3,2)
     +   - A(3,1)*A(2,2)*A(1,3)
     +   - A(3,2)*A(2,3)*A(1,1)
     +   - A(3,3)*A(2,1)*A(1,2)
   
      RETURN
      END
C *******************************************************************

C ******************************************************************
      SUBROUTINE ONEM(A)
C
C   THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C   3 BY 3 MATRIX [A]
C ------------------------------------------------------------------
         IMPLICIT REAL*8(A-H,O-Z)
              REAL*8 A(3,3)
      
         DO I = 1,3
           DO J = 1,3
              IF (I .EQ. J) THEN      
                  A(I,J)      = 1.0D0
              ELSE
                  A(I,J)      = 0.0D0
              END IF
            ENDDO
         ENDDO      
         RETURN
         END
C ******************************************************************

      SUBROUTINE VOT2IN(A, B)
        REAL*8 A(6), B(3,3)
        
        B(1,1) = A(1); B(1,2) = A(4); B(1,3) = A(5);
        B(2,1) = A(4); B(2,2) = A(2); B(2,3) = A(6);
        B(3,1) = A(5); B(3,2) = A(6); B(3,3) = A(3);

      RETURN
      END SUBROUTINE

C ******************************************************************

      SUBROUTINE IN2VOT(A, B)
        REAL*8 A(3,3), B(6)
        
        B(1) = A(1,1)
        B(2) = A(2,2)
        B(3) = A(3,3)
        B(4) = A(1,2)
        B(5) = A(1,3)
        B(6) = A(2,3)

      RETURN
      END SUBROUTINE

C ******************************************************************

      SUBROUTINE MATINV(A,N,NP,INDX,Y)

C
C	Given an NxN matrix [A], with physical dimension NP, this 
C	routine replaces it by the LU decomposition of a rowwise 
C	permutation of itself. [A] and N are input. [A] is output, 
C	arranged in LU form. INDX is an output vector which records
C	the row permutation effected by the partial pivoting; 
C	D is output as +1 or -1 depending on wheter the number of
C	row interchanges was even or odd, respectively.
C
C	Once the LU decomposition is performed, this routine 
C	calculates the inverse of [A] by using subroutine LUBKSB.
C	Note that INDX is input as the permutation vector returned by 	
C	LUDCMP. {B} is input as the right-hand side vector {B}, and 	
C	returns with the solution vector {X}. [A], N, NP, INDX are not 
C	modified by this routine, and are left in place
C	for succesive calls with different right-hand sides {B}.
C	This routine takes into account that {B} will begin with
C	many zero elements, so it is efficient for use in matrix 
C	inversion.
C
C	The inverse of [A] is calculated using as many unit vectors
C	{B} needed as right hand side vectors. The result is
C	returned as the matrix [Y].
C
	      IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(NP,NP), Y(NP,NP), INDX(NP)
C
C	Set up the identity matrix  
C
	      DO I = 1,N
	        DO J = 1,N
	          Y(I,J) = 0.D0
          ENDDO
        Y(I,I) = 1.D0
        ENDDO
C
C	Decompose the matrix just once
C
	      CALL LUDCMP(A,N,NP,INDX,D)
C
C	Find the inverse by columns. It is necessary to recognize
C	that FORTRAN stores two dimensional matrices by column, so
C	so that Y(1,J) is the address of the Jth column of Y.
C
	      DO J=1,N
	        CALL LUBKSB(A,N,NP,INDX,Y(1,J))
        ENDDO

	    RETURN
      END SUBROUTINE

C ******************************************************************
      
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
C
C	Given an NxN matrix [A], with physical dimension NP, this 
C	routine replaces it by the LU decomposition of a row-wise 
C	permutation of itself. [A] and N are input. [A] is output, 
C	arranged in LU form. INDX is an output vector which records
C	the row permutation effected by the partial pivoting; 
C	D is output as +1 or -1 depending on wheter the nuber of
C	row interchanges was even or odd, respectively. This routine
C	is used in combination with LUBKSB to solve linear equations 
C	or invert a matrix.
C
        IMPLICIT REAL*8 (A-H,O-Z)
        
        PARAMETER (NMAX=100,TINY=1.0E-20)
        DIMENSION A(NP,NP),INDX(N),VV(NMAX)
        D=1.D0
        DO I=1,N
          AAMAX=0.D0
          DO J=1,N
            IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
          ENDDO
          IF (AAMAX.EQ.0.) THEN
            PAUSE 'Singular matrix.'
          ENDIF
          VV(I)=1.D0/AAMAX
        ENDDO
        DO J=1,N
          IF (J.GT.1) THEN
            DO I=1,J-1
              SUM=A(I,J)
              IF (I.GT.1)THEN
                DO K=1,I-1
                  SUM=SUM-A(I,K)*A(K,J)
                ENDDO
                A(I,J)=SUM
              ENDIF
            ENDDO
          ENDIF
          AAMAX=0.D0
          DO I=J,N
            SUM=A(I,J)
            IF (J.GT.1)THEN
              DO K=1,J-1
                SUM=SUM-A(I,K)*A(K,J)
              ENDDO
              A(I,J)=SUM
            ENDIF
            DUM=VV(I)*DABS(SUM)
            IF (DUM.GE.AAMAX) THEN
              IMAX=I
              AAMAX=DUM
            ENDIF
          ENDDO
          IF (J.NE.IMAX)THEN
            DO K=1,N
              DUM=A(IMAX,K)
              A(IMAX,K)=A(J,K)
              A(J,K)=DUM
            ENDDO
            D=-D
            VV(IMAX)=VV(J)
          ENDIF
          INDX(J)=IMAX
          IF(J.NE.N)THEN
            IF(A(J,J).EQ.0.) THEN 
              A(J,J)=TINY
            ENDIF
            DUM=1./A(J,J)
            DO I=J+1,N
              A(I,J)=A(I,J)*DUM
            ENDDO
          ENDIF
        ENDDO
        IF(A(N,N).EQ.0.) THEN
          A(N,N)=TINY
        ENDIF
      RETURN
      END

C ******************************************************************
C ******************************************************************
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
C
C	Solves the set of N linear equations [A]{X} = {B}. 
C	Here [A] is input, not as the matrix [A], but as its LU 
C	decomposition, determined by the routine LUDCMP. INDX
C	is input as the permutation vector returned by LUDCMP. {B}
C	is input as the right-hand side vector {B}, and returns
C	with the solution vector {X}. [A], N, NP, INDX are not 
C	modified by this routine, and can be left in place
C	for succesive calls with different right-hand sides {B}.
C	This routine takes into account that {B} will begin with
C	many zero elements, so it is efficient for use in matrix 
C	inversion.
C
	      IMPLICIT REAL*8 (A-H,O-Z)
	      DIMENSION A(NP,NP),INDX(N),B(N)
        II=0
        DO I=1,N
          LL=INDX(I)
          SUM=B(LL)
          B(LL)=B(I)
          IF (II.NE.0)THEN
            DO J=II,I-1
              SUM=SUM-A(I,J)*B(J)
            ENDDO
          ELSE IF (SUM.NE.0.) THEN
            II=I
          ENDIF
          B(I)=SUM
        ENDDO
        DO I=N,1,-1
         SUM=B(I)
         IF(I.LT.N)THEN
           DO J=I+1,N
             SUM=SUM-A(I,J)*B(J)
           ENDDO
         ENDIF
         B(I)=SUM/A(I,I)
        ENDDO

      RETURN
      END

C ******************************************************************
C ******************************************************************
      ! SUBROUTINE DOT_PRODUCT(A, B, C)
      !   REAL*8 A(6), B(6), C

      !   C = 0.D0
      !   DO I = 1,6
      !     C = C + A(I)*B(I)
      !   ENDDO

      ! RETURN
      ! END SUBROUTINE