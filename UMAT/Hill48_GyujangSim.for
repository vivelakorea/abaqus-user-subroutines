CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C writen by : Gyujang Sim
C date      : 2022/12/30
C institute : Seoul National University 
C             Mechanics of Materials Laboratory
C             https://sites.google.com/view/snumml/home
C reference : [1] M. A. Crisfield, "Non-linear Finite Element Analysis
C                 of Solids and Structures"
C             [2] J. C. Simo, T. J. R. Hudges, "Computational 
C                 Inelasticity"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INCLUDE 'utils.for'

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C!!!!!!!!                 START OF UMAT                       !!!!!!!!!!

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
  
C=======================================================================
C========             VARIABLE DECLARATION                   ===========
C=======================================================================

C     DIMENSION DEFINITION
      REAL*8 E_MOD, E_NU, HARD_PARAM(3), HILL_PARAM(6), 
     1       CE(6,6), CE_INV(6,6), eqstress, eqstress0, eqstrain, 
     2       ystress, hard_mod, dfds(6), ddfdss(6,6), f, dlam, ddlam, 
     3       stress0(6), cur_stress(6), residual(6), Qmat(6,6),
     4       Qmat_INV(6,6), dstress(6), stress_tolerance, ALGMOD(6,6),
     5       UNITFOURTH(6,6), tmp6(6), tmp66(6,6), tmp3(3), tmp33(3,3)

      INTEGER itmp6(6), iter

C     PARAMETER
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0, 
     1 NEWTON=100, TOLER=1.0D-6)

      do i = 1,6
        do j = 1,6
          if (i .eq. j) then
            UNITFOURTH(i,j) = ONE
          else
            UNITFOURTH(i,j) = ZERO
          endif
        enddo
      enddo

      iter = 0
C=======================================================================
C========              VARIABLE READING                       ==========
C=======================================================================

C     ELASTIC PROPERTIES
      E_MOD = PROPS(1)
      E_NU = PROPS(2)

C     HARDENING PARAMETERS (SIGMA = H1 + H2(1 - exp(-H3*epsilon)))
      HARD_PARAM(1:3) = PROPS(3:5)
      
C     HILL'S YIELD FUNCTION
      HILL_PARAM(1:6) = PROPS(6:11)

      ! if (KINC .eq. 1) then
      !   STATEV(1)=0.D0
      ! endif

C     EQUIVALENT PLASTIC STRAIN
      eqstrain = STATEV(1)

C=======================================================================
C========                ELASTIC GUESS                        ==========
C=======================================================================

C     ELASTIC STIFFNESS

      do i=1,6
        do j=1,6
          if ((i .le. 3) .and. (j .le. 3)) then
            !    o o o x x x
            !    o o o x x x
            !    o o o x x x
            !    x x x x x x
            !    x x x x x x
            !    x x x x x x
            if (i .eq. j) then
              CE(i,j) = E_MOD*(ONE - E_NU)/((1 + E_NU)*(1 - TWO*E_NU))
              CE_INV(i,j) = ONE/E_MOD
            else
              CE(i,j)     = E_MOD*E_NU/((1 + E_NU)*(1 - TWO*E_NU))
              CE_INV(i,j) = -E_NU/E_MOD
            endif
          elseif ((i .ge. 4) .and. (j .ge. 4)) then
            !    x x x x x x
            !    x x x x x x
            !    x x x x x x
            !    x x x o o o
            !    x x x o o o
            !    x x x o o o
            if (i .eq. j) then
              CE(i,j)     = E_MOD/(TWO + TWO*E_NU)
              CE_INV(i,j) = (TWO + TWO*E_NU)/E_MOD
            else
              CE(i,j) = ZERO
              CE_INV(i,j) = ZERO
            endif
          else
          !    x x x o o o
          !    x x x o o o
          !    x x x o o o
          !    o o o x x x
          !    o o o x x x
          !    o o o x x x
            CE(i,j) = ZERO
            CE_INV(i,j) = ZERO
          endif
        enddo
      enddo

C     elastic guess
      STRESS(1:6) = STRESS + matmul(CE, DSTRAN)
      DDSDDE(1:6,1:6) = CE
      
C=======================================================================
C========              PLASTICITY CHECK                       ==========
C=======================================================================

      call HILL48(stress, HILL_PARAM, dfds, ddfdss, eqstress)
      call VOCE(EQSTRAIN, HARD_PARAM, ystress, hard_mod)

      f = eqstress - ystress
    
      stress_tolerance = TOLER*eqstress

      if (f .gt. stress_tolerance) then
        iter = 1
C       ================================================================
C       ========              RETURN MAPPING  [1]             ==========
C       ================================================================
        stress0(1:6) = stress(1:6)
        eqstress0 = eqstress  
 
        dlam = ZERO
        cur_stress(1:6) = stress(1:6)

        
        do while ( 
     $             (abs(f) .gt. stress_tolerance) 
     $             .or. 
     $             (.not. (all(abs(residual) .lt. stress_tolerance)))
     $           )
          
          call HILL48(cur_stress, HILL_PARAM, dfds, ddfdss, eqstress)
          
  
          Qmat = UNITFOURTH + dlam*matmul(CE, ddfdss)
  
          call MATINV(Qmat, 6, 6, itmp6, Qmat_INV)
  
          call VOCE(EQSTRAIN+dlam, HARD_PARAM, ystress, hard_mod)
  
          f = eqstress - ystress
  
          residual(1:6) = cur_stress - 
     $                       (stress0 - dlam*matmul(CE, dfds))
  
          ddlam = (f - dot_product(dfds, matmul(Qmat_INV, residual)))
     $       /(dot_product(dfds, matmul(matmul(Qmat_INV, CE), dfds) + 
     $                                                       hard_mod))

          dstress = -matmul(Qmat_INV, residual) - 
     $                 ddlam*matmul(matmul(Qmat_INV, CE), dfds)

          
          dlam = dlam + ddlam
          
          cur_stress = cur_stress + dstress

          ! write(*,*) 'new iter'
          ! write(*,*) iter
          ! write(*,*) stress0 
          ! write(*,*) cur_stress          
          ! write(*,*) dlam
          ! write(*,*) ddlam
          ! write(*,*) f
          ! write(*,*) residual

          if ((iter .gt. NEWTON) .or. (abs(f) .gt. eqstress0*1.D1)) then
            PNEWDT = ONE/TWO
            ! write(*,*) 'reduce time increment'
            goto 666
          endif

          iter = iter + 1
        enddo

        eqstrain = eqstrain + dlam
      
        stress(1:6) = cur_stress(1:6)

C       ================================================================
C       ======                    DDSDDE [2]                       =====
C       ================================================================

        call HILL48(stress, HILL_PARAM, dfds, ddfdss, eqstress)
        
        CALL ALGO_MOD(CE_INV, dlam, ddfdss, algmod)
        tmp6 = matmul(algmod, dfds)
        call DYADIC(tmp6, tmp6, tmp66)
        DDSDDE(1:6,1:6) = algmod - (tmp66)/
     &                           (dot_product(dfds, tmp6))
       
      endif
C=====================================================================
C========                PASS STATE VARIABLES                 ========
C=====================================================================
      
      STATEV(1) = eqstrain
      STATEV(2) = iter
  
  666 RETURN
      END


C!!!!!!!!                 END OF UMAT                         !!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










