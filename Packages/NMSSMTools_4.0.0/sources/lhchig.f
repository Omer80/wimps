      SUBROUTINE LHCHIG(PAR,PROB)

*   Subroutine to check LHC constraints

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION PAR(*),PROB(*),SIG(3,8)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,BRBBSM
      DOUBLE PRECISION BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM,LHC_TBH
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH

      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     .      BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     .       HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     .       HCBRSUSY,HCWIDTH
      COMMON/LHCSIG/SIG

* Loop over H1, H2, H3

      DO I=1,3


       DO J=1,8
        SIG(I,J)=0d0
       ENDDO

       CALL HDECAY(SMASS(I),BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .      BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   H -> tautau
* VH:
       IF(BRLLSM.NE.0d0)SIG(I,1)=CV(I)**2*BRLL(I)/BRLLSM
* ggF:
       IF(BRLLSM.NE.0d0)SIG(I,2)=CJ(I)**2*BRLL(I)/BRLLSM
       
*   H -> bb
* VH:
       IF(BRBBSM.NE.0d0)SIG(I,3)=CV(I)**2*BRBB(I)/BRBBSM
* ttH:
       IF(BRGGSM.NE.0d0)SIG(I,4)=CU(I)**2*BRBB(I)/BRBBSM

*   H -> ZZ/WW
* VBF=VH:
       IF(BRZZSM.NE.0d0)SIG(I,5)=CV(I)**2*BRZZ(I)/BRZZSM
* ggF:
       IF(BRZZSM.NE.0d0)SIG(I,6)=CJ(I)**2*BRZZ(I)/BRZZSM
       
*   H -> gammagamma
* VBF=VH:
       IF(BRGGSM.NE.0d0)SIG(I,7)=CV(I)**2*BRGG(I)/BRGGSM
* ggF:
       IF(BRGGSM.NE.0d0)SIG(I,8)=CJ(I)**2*BRGG(I)/BRGGSM

      ENDDO

* Bound on Br(t->bH+)*BR(H+->tau nu)

      PROB(45)=DDIM(brtopbh*HCBRL/LHC_TBH(CMASS),1d0)

      END


      DOUBLE PRECISION FUNCTION LHC_TBH(M)

* ATLAS constraints on BR(t->bH+)*BR(H+->taunu), ATLAS-CONF-2011-151 tab.5

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=8)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,100d0,110d0,120d0,130d0,140d0,150d0,160d0/ 
      DATA Y/.104d0,.098d0,.095d0,.077d0,.066d0,.071d0,.052d0,.141d0/ 

      LHC_TBH=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_TBH=(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO

      END
