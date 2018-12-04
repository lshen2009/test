program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE gckpp_Function
  USE gckpp_Jacobian
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE
  

  REAL(kind=dp) :: JVS(LU_NONZERO_11),JVS_orig(LU_NONZERO)
  REAL(kind=dp) :: X_selected(LU_NSEL_11),X_deleted(LU_NDEL_11),X(NVAR)
  REAL(kind=dp) :: JVS1(LU_NONZERO_11),JVS2(LU_NONZERO)
  REAL(kind=dp) :: F(NFIX),RCT(NREACT),Vdot_SEL(LU_NSEL_11),Vdot(NVAR),diff(LU_NSEL_11)
  REAL(kind=dp) :: diff2(LU_NONZERO_11)
  REAL(kind=dp) :: Jcb1(LU_NSEL_11,LU_NSEL_11),Jcb2(NVAR,NVAR),Jcb3(LU_NSEL_11,LU_NSEL_11)
  REAL(kind=dp) :: JacDiff(LU_NSEL_11,LU_NSEL_11)
  INTEGER:: IER,i,j,k

  CALL initialize_1D(JVS_orig)
  CALL initialize_1D(X)
  CALL initialize_1D(F)
  CALL initialize_1D(RCT)
  Vdot=0
  Vdot_SEL=0
  Jcb1=0
  Jcb2=0
  
  X_selected=X(select_ind_11)
  X_deleted=X(delete_ind_11)
  !First check the Function
  CALL Fun_11 ( X_selected,X_deleted, F, RCT, Vdot_SEL)
  CALL Fun_13 ( X, F, RCT, Vdot, NVAR )
  diff=Vdot(select_ind_11)-Vdot_SEL  
  print *, "------Main_11--------"
  print *, "gckpp_Function",SUM(ABS(diff))
  

  CALL Jac_SP_11(X_selected,X_deleted, F, RCT, JVS1)
  DO i=1,LU_NONZERO_11
     Jcb1(LU_IROW_11(i),LU_ICOL_11(i)) = JVS1(i)
  END DO  
  
  CALL Jac_SP_13(X, F, RCT, JVS2)
  DO i=1,LU_NONZERO
     Jcb2(LU_IROW_13(i),LU_ICOL_13(i)) = JVS2(i)
  END DO

  DO i=1,LU_NSEL_11
    DO j=1,LU_NSEL_11
	  Jcb3(i,j)=Jcb2(select_ind_11(i),select_ind_11(j))
	END DO
  END DO
  
  JacDiff=Jcb3-Jcb1
  print *,"gckpp_Jac_SP", SUM(ABS(JacDiff))



end program
