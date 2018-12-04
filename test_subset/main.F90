program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE gckpp_Function
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE
  

  REAL(kind=dp) :: JVS(LU_NONZERO_2),JVS_orig(LU_NONZERO_2)
  REAL(kind=dp) :: X_selected(LU_NSEL_2),X_deleted(LU_NDEL_2),X(NVAR)
  REAL(kind=dp) :: JVS1(LU_NONZERO_2),JVS2(LU_NONZERO_2)
  REAL(kind=dp) :: F(NFIX),RCT(NREACT),Vdot_SEL(LU_NSEL_2),Vdot(NVAR),diff(LU_NSEL_2)
  INTEGER:: IER,i,j,k

  CALL initialize_1D(JVS_orig)
  CALL initialize_1D(X)
  CALL initialize_1D(F)
  CALL initialize_1D(RCT)
  Vdot=0
  Vdot_SEL=0
  
  X_selected=X(select_ind_2)
  X_deleted=X(delete_ind_2)
  !First check the Function
  CALL Fun_2 ( X_selected,X_deleted, F, RCT, Vdot_SEL)
  CALL Fun_13 ( X, F, RCT, Vdot, NVAR )
  diff=Vdot(select_ind_2)-Vdot_SEL  
  print *, SUM(ABS(diff))

end program
