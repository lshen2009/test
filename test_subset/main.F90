program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE
  

  REAL(kind=dp) :: JVS(LU_NONZERO_12),JVS_orig(LU_NONZERO_12)
  REAL(kind=dp) :: VAR_Selected(LU_NSEL_12),VAR_Selected(LU_NDEL_12),X(NVAR)
  REAL(kind=dp) :: JVS1(LU_NONZERO_12),JVS2(LU_NONZERO_12)
  INTEGER:: IER,i,j,k

  CALL initialize_1D(JVS_orig)
  CALL initialize_1D(X)
  
  !First check the Function
    

end program
