program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE gckpp_Function
  USE gckpp_Function2
  USE gckpp_Jacobian
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE
  

  REAL(kind=dp) :: X(NVAR)
  REAL(kind=dp) :: F(NFIX),RCT(NREACT),Vdot(NVAR),diff(NVAR)
  REAL(kind=dp) :: Prate(NVAR),Lrate(NVAR),Lrate2(NVAR),Vdot2(NVAR),diff(NVAR)
  INTEGER:: IER,i,j,k

  CALL initialize_1D(X)
  CALL initialize_1D(F)
  CALL initialize_1D(RCT)
  Vdot=0
  diff=99
  
  !First check the Function  
  CALL Fun ( X, F, RCT, Vdot)
  CALL Fun2 ( X, F, RCT, Vdot2, Prate, Lrate, Lrate2)
  diff=Vdot2-Vdot
  print *, diff  

end program
