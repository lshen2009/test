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
  REAL(kind=dp) :: F(NFIX),RCT(NREACT),Vdot(NVAR),Vdot2(NVAR),diff(NVAR)
  REAL(kind=dp) :: Prate(NVAR),Lrate(NVAR),Lrate2(NVAR)
  INTEGER:: IER,i,j,k

  CALL initialize_1D(X)
  CALL initialize_1D(F)
  CALL initialize_1D(RCT)
  Vdot=0
  diff=99
  X=X*9999
  F=F*9999
  RCT=RCT*9999
  
  !First check the Function  
  CALL Fun ( X, F, RCT, Vdot)
  CALL Fun2 ( X, F, RCT, Vdot2, Prate, Lrate, Lrate2)
  diff=Prate+Lrate-Vdot
  print *, SUM(ABS(diff  ))
  diff=Lrate-Lrate2*X
  print *, SUM(ABS(diff  ))  

end program