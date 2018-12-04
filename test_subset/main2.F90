program main
  USE gckpp_Parameters
  USE gckpp_Function2 
  USE initialize
 

  IMPLICIT NONE
  

  REAL(kind=dp) :: X(NVAR)
  REAL(kind=dp) :: F(NFIX),RCT(NREACT),Vdot(NVAR),Vdot2(NVAR),diff(NVAR)
  REAL(kind=dp) :: Prate(NVAR),Lrate(NVAR),Lrate2(NVAR)
  INTEGER:: IER,i,j,k

  CALL initialize_1D(X)
  CALL initialize_1D(F)
  CALL initialize_1D(RCT)
  Vdot=0
  Vdot2=0
  Prate=0
  Lrate=0
  Lrate2=0
  
  diff=0
  X=X*999
  F=F*999
  RCT=RCT*999
  
  !First check the Function  
  CALL Fun ( X, F, RCT, Vdot)
  CALL Fun2 ( X, F, RCT, Vdot2, Prate, Lrate, Lrate2)
  print *,'----------'
  diff=0
  diff=Vdot-Vdot2
  print *, SUM(ABS(diff  ))
  print *,'----------'
  diff=0
  diff=Prate+Lrate-Vdot2
  DO i=1,234
     print *,i,diff(i)
  END DO
  print *,'---------'
  print *, SUM(ABS(diff  ))
  diff=0  
  diff=Lrate-Lrate2*X
  print *, SUM(ABS(diff  ))  

end program
