program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE

  REAL(kind=dp) :: JVS(LU_NONZERO_12),JVS1(LU_NONZERO_12),JVS2(LU_NONZERO_12)
  REAL(kind=dp) :: JVS_orig(LU_NONZERO_12)
  REAL(kind=dp) :: X(LU_NSEL_12),t1,t2,deltaT1,deltaT2
  LOGICAL :: ind(LU_NONZERO_12)
  INTEGER:: IER,i,j,k

  CALL initialize_1D(JVS_orig)
  CALL initialize_1D(X)
  JVS1=JVS_orig
  CALL KppDecomp(JVS1,IER,LU_NSEL_12,LU_NONZERO_12,LU_CROW_12,LU_DIAG_12,LU_ICOL_12)
  JVS2=JVS_orig
  CALL KppDecomp_12 (JVS2,IER)
  JVS=JVS2-JVS1
  print *,"difference:", SUM(ABS(JVS))
 
  JVS1=JVS_orig
  IER=0
  CALL CPU_TIME(time=t1)
  DO i=1,1000
  CALL KppDecomp(JVS1,IER,LU_NSEL_12,LU_NONZERO_12,LU_CROW_12,LU_DIAG_12,LU_ICOL_12) 
  ENDDO
  CALL CPU_TIME(time=t2)
  deltaT1=t2-t1
  print *, "KppDecomp: ",deltaT1

  JVS2=JVS_orig
  CALL CPU_TIME(time=t1)
  DO i=1,1000
  CALL KppDecomp_12 (JVS2,IER)
  END DO
  CALL CPU_TIME(time=t2)
  deltaT2=t2-t1
  print *, "KppDecomp_12: ",deltaT2
  !print *, JVS1 - JVS2
  print *, deltaT1/deltaT2
end program
