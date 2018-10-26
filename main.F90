program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE

  REAL(kind=dp) :: JVS1(LU_NONZERO),JVS2(LU_NONZERO)
  REAL(kind=dp) :: JVS_orig(LU_NONZERO)
  REAL(kind=dp) :: X1(NVAR),X2(NVAR),X_orig(NVAR)
  LOGICAL :: ind(LU_NONZERO)
  INTEGER:: IER,i,j,k,ind2(LU_NONZERO)
  REAL, ALLOCATABLE :: new(:)
  !INTEGER, ALLOCATABLE :: LU_IROW2(:),LU_ICOL2(:),LU_DIAG2(:)
  INTEGER :: LU_DIAG2(NVAR+1),LU_CROW2(NVAR+1),LU_ICOL2(LU_NONZERO)

  CALL initialize_1D(JVS_orig)
  CALL initialize_1D(X_orig)
  IER=0

  JVS1=JVS_orig
  CALL KppDecomp(JVS1,IER)
  JVS2=JVS_orig
  CALL KppDecomp2(JVS2,IER)
  !print *, JVS1-JVS2

  
  OPEN(unit=1101,file='text_JVS.txt')
  write (1101,*),JVS1
  close(1101)

  OPEN(unit=1101,file='text_JVS_orig.txt')
  write (1101,*),JVS_orig
  close(1101)

  JVS1=JVS_orig
  X1=X_orig
  CALL KppSolveIndirect(JVS1,X1)

  OPEN(unit=1101,file='text_X.txt')
  write (1101,*),X1
  close(1101)

  OPEN(unit=1101,file='text_X_orig.txt')
  write (1101,*),X_orig
  close(1101) 


end program
