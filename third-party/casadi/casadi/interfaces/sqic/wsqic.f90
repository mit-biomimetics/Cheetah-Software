module SQICModule
 use snModulePrecision, only : ip, rp
 use SQIC,              only : qpProb
 type(qpProb)              :: QP

 character(8), allocatable :: Names(:)
 real(rp),     allocatable :: cObj(:)

end module

subroutine wsqic (m, n, nnzA, indA, locA, valA, bl, bu, hEtype, hs, x, pi, rc, &
                  nnzH, indH, locH, valH) bind ( C, name="sqic" )

  use SQICModule

  implicit none

  integer                   :: INFO
  integer(ip)               :: Errors, iObj, iPrint, iSumm, ncObj, &
                               m, n, n_inf, nnH, nnzH, nNames, nnzA, nS, iSpecs

  real(rp)                  :: ObjAdd, sInf


  real(rp):: bl(n+m), bu(n+m), x(n+m), valA(nnzA), valH(nnzH) ,pi(m), rc(n+m)
  integer(ip):: indA(nnzA), locA(n+1), indH(nnzH), locH(n+1), hEtype(n+m), hs(n+m)

  character(8)              :: probName


  Errors   = 0

  ncObj    =  0
  nnH      = n
  nNames   =  1
  ObjAdd = 0.0

  iSumm    = 0;  iPrint = 6;

  if ( allocated(cObj) )   deallocate ( cObj )
  if ( allocated(Names) )  deallocate ( Names )

  ! Allocate space for problem.
  allocate ( cObj(ncObj) )
  allocate ( Names(nNames) )

  probName = 'QP'
  call QP%load ( probName, m, n, nnH, m, ObjAdd, &
                 nnzA, indA, locA, valA, bl, bu, ncObj, cObj, &
                 nNames, Names, hEtype, hs, x, pi, rc, nnzH, indH, locH, valH )

  ! Initialize SQIC.
  call QP%begin ( iPrint, iSumm )

  !print *, 'Reading params.spc'
  !iSpecs = 4
  !open ( iSpecs, file='params.spc',   status='unknown'     )

  ! Read in options from a file.
  !call QP%specs ( iSpecs, INFO )

  ! Set options.
  call qp%set   ( 'Print level        1', iPrint, iSumm, Errors )
  call qp%set   ( 'Summary level      1', iPrint, iSumm, Errors )
  call QP%set   ( 'Print frequency    1', iPrint, iSumm, Errors )
  call QP%set   ( 'Summary frequency  1', iPrint, iSumm, Errors )

end subroutine wsqic

subroutine sqicSolve (Obj)  bind ( C, name="sqicSolve" )
  use snModulePrecision, only : ip, rp
  use SQICModule

  implicit none

  integer                   :: INFO
  integer(ip)               :: nS, n_inf

  real(rp)                  :: Obj, sInf


  ! Solve the QP.
  call QP%solve ( 'Cold', INFO, nS, n_inf, sInf, Obj )

end subroutine sqicSolve


subroutine sqicSolveStabilized (Obj,mu,lenpi,piE)  bind ( C, name="sqicSolveStabilized" )
  use snModulePrecision, only : ip, rp
  use SQICModule

  implicit none

  integer                   :: INFO
  integer(ip)               :: nS, n_inf, lenpi

  real(rp)                  :: Obj, sInf, mu, piE(lenpi)

  ! Solve the QP.
  call QP%solveR ( 'Cold', mu, lenpi, piE, INFO, nS, n_inf, sInf, Obj )

end subroutine sqicSolveStabilized


subroutine sqicDestroy ()  bind ( C, name="sqicDestroy" )
  use SQICModule

  implicit none

  call QP%end

  if ( allocated(cObj) )   deallocate ( cObj )
  if ( allocated(Names) )  deallocate ( Names )

  !close       ( iPrint )

end subroutine sqicDestroy
