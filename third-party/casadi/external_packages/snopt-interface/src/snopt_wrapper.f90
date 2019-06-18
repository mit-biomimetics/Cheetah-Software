!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File  c_snopt.f90 --- C/C++ wrapper for SNOPT7
!
! 2014 Jul 07: First version.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module snopt_wrapper
  use  iso_c_binding
  implicit none

  public

  external :: sninit, snspec, snopta, snkera, snmema, snjac, &
              snoptb, snkerb, snoptc, snkerc, snmem,         &
              npopt,  npkern,                                &
              sngeti, sngetr, sngetc, snseti, snsetr, snset, &
              snFileOpenRead, snFileOpenAppend, snFileClose

  public  :: f_sninit, f_snsetprint, f_snspec, &
             f_snmema, f_snopta, f_snkera, f_snjac,  &
             f_snmem,  f_snoptb, f_snkerb, f_snoptc, f_snkerc, &
             f_snset,  f_snseti, f_snsetr, &
             f_sngetc, f_sngeti, f_sngetr, &
             f_snend
  private :: newunit

  !-----------------------------------------------------------------------------

  interface
     ! Interface for user-defined subroutines.
     subroutine snLog(iAbort, info, HQNType, KTcond, MjrPrt, minimz,    &
                      n, nb, nnCon0, nS, itn, nMajor, nMinor, nSwap,    &
                      condHz, iObj, sclObj, ObjAdd, fMrt, PenNrm, step, &
                      prInf, duInf, vimax, virel, hs,                   &
                      ne, nlocJ, locJ, indJ, Jcol,                      &
                      Ascale, bl, bu, fCon, yCon, x,                    &
                      cu, lencu, iu, leniu, ru, lenru,                  &
                      cw, lencw, iw, leniw, rw, lenrw)
       logical, intent(in) :: KTcond(2)
       integer, intent(in) :: iObj, info(6), HQNType,                   &
                              lencu, lencw, leniu, leniw, lenru, lenrw, &
                              MjrPrt, minimz, n, ne, nb, nlocJ, nnCon0, &
                              nS, itn, nMajor, nMinor, nSwap,           &
                              hs(nb), locJ(nlocJ), indJ(ne)

       double precision, intent(in) :: condHz, sclObj, ObjAdd, fMrt,             &
                                       PenNrm, virel, vimax, step, prInf, duInf, &
                                       Ascale(nb), bl(nb), bu(nb), fCon(nnCon0), &
                                       Jcol(ne), yCon(nnCon0), x(nb)

       integer,          intent(inout) :: iu(leniu), iw(leniw)
       double precision, intent(inout) :: ru(lenru), rw(lenrw)
       character,        intent(inout) :: cu(lencu)*8, cw(lencw)*8

       integer,          intent(out)   :: iAbort
     end subroutine snLog

     !--------------------------------------------------------------------------

     subroutine snLog2(Prob, ProbTag, Elastc, gotR, jstFea, feasbl, &
                       m, mBS, nnH, nS, jSq, jBr, jSr,              &
                       linesP, linesS, itn, itQP, kPrc, lvlInf,     &
                       pivot, step, nInf, sInf, wtInf,              &
                       ObjPrt, condHz, djqPrt, rgNorm, kBS, xBS,    &
                       iw, leniw)
       character, intent(in)    :: ProbTag*20
       logical,   intent(in)    :: Elastc, gotR, jstFea, feasbl
       integer,   intent(in)    :: Prob, m, mBS, nnH, nS, jSq, jBr, jSr,    &
                                   itn, itQP, kPrc, linesP, linesS, lvlInf, &
                                   nInf, kBS(mBS), leniw
       integer,   intent(inout) :: iw(leniw)

       double precision, intent(in) :: condHz, djqPrt, ObjPrt, pivot, rgNorm, &
                                       step, sInf, wtInf, xBS(mBS)
     end subroutine snLog2

     !--------------------------------------------------------------------------

     subroutine snSTOP(iAbort, info, HQNType, KTcond, MjrPrt, minimz,    &
                       m, maxS, n, nb, nnCon0, nnCon, nnObj0, nnObj, nS, &
                       itn, nMajor, nMinor, nSwap,                       &
                       condHz, iObj, sclObj, ObjAdd, fMrt, PenNrm, step, &
                       prInf, duInf, vimax, virel, hs,                   &
                       ne, nlocJ, locJ, indJ, Jcol, negCon,              &
                       Ascale, bl, bu, fCon, gCon, gObj,                 &
                       yCon, pi, rc, rg, x,             &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw)
       logical, intent(in) :: KTcond(2)
       integer, intent(in) :: HQNType, info(6), iObj, itn,                &
                              lencu, lencw, leniu, leniw, lenru, lenrw,   &
                              MjrPrt, minimz, m, maxS, n, nb, ne, negCon, &
                              nlocJ, nnCon0, nnCon, nnObj0, nnObj,        &
                              nMajor, nMinor, nS, nSwap,                  &
                              hs(nb), locJ(nlocJ), indJ(ne)
       double precision, intent(in) :: &
            condHz, sclObj, ObjAdd, fMrt, PenNrm, virel, vimax, step, &
            prInf, duInf, Ascale(nb), bl(nb), bu(nb), fCon(nnCon0),   &
            gCon(negCon), gObj(nnObj0), Jcol(ne), pi(m),              &
            rc(nb), rg(maxS), yCon(nnCon0), x(nb)

       integer,          intent(inout) :: iu(leniu), iw(leniw)
       double precision, intent(inout) :: ru(lenru), rw(lenrw)
       character,        intent(inout) :: cu(lencu)*8, cw(lencw)*8

       integer, intent(out) :: iAbort

     end subroutine snSTOP

     !--------------------------------------------------------------------------

     subroutine sqLog(Prob, ProbTag, Elastc, gotR, jstFea, feasbl, &
                      m, mBS, nnH, nS, jSq, jBr, jSr,              &
                      linesP, linesS, itn, itQP, kPrc, lvlInf,     &
                      pivot, step, nInf, sInf, wtInf,              &
                      ObjPrt, condHz, djqPrt, rgNorm, kBS, xBS,    &
                      iw, leniw)
       character,        intent(in) :: ProbTag*20
       logical,          intent(in) :: Elastc, gotR, jstFea, feasbl
       integer,          intent(in) :: Prob, m, mBS, nnH, nS, jSq, jBr, jSr, &
                                       itn, itQP, kPrc, linesP, linesS,      &
                                       lvlInf, nInf, kBS(mBS), leniw
       double precision, intent(in) :: condHz, djqPrt, ObjPrt, pivot, rgNorm, &
                                       step, sInf, wtInf, xBS(mBS)
       integer,       intent(inout) :: iw(leniw)
     end subroutine sqLog

     !--------------------------------------------------------------------------
     ! SNOPTA:
     subroutine iusrfunA(Status, n, x, needf, nF, F, needG, lenG, G, &
                         cu, lencu, iu, leniu, ru, lenru)
       integer,          intent(in) :: n, nF, needF, needG, lenG, &
                                       lencu, leniu, lenru
       double precision, intent(in) :: x(n)

       integer,          intent(inout) :: Status, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       character,        intent(inout) :: cu(lencu)*8

       double precision, intent(out)   :: F(nF), G(lenG)

     end subroutine iusrfunA

     !--------------------------------------------------------------------------
     ! SNOPTB:
     subroutine ifuncon(mode, nnCon, nnJac, negCon, &
                        x, fCon, gCon, Status, &
                        cu, lencu, iu, leniu, ru, lenru)
       integer,          intent(in)    :: mode, nnCon, nnJac, negCon, &
                                          lencu, leniu, lenru
       double precision, intent(in)    :: x(nnJac)

       character,        intent(inout) :: cu(lencu)*8
       integer,          intent(inout) :: Status, iu(leniu)
       double precision, intent(inout) :: ru(lenru)

       double precision, intent(out)   :: fCon(nnCon), gCon(negCon)

     end subroutine ifuncon

     subroutine ifunobj(mode, nnObj, x, fObj, gObj, Status, &
                        cu, lencu, iu, leniu, ru, lenru)
       integer,          intent(in)    :: mode, nnObj, lencu, leniu, lenru
       double precision, intent(in)    :: x(nnObj)

       character,        intent(inout) :: cu(lencu)*8
       integer,          intent(inout) :: Status, iu(leniu)
       double precision, intent(inout) :: ru(lenru)

       double precision, intent(out)   :: fObj, gObj(nnObj)

     end subroutine ifunobj

     !--------------------------------------------------------------------------
     ! SNOPTC:
     subroutine iusrfunC(mode, nnObj, nnCon, nnJac, nnL, negCon, &
                         x, fObj, gObj, fCon, gCon, Status, &
                         cu, lencu, iu, leniu, ru, lenru)
       integer,          intent(in)    :: mode, nnObj, nnCon, nnJac, nnL, &
                                          negCon, lencu, leniu, lenru
       double precision, intent(in)    :: x(nnL)

       integer,          intent(inout) :: Status, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       character,        intent(inout) :: cu(lencu)*8

       double precision, intent(out)   :: fObj, fCon(nnCon), &
                                          gObj(nnObj), gCon(negCon)
     end subroutine iusrfunC
  end interface

  !-----------------------------------------------------------------------------

  ! Character arrays don't work well with Fortran/C/C++ so have a dummy one here.
  integer, parameter :: lencw = 500
  character*8        :: cw(lencw)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sninit(name, len, summOn, iw, leniw, rw, lenrw) &
       bind(C,name="f_sninit")

    integer(c_int),    intent(in), value :: len, summOn, leniw, lenrw
    character(c_char), intent(in)        :: name(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    ! Call snInit.  If a name is provided, use it as the print file name.
    ! Else assume no print file required (for now).
    !
    ! 07 Jul 2014: First version.
    !===========================================================================
    character(len) :: file
    integer        :: j, iPrt, iSum

    if (len == 0) then
       iPrt = 0
    else
       if (name(1) == c_null_char) then
          iPrt = 0
       else
          file  = ''
          do j = 1, len
             if (name(j) == c_null_char) exit
             file(j:j) = name(j)
          end do
          iPrt = newunit()
          call snFileOpenAppend(iPrt, trim(file))
       end if
    end if

    if (summOn == 0) then
       iSum = 0
    else
       iSum = 6
    end if

    call snInit(iPrt, iSum, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sninit

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snsetprint(name, len, iw, leniw, rw, lenrw) &
       bind(C,name="f_snsetprint")

    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: name(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    ! Set print file name and unit.
    !===========================================================================
    integer        :: Errors, j, iPrt
    character(len) :: prtfile

    prtfile = ''
    do j = 1, len
       if (name(j) == c_null_char) exit
       prtfile(j:j) = name(j)
    end do

    if (prtfile /= '') then
       iPrt = newunit()
       call snFileOpenAppend(iPrt,trim(prtfile))
       call snSeti('Print file', iPrt, 0, 0, Errors, &
                    cw, lencw, iw, leniw, rw, lenrw)
    end if

  end subroutine f_snsetprint

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snspec(name, len, inform, iw, leniw, rw, lenrw) &
       bind(C,name="f_snspec")

    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: name(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: inform

    !===========================================================================
    ! Read options from the given specifications file.
    !===========================================================================
    integer        :: iSpec, j
    character(len) :: spcfile

    inform  = 0
    iSpec   = 4

    ! Get specs file name.
    spcfile = ''
    do j = 1, len
       if (name(j) == c_null_char) exit
       spcfile(j:j) = name(j)
    end do

    ! If we have a file, try to read it.
    if (spcfile /= '') then
       call snFileOpenRead(iSpec,trim(spcfile))
       call snSpec(iSpec, inform, cw, lencw, iw, leniw, rw, lenrw)
       close(iSpec)
    end if

  end subroutine f_snspec

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snmema(inform, nF, n, lenA, lenG, miniw, minrw, &
                      iw, leniw, rw, lenrw) bind(C,name="f_snmema")

    integer(c_int), intent(in), value :: n, nF, lenA, lenG, leniw, lenrw
    integer(c_int), intent(inout)     :: iw(leniw)
    real(c_double), intent(inout)     :: rw(lenrw)
    integer(c_int), intent(out)       :: inform, miniw, minrw

    !===========================================================================
    ! Call snMemA to get workspace estimate.
    !===========================================================================
    integer :: mincw, nxname, nFname

    nxname = 1
    nFname = 1

    call snMemA(inform, nF, n, nxname, nfname, &
                lenA, lenG, mincw, miniw, minrw, &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snmema

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snjac(inform, nF, n, usrfunC, x, xlow, xupp, &
                     iAfun, jAvar, lenA, neA, A,            &
                     iGfun, jGvar, lenG, neG,               &
                     miniw, minrw, iu, leniu, ru, lenru,    &
                     iw, leniw, rw, lenrw) bind(C,name="f_snjac")

    integer(c_int), intent(in), value :: n, nF, lenA, lenG, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in)        :: xlow(n), xupp(n)

    integer(c_int), intent(inout)     :: iAfun(lenA), jAvar(lenA), &
                                         iGfun(lenG), jGvar(lenG), &
                                         iw(leniw), iu(leniu)
    real(c_double), intent(inout)     :: x(n), A(lenA), rw(lenrw), ru(lenru)

    integer(c_int), intent(out)       :: inform, neA, neG, miniw, minrw
    type(c_funptr), value             :: usrfunc

    !===========================================================================
    ! Call snJac to get Jacobian structure for SNOPTA.
    !===========================================================================
    integer :: mincw
    procedure(iusrfunA), pointer :: usrfun

    call c_f_procpointer(usrfunC, usrfun)

    call snJac(inform, nF, n, usrfun,              &
               iAfun, jAvar, lenA, neA, A,         &
               iGfun, jGvar, lenG, neG,            &
               x, xlow, xupp, mincw, miniw, minrw, &
               cw, lencw, iu, leniu, ru, lenru,    &
               cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snjac

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snopta(Start, Prob, nF, n, ObjAdd, ObjRow, c_usrfun, &
                      iAfun, jAvar, neA, A, iGfun, jGvar, neG,      &
                      xlow, xupp, Flow, Fupp, x, xstate, xmul,      &
                      F, Fstate, Fmul, INFO, nS, nInf, sInf,        &
                      miniw, minrw,                                 &
                      iu, leniu, ru, lenru,                         &
                      iw, leniw, rw, lenrw) bind(C,name="f_snopta")

    integer(c_int), intent(in), value :: Start, n, nF, ObjRow, neA, neG, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    integer(c_int),    intent(in) :: iAfun(neA), jAvar(neA), &
                                     iGfun(neG), jGvar(neG)
    real(c_double),    intent(in) :: Flow(nF), Fupp(nF), &
                                     xlow(n), xupp(n), A(neA)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, nS, xstate(n), Fstate(nF)
    real(c_double), intent(inout) :: sInf, F(nF), Fmul(nF), x(n), xmul(n)

    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    type(c_funptr), value         :: c_usrfun

    !===========================================================================
    ! Solve the problem with SNOPTA.
    !===========================================================================
    integer      :: j, nxname, nFname, mincw
    character(8) :: pname, xnames(1), Fnames(1)

    procedure(iusrfunA), pointer :: usrfun

    nxname = 1
    nFname = 1

    call c_f_procpointer(c_usrfun, usrfun)

    pname  = ''
    do j = 1, 8
       if (Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call snoptA(Start, nF, n, nxname, nFname,     &
                ObjAdd, ObjRow, pname, usrfun,    &
                iAfun, jAvar, neA, neA, A,        &
                iGfun, jGvar, neG, neG,           &
                xlow, xupp, xnames,               &
                Flow, Fupp, Fnames,               &
                x, xstate, xmul, F, Fstate, Fmul, &
                INFO, mincw, miniw, minrw,        &
                nS, nInf, sInf,                   &
                cw, lencw, iu, leniu, ru, lenru,  &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snopta

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snkera(Start, Prob, nF, n, ObjAdd, ObjRow,             &
                      c_usrfun, c_snLog, c_snLog2, c_sqLog, c_snSTOP, &
                      iAfun, jAvar, neA, A, iGfun, jGvar, neG,        &
                      xlow, xupp, Flow, Fupp, x, xstate, xmul,        &
                      F, Fstate, Fmul, INFO, nS, nInf, sInf,          &
                      miniw, minrw,                                   &
                      iu, leniu, ru, lenru,                           &
                      iw, leniw, rw, lenrw) bind(C,name="f_snkera")

    integer(c_int), intent(in), value :: Start, n, nF, ObjRow, neA, neG, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    integer(c_int),    intent(in) :: iAfun(neA), jAvar(neA), &
                                     iGfun(neG), jGvar(neG)
    real(c_double),    intent(in) :: Flow(nF), Fupp(nF), &
                                     xlow(n), xupp(n), A(neA)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, nS, xstate(n), Fstate(nF)
    real(c_double), intent(inout) :: sInf, F(nF), Fmul(nF), x(n), xmul(n)

    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    type(c_funptr), value         :: c_usrfun, c_snLog, c_snLog2, &
                                     c_sqLog, c_snSTOP

    !===========================================================================
    ! Advanced solve with SNOPTA.
    !===========================================================================
    integer      :: j, nxname, nFname, mincw
    character(8) :: pname, xnames(1), Fnames(1)

    procedure(iusrfunA), pointer :: usrfun
    procedure(snSTOP),   pointer :: mySTOP
    procedure(snLog),    pointer :: myLog
    procedure(snLog2),   pointer :: myLog2
    procedure(sqLog),    pointer :: myLogQ

    nxname = 1
    nFname = 1

    usrfun => null()
    mySTOP => null()
    myLog  => null()
    myLog2 => null()
    myLogQ => null()

    call c_f_procpointer(c_usrfun, usrfun)
    call c_f_procpointer(c_snLog,  myLog )
    call c_f_procpointer(c_snLog2, myLog2)
    call c_f_procpointer(c_sqLog,  myLogQ)
    call c_f_procpointer(c_snSTOP, mySTOP)

    if (.not. associated(myLog) ) myLog  => snLog
    if (.not. associated(myLog2)) myLog2 => snLog2
    if (.not. associated(myLogQ)) myLogQ => sqLog
    if (.not. associated(mySTOP)) mySTOP => snSTOP


    pname  = ''
    do j = 1, 8
       if (Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call snKerA(Start, nF, n, nxname, nFname,          &
                ObjAdd, ObjRow, pname,                 &
                usrfun, myLog, myLog2, myLogQ, mySTOP, &
                iAfun, jAvar, neA, neA, A,             &
                iGfun, jGvar, neG, neG,                &
                xlow, xupp, xnames,                    &
                Flow, Fupp, Fnames,                    &
                x, xstate, xmul, F, Fstate, Fmul,      &
                INFO, mincw, miniw, minrw,             &
                nS, nInf, sInf,                        &
                cw, lencw, iu, leniu, ru, lenru,       &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snkera

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snmem(info, m, n, ne, negCon, nnCon, nnObj, nnJac, &
                     miniw, minrw, iw, leniw, rw, lenrw) &
                     bind(C,name="f_snmem")

    integer(c_int), intent(in), value :: m, n, ne, negCon, nnCon, nnObj, nnJac, &
                                         leniw, lenrw
    integer(c_int), intent(inout)     :: iw(leniw)
    real(c_double), intent(inout)     :: rw(lenrw)

    integer(c_int), intent(out)       :: INFO, miniw, minrw

    !===========================================================================
    ! Estimate workspace for SNOPTB/C.
    !===========================================================================
    integer :: mincw

    call snMem(INFO, m, n, ne, negCon, nnCon, nnJac, nnObj, &
               mincw, miniw, minrw, &
               cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snmem

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snoptb(Start, Prob, m, n, ne, nnCon, nnObj, nnJac, &
                      iObj, ObjAdd, c_funcon, c_funobj,           &
                      Jval, indJ, locJ, bl, bu, hs,               &
                      x, pi, rc, INFO, nS, nInf, sInf, Obj,       &
                      miniw, minrw, iu, leniu, ru, lenru,         &
                      iw, leniw, rw, lenrw) &
                      bind(C,name="f_snoptb")

    integer(c_int), intent(in), value :: Start, m, n, ne, &
                                         nnCon, nnObj, nnJac, iObj, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    real(c_double),    intent(in) :: Jval(ne), bl(n+m), bu(n+m)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, nS, indJ(ne), locJ(n+1), &
                                     hs(n+m)
    real(c_double), intent(inout) :: sInf, x(n+m), pi(m), rc(n+m)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: Obj
    type(c_funptr), value         :: c_funcon, c_funobj

    !===========================================================================
    ! Solve the problem with SNOPTB.
    !===========================================================================
    integer      :: j, nName, mincw
    character(8) :: pname, Names(1)
    character(4) :: nStart

    procedure(ifuncon), pointer :: funcon
    procedure(ifunobj), pointer :: funobj

    nName = 1

    if      (Start == 1) then
       nStart = 'Warm'
    else if (Start == 2) then
       nStart = 'Hot'
    else
       nStart = 'Cold'
    end if

    call c_f_procpointer(c_funcon, funcon)
    call c_f_procpointer(c_funobj, funobj)

    pname  = ''
    do j = 1, 8
       if (Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call snOptB(nStart, m, n, ne, nName,         &
                nnCon, nnObj, nnJac,             &
                iObj, ObjAdd, pname,             &
                funcon, funobj,                  &
                Jval, indJ, locJ, bl, bu, Names, &
                hs, x, pi, rc,                   &
                INFO, mincw, miniw, minrw,       &
                nS, nInf, sInf, Obj,             &
                cw, lencw, iu, leniu, ru, lenru, &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snoptb

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snkerb(Start, Prob, m, n, ne, nnCon, nnObj, nnJac, &
                      iObj, ObjAdd, c_funcon, c_funobj,           &
                      c_snLog, c_snLog2, c_sqLog, c_snSTOP,       &
                      Jval, indJ, locJ, bl, bu, hs,               &
                      x, pi, rc, INFO, nS, nInf, sInf, Obj,       &
                      miniw, minrw, iu, leniu, ru, lenru,         &
                      iw, leniw, rw, lenrw) bind(C,name="f_snkerb")

    integer(c_int), intent(in), value :: Start, m, n, ne, &
                                         nnCon, nnObj, nnJac, iObj, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    real(c_double),    intent(in) :: Jval(ne), bl(n+m), bu(n+m)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, nS, indJ(ne), locJ(n+1), &
                                     hs(n+m)
    real(c_double), intent(inout) :: sInf, x(n+m), pi(m), rc(n+m)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: Obj
    type(c_funptr), value         :: c_funcon, c_funobj
    type(c_funptr), value         :: c_snLog, c_snLog2, c_sqLog, c_snSTOP

    !===========================================================================
    ! Advanced solve with SNOPTB.
    !===========================================================================
    integer      :: j, nName, mincw
    character(8) :: pname, Names(1)
    character(4) :: nStart

    procedure(ifuncon), pointer :: funcon
    procedure(ifunobj), pointer :: funobj
    procedure(snSTOP),  pointer :: mySTOP
    procedure(snLog),   pointer :: myLog
    procedure(snLog2),  pointer :: myLog2
    procedure(sqLog),   pointer :: myLogQ

    nName = 1

    if      (Start == 1) then
       nStart = 'Warm'
    else if (Start == 2) then
       nStart = 'Hot'
    else
       nStart = 'Cold'
    end if

    funcon => null()
    funobj => null()
    mySTOP => null()
    myLog  => null()
    myLog2 => null()
    myLogQ => null()

    call c_f_procpointer(c_funcon, funcon)
    call c_f_procpointer(c_funobj, funobj)
    call c_f_procpointer(c_snLog,  myLog )
    call c_f_procpointer(c_snLog2, myLog2)
    call c_f_procpointer(c_sqLog,  myLogQ)
    call c_f_procpointer(c_snSTOP, mySTOP)

    if (.not. associated(myLog) ) myLog  => snLog
    if (.not. associated(myLog2)) myLog2 => snLog2
    if (.not. associated(myLogQ)) myLogQ => sqLog
    if (.not. associated(mySTOP)) mySTOP => snSTOP


    pname  = ''
    do j = 1, 8
       if (Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call snKerB(nStart, m, n, ne, nName,         &
                nnCon, nnObj, nnJac,             &
                iObj, ObjAdd, pname,             &
                funcon, funobj,                  &
                myLog, myLog2, myLogQ, mySTOP,   &
                Jval, indJ, locJ, bl, bu, Names, &
                hs, x, pi, rc,                   &
                INFO, mincw, miniw, minrw,       &
                nS, nInf, sInf, Obj,             &
                cw, lencw, iu, leniu, ru, lenru, &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snkerb

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snoptc(Start, Prob, m, n, ne, nnCon, nnObj, nnJac,       &
                      iObj, ObjAdd, c_usrfun, Jval, indJ, locJ,         &
                      bl, bu, hs, x, pi, rc, INFO, nS, nInf, sInf, Obj, &
                      miniw, minrw, iu, leniu, ru, lenru,               &
                      iw, leniw, rw, lenrw) bind(C,name="f_snoptc")

    integer(c_int), intent(in), value :: Start, m, n, ne, &
                                         nnCon, nnObj, nnJac, iObj, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    real(c_double),    intent(in) :: Jval(ne), bl(n+m), bu(n+m)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, nS, indJ(ne), locJ(n+1), &
                                     hs(n+m)
    real(c_double), intent(inout) :: sInf, x(n+m), pi(m), rc(n+m)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: Obj
    type(c_funptr), value         :: c_usrfun

    !===========================================================================
    ! Solve the problem with SNOPTC.
    !===========================================================================
    integer      :: j, nName, mincw
    character(8) :: pname, Names(1)
    character(4) :: nStart

    procedure(iusrfunC), pointer :: usrfun

    nName = 1

    if      (Start == 1) then
       nStart = 'Warm'
    else if (Start == 2) then
       nStart = 'Hot'
    else
       nStart = 'Cold'
    end if

    call c_f_procpointer (c_usrfun, usrfun)

    pname  = ''
    do j = 1, 8
       if (Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call snOptC(nStart, m, n, ne, nName,         &
                nnCon, nnObj, nnJac,             &
                iObj, ObjAdd, pname, usrfun,     &
                Jval, indJ, locJ, bl, bu, Names, &
                hs, x, pi, rc,                   &
                INFO, mincw, miniw, minrw,       &
                nS, nInf, sInf, Obj,             &
                cw, lencw, iu, leniu, ru, lenru, &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snoptc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snkerc(Start, Prob, m, n, ne,                          &
                      nnCon, nnObj, nnJac, iObj, ObjAdd,              &
                      c_usrfun, c_snLog, c_snLog2, c_sqLog, c_snSTOP, &
                      Jval, indJ, locJ, bl, bu, hs,                   &
                      x, pi, rc, INFO, nS, nInf, sInf, Obj,           &
                      miniw, minrw, iu, leniu, ru, lenru,             &
                      iw, leniw, rw, lenrw) bind(C,name="f_snkerc")

    integer(c_int), intent(in), value :: Start, m, n, ne, &
                                         nnCon, nnObj, nnJac, iObj, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    real(c_double),    intent(in) :: Jval(ne), bl(n+m), bu(n+m)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, nS, indJ(ne), locJ(n+1), &
                                     hs(n+m)
    real(c_double), intent(inout) :: sInf, x(n+m), pi(m), rc(n+m)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: Obj
    type(c_funptr), value         :: c_usrfun
    type(c_funptr), value         :: c_snLog, c_snLog2, c_sqLog, c_snSTOP

    !===========================================================================
    ! Advanced solve with SNOPTC.
    !===========================================================================
    integer      :: j, nName, mincw
    character(8) :: pname, Names(1)
    character(4) :: nStart

    procedure(iusrfunC), pointer :: usrfun
    procedure(snSTOP),   pointer :: mySTOP
    procedure(snLog),    pointer :: myLog
    procedure(snLog2),   pointer :: myLog2
    procedure(sqLog),    pointer :: myLogQ

    nName = 1

    if      (Start == 1) then
       nStart = 'Warm'
    else if (Start == 2) then
       nStart = 'Hot'
    else
       nStart = 'Cold'
    end if

    usrfun => null()
    mySTOP => null()
    myLog  => null()
    myLog2 => null()
    myLogQ => null()

    call c_f_procpointer(c_usrfun, usrfun)
    call c_f_procpointer(c_snLog,  myLog )
    call c_f_procpointer(c_snLog2, myLog2)
    call c_f_procpointer(c_sqLog,  myLogQ)
    call c_f_procpointer(c_snSTOP, mySTOP)

    if (.not. associated(myLog) ) myLog  => snLog
    if (.not. associated(myLog2)) myLog2 => snLog2
    if (.not. associated(myLogQ)) myLogQ => sqLog
    if (.not. associated(mySTOP)) mySTOP => snSTOP


    pname  = ''
    do j = 1, 8
       if (Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call snKerC(nStart, m, n, ne, nName,               &
                nnCon, nnObj, nnJac,                   &
                iObj, ObjAdd, pname,                   &
                usrfun, myLog, myLog2, myLogQ, mySTOP, &
                Jval, indJ, locJ, bl, bu, Names,       &
                hs, x, pi, rc,                         &
                INFO, mincw, miniw, minrw,             &
                nS, nInf, sInf, Obj,                   &
                cw, lencw, iu, leniu, ru, lenru,       &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snkerc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_npopt(n, nclin, ncnln, ldA, ldcJ, ldH,                     &
                     A, bl, bu, c_funcon, c_funobj, INFO, majIts, iState, &
                     cCon, cJac, cMul, Objf, grad, Hess, x,               &
                     iw, leniw, rw, lenrw) bind(C,name="f_npopt")

    integer(c_int), intent(in), value :: n, nclin, ncnln, ldA, ldcJ, ldH, &
                                         leniw, lenrw
    real(c_double), intent(in)    :: bl(n+nclin+ncnln), bu(n+nclin+ncnln), &
                                     A(ldA,*)

    integer(c_int), intent(inout) :: iState(n+nclin+ncnln)
    real(c_double), intent(inout) :: cCon(*), grad(n), Hess(ldH,*), x(n), &
                                     cJac(ldcJ,*), cMul(n+nclin+ncnln)
    integer(c_int), intent(inout) :: iw(leniw)
    real(c_double), intent(inout) :: rw(lenrw)

    integer(c_int), intent(out)   :: INFO, majIts
    real(c_double), intent(out)   :: Objf
    type(c_funptr), value         :: c_funcon, c_funobj

    !===========================================================================
    ! Call SNOPT via NPOPT interface.
    !===========================================================================
    procedure(ifuncon), pointer :: funcon
    procedure(ifunobj), pointer :: funobj

    call c_f_procpointer (c_funcon, funcon)
    call c_f_procpointer (c_funobj, funobj)

    call npopt(n, nclin, ncnln, ldA, ldcJ, ldH,       &
               A, bl, bu, funcon, funobj,             &
               INFO, majIts, iState,                  &
               cCon, cJac, cMul, Objf, grad, Hess, x, &
               iw, leniw, rw, lenrw)

  end subroutine f_npopt

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_npkern(n, nclin, ncnln, ldA, ldcJ, ldH, A, bl, bu,      &
                      c_funcon, c_funobj, c_snLog, c_snLog2, c_snSTOP, &
                      INFO, majIts, iState,                            &
                      cCon, cJac, cMul, Objf, grad, Hess, x,           &
                      iw, leniw, rw, lenrw) bind(C,name="f_npkern")

    integer(c_int), intent(in), value :: n, nclin, ncnln, ldA, ldcJ, ldH, &
                                         leniw, lenrw

    real(c_double), intent(in)    :: bl(n+nclin+ncnln), bu(n+nclin+ncnln), &
                                     A(ldA,*)
    integer(c_int), intent(inout) :: iState(n+nclin+ncnln)
    real(c_double), intent(inout) :: cCon(*), grad(n), Hess(ldH,*), x(n), &
                                     cJac(ldcJ,*), cMul(n+nclin+ncnln)
    integer(c_int), intent(inout) :: iw(leniw)
    real(c_double), intent(inout) :: rw(lenrw)

    integer(c_int), intent(out)   :: INFO, majIts
    real(c_double), intent(out)   :: Objf
    type(c_funptr), value         :: c_funcon, c_funobj
    type(c_funptr), value         :: c_snLog, c_snLog2, c_snSTOP

    !===========================================================================
    ! Advanced solve with NPOPT interface.
    !===========================================================================
    procedure(ifuncon),  pointer :: funcon
    procedure(ifunobj),  pointer :: funobj
    procedure(snSTOP),   pointer :: mySTOP
    procedure(snLog),    pointer :: myLog
    procedure(snLog2),   pointer :: myLog2

    funcon => null()
    funobj => null()
    mySTOP => null()
    myLog  => null()
    myLog2 => null()

    call c_f_procpointer(c_funcon, funcon)
    call c_f_procpointer(c_funobj, funobj)
    call c_f_procpointer(c_snLog,  myLog )
    call c_f_procpointer(c_snLog2, myLog2)
    call c_f_procpointer(c_snSTOP, mySTOP)

    if (.not. associated(myLog) ) myLog  => snLog
    if (.not. associated(myLog2)) myLog2 => snLog2
    if (.not. associated(mySTOP)) mySTOP => snSTOP

    call npKerN(n, nclin, ncnln, ldA, ldcJ, ldH,       &
                A, bl, bu,                             &
                funcon, funobj, myLog, myLog2, mySTOP, &
                INFO, majIts, iState,                  &
                cCon, cJac, cMul, Objf, grad, Hess, x, &
                iw, leniw, rw, lenrw)

  end subroutine f_npkern

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snend(iw, leniw, rw, lenrw) bind(C,name="f_snend")
    integer(c_int),    intent(in), value :: leniw, lenrw
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    ! Finish up.
    !===========================================================================
    close(iw(12))  ! print file
    if (iw(13) /= 6) close(iw(13))  ! summary file

  end subroutine f_snend

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snset(option, len, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_snset")
    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    ! Set option via string.
    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''
    do j = 1, len
       if (option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call snSet(buffer, 0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snset

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snseti(option, len, ivalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_snseti")
    integer(c_int),    intent(in), value :: len, ivalue, leniw, lenrw
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    ! Set option with integer value.
    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''
    do j = 1, len
       if (option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call snSeti(buffer, ivalue, 0, 0, Errors, &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snseti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snsetr(option, len, rvalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_snsetr")
    integer(c_int),    intent(in), value :: len, leniw, lenrw
    real(c_double),    intent(in), value :: rvalue
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    ! Set option with real value.
    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''

    do j = 1, len
       if (option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call snSetr(buffer, rvalue, 0, 0, Errors, &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snsetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sngetc(option, lin, cvalue, lout, Errors, &
                      iw, leniw, rw, lenrw) bind(C,name="f_sngetc")
    integer(c_int),    intent(in), value :: lin, lout, leniw, lenrw
    character(c_char), intent(in)        :: option(lin)
    character(c_char), intent(inout)     :: cvalue(lout)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    ! Get option value via string.
    !===========================================================================
    character(lin)  :: buffer
    character(lout) :: buffout
    integer         :: j

    errors = 0
    buffer = ''
    do j = 1, lin
       if (option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call snGetC(buffer, buffout, Errors, cw, lencw, iw, leniw, rw, lenrw)

    do j = 1, lout-1
       cvalue(j) = buffout(j:j)
    end do
    cvalue(lout) = c_null_char

  end subroutine f_sngetc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sngeti(option, len, ivalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_sngeti")
    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: ivalue, Errors

    !===========================================================================
    ! Get integer option value.
    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''
    do j = 1, len
       if (option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call snGetI(buffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sngeti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sngetr(option, len, rvalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_sngetr")
    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: option(len)
    real(c_double),    intent(out)       :: rvalue
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    ! Get real option value.
    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''
    do j = 1, len
       if (option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call snGetR(buffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sngetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function newunit()
    !===========================================================================
    integer, parameter :: unit_min = 10, unit_max = 1000
    logical :: opened
    integer :: j

    newunit = -1
    do j = unit_min, unit_max
       inquire(unit=j,opened=opened)
       if (.not. opened) then
          newunit = j
          exit
       end if
    end do

  end function newunit

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module snopt_wrapper
