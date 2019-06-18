!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File  sqopt_wrapper.f90 --- C/C++ wrapper for SQOPT
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module sqopt_wrapper
  use  iso_c_binding
  implicit none

  public

  external :: sqinit, sqspec,                                &
              sqmem,  snoptq, snkerq,                        &
              sqgeti, sqgetr, sqgetc, sqseti, sqsetr, sqset, &
              snFileOpenRead, snFileOpenAppend, snFileClose

  public  :: f_sqinit, f_sqsetprint, f_sqspec, &
             f_sqmem,  f_sqopt,  f_snkerq, &
             f_sqset,  f_sqseti, f_sqsetr, &
             f_sqgetc, f_sqgeti, f_sqgetr, &
             f_sqend
  private :: newunit

  !-----------------------------------------------------------------------------

  interface
     ! Interface for user-defined subroutines.

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

     subroutine iusrHx(nnH, x, Hx, nState, cu, lencu, iu, leniu, ru, lenru)
       integer,          intent(in)    :: nnH, nState, lencu, leniu, lenru
       double precision, intent(in)    :: x(nnH)

       character,        intent(inout) :: cu(lencu)*8
       integer,          intent(inout) :: iu(leniu)
       double precision, intent(inout) :: ru(lenru)

       double precision, intent(out)   :: Hx(nnH)
     end subroutine iusrHx

  end interface

  !-----------------------------------------------------------------------------

  ! Character arrays don't work well with Fortran/C/C++ so have a dummy one here.
  integer, parameter :: lencw = 500
  character*8        :: cw(lencw)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqinit(name, len, summOn, iw, leniw, rw, lenrw) &
       bind(C,name="f_sqinit")

    integer(c_int),    intent(in), value :: len, summOn, leniw, lenrw
    character(c_char), intent(in)        :: name(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    ! Call sqInit.  If a name is provided, use it as the print file name.
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

    call sqInit(iPrt, iSum, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sqinit

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqsetprint(name, len, iw, leniw, rw, lenrw) &
       bind(C,name="f_sqsetprint")

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
       call sqSeti('Print file', iPrt, 0, 0, Errors, &
                    cw, lencw, iw, leniw, rw, lenrw)
    end if

  end subroutine f_sqsetprint

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqspec(name, len, inform, iw, leniw, rw, lenrw) &
       bind(C,name="f_sqspec")

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
       call sqSpec(iSpec, inform, cw, lencw, iw, leniw, rw, lenrw)
       close(iSpec)
    end if

  end subroutine f_sqspec

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqmem(info, m, n, neA, ncObj, nnH, &
                     miniw, minrw, iw, leniw, rw, lenrw) &
                     bind(C,name="f_sqmem")

    integer(c_int), intent(in), value :: m, n, neA, ncObj, nnH, &
                                         leniw, lenrw
    integer(c_int), intent(inout)     :: iw(leniw)
    real(c_double), intent(inout)     :: rw(lenrw)

    integer(c_int), intent(out)       :: INFO, miniw, minrw

    !===========================================================================
    ! Estimate workspace for SQOPT.
    !===========================================================================
    integer :: mincw

    call sqMem(INFO, m, n, neA, ncObj, nnH, &
               mincw, miniw, minrw, &
               cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sqmem

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqopt(Start, c_qpHx, m, n, neA, ncObj, nnH, &
                     iObj, ObjAdd, Prob,                   &
                     valA, indA, locA, bl, bu, cObj,       &
                     eType, hs, x, pi, rc,                 &
                     INFO, nS, nInf, sInf, Obj,            &
                     miniw, minrw,                         &
                     iu, leniu, ru, lenru,                 &
                     iw, leniw, rw, lenrw) bind(C,name='f_sqopt')
    integer(c_int), intent(in), value :: Start, m, n, iObj, neA, ncObj, nnH, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    integer(c_int),    intent(in) :: eType(n+m), indA(neA), locA(n+1)
    real(c_double),    intent(in) :: cObj(ncObj), bl(n+m), bu(n+m), valA(neA)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, nS, hs(n)
    real(c_double), intent(inout) :: sInf, x(n), pi(m), rc(n+m)

    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: Obj
    type(c_funptr), value         :: c_qpHx

    !===========================================================================
    ! Solve the problem with SQOPT.
    !===========================================================================
    integer      :: j, nNames, mincw
    character(8) :: pname, Names(1), Fnames(1)
    character(4) :: nStart

    procedure(iusrHx), pointer :: qpHx

    nNames = 1

        if      (Start == 1) then
       nStart = 'Warm'
    else if (Start == 2) then
       nStart = 'Hot'
    else
       nStart = 'Cold'
    end if

    call c_f_procpointer(c_qpHx, qpHx)

    pname  = ''
    do j = 1, 8
       if (Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call sqOpt(nStart, qpHx, m, n, neA, nNames,       &
               ncObj, nnH, iObj, ObjAdd, pname,       &
               valA, indA, locA, bl, bu, cObj, Names, &
               eType, hs, x, pi, rc,                  &
               INFO, mincw, miniw, minrw,             &
               nS, nInf, sInf, Obj,                   &
               cw, lencw, iu, leniu, ru, lenru,       &
               cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sqopt

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_snkerq(Start, c_qpHx, c_sqLog,         &
                      m, n, neA, ncObj, nnH,          &
                      iObj, ObjAdd, Prob,             &
                      valA, indA, locA, bl, bu, cObj, &
                      eType, hs, x, pi, rc,           &
                      INFO, nS, nInf, sInf, Obj,      &
                      miniw, minrw,                   &
                      iu, leniu, ru, lenru,           &
                      iw, leniw, rw, lenrw) bind(C,name='f_snkerq')
    integer(c_int), intent(in), value :: Start, m, n, iObj, neA, ncObj, nnH, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    integer(c_int),    intent(in) :: eType(n+m), indA(neA), locA(n+1)
    real(c_double),    intent(in) :: cObj(ncObj), bl(n+m), bu(n+m), valA(neA)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, nS, hs(n)
    real(c_double), intent(inout) :: sInf, x(n), pi(m), rc(n+m)

    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: Obj
    type(c_funptr), value         :: c_qpHx, c_sqLog

    !===========================================================================
    ! Solve the problem with SQOPT.
    !===========================================================================
    integer      :: j, nNames, mincw
    character(8) :: pname, Names(1), Fnames(1)
    character(4) :: nStart

    procedure(iusrHx), pointer :: qpHx
    procedure(sqLog),  pointer :: myLogQ

    nNames = 1

    if      (Start == 1) then
       nStart = 'Warm'
    else if (Start == 2) then
       nStart = 'Hot'
    else
       nStart = 'Cold'
    end if

    myLogQ => null()

    call c_f_procpointer(c_qpHx, qpHx)
    call c_f_procpointer(c_sqLog, myLogQ)

    if (.not. associated(myLogQ)) myLogQ => sqLog

    pname  = ''
    do j = 1, 8
       if (Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call snKerQ(nStart, qpHx, myLogQ,                  &
                m, n, neA, nNames,                     &
                ncObj, nnH, iObj, ObjAdd, pname,       &
                valA, indA, locA, bl, bu, cObj, Names, &
                eType, hs, x, pi, rc,                  &
                INFO, mincw, miniw, minrw,             &
                nS, nInf, sInf, Obj,                   &
                cw, lencw, iu, leniu, ru, lenru,       &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_snkerq

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqend(iw, leniw, rw, lenrw) bind(C,name="f_sqend")
    integer(c_int),    intent(in), value :: leniw, lenrw
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    ! Finish up.
    !===========================================================================
    close(iw(12))  ! print file
    if (iw(13) /= 6) close(iw(13))  ! summary file

  end subroutine f_sqend

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqset(option, len, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_sqset")
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

    call sqSet(buffer, 0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sqset

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqseti(option, len, ivalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_sqseti")
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

    call sqSeti(buffer, ivalue, 0, 0, Errors, &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sqseti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqsetr(option, len, rvalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_sqsetr")
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

    call sqSetr(buffer, rvalue, 0, 0, Errors, &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sqsetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqgetc(option, lin, cvalue, lout, Errors, &
                      iw, leniw, rw, lenrw) bind(C,name="f_sqgetc")
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

    call sqGetC(buffer, buffout, Errors, cw, lencw, iw, leniw, rw, lenrw)

    do j = 1, lout-1
       cvalue(j) = buffout(j:j)
    end do
    cvalue(lout) = c_null_char

  end subroutine f_sqgetc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqgeti(option, len, ivalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_sqgeti")
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

    call sqGetI(buffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sqgeti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_sqgetr(option, len, rvalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_sqgetr")
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

    call sqGetR(buffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_sqgetr

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

end module sqopt_wrapper
