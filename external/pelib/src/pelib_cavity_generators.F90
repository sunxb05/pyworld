module pelib_cavity_generators

    use pelib_precision

    implicit none

    private

    public :: fixtes

    ! Maximum number of tessera
    integer(ip), save :: MXFFTS = 0
    ! Number of centers
    integer(ip), save :: NCENTS = 0
    ! Number of tessera per centers
    integer(ip), public, save :: NTSATM = 60
    ! Number of tessera
    integer(ip), public, save :: NFFTS = 0

contains

!------------------------------------------------------------------------------
!
! The following code has been implemented into the gamess program by Hui Li
! and adapted for use in dalton(pelib)
!
! The tesselation is FIXPVA2
!     Nandun M. Thellamurege and Hui Li THE JOURNAL OF CHEMICAL PHYSICS 137, 246101 (2012)
!
!------------------------------------------------------------------------------
subroutine fixtes(all_centers, all_charges)

    use pelib_constants
    use pelib_options

    real(rp), dimension(:,:), intent(in) :: all_centers
    real(rp), dimension(:), intent(in) :: all_charges

    real(rp), dimension(960) :: AST
    real(rp), dimension(360) :: IDUM
    real(rp), dimension(24) :: THEV
    real(rp), dimension(24) :: FIV
    real(rp), dimension(3,960) :: CDTST
    real(rp), dimension(3,20) :: VT20
    real(rp), dimension(6,60) :: JVT1
    real(rp), dimension(122,3) :: CV
    real(rp), dimension(3,40,960) :: DAIT
    integer(ip), dimension(41,960) :: IDTMP
    real(rp), dimension(:), allocatable :: XTS, YTS, ZTS
    real(rp), dimension(:), allocatable :: AFIX, RFIX
    real(rp), dimension(:,:,:), allocatable :: DAI
    integer(ip), dimension(:), allocatable :: IDATOM
    integer(ip), dimension(:,:), allocatable :: IDDAI
    real(rp), dimension(:,:), allocatable :: TMPTS
    integer(ip) :: II, I, J, INUC, NFFTS, IFFAT
    integer(ip) :: KFFTS, ITS, JJJ
    real(rp) :: TH, FI, CTH, STH, FIR
    real(rp), parameter :: GOLD = 1.618033988749895
    real(rp), parameter :: ONEGOLD = 1.0 / GOLD
    real(rp), parameter :: SQRT13 = 0.577350269189626

    EQUIVALENCE (IDUM(1),JVT1(1,1))

    NCENTS = size(all_charges)
    MXFFTS = NCENTS * NTSATM
    ALLOCATE(XTS(MXFFTS))
    ALLOCATE(YTS(MXFFTS))
    ALLOCATE(ZTS(MXFFTS))
    ALLOCATE(AFIX(MXFFTS))
    ALLOCATE(RFIX(MXFFTS))
    ALLOCATE(IDATOM(MXFFTS))
    ALLOCATE(IDDAI(41,MXFFTS))
    ALLOCATE(DAI(3,40,MXFFTS))
    ALLOCATE(TMPTS(3,MXFFTS))

    THEV  = (/0.6523581398,1.107148718,1.382085796,  &
            &  1.759506858,2.034443936,2.489234514,  &
            &                 0.3261790699,0.5535743589, &
            &  0.8559571251,0.8559571251,1.017221968,&
            &  1.229116717,1.229116717,1.433327788,  &
            &  1.570796327,1.570796327,1.708264866,  &
            &  1.912475937,1.912475937,2.124370686,  &
            &  2.285635528,2.285635528,2.588018295,  &
            &  2.815413584/)

    FIV   = (/ 0.6283185307,0.0000000000,                 &
            &  0.6283185307,0.0000000000,0.6283185307,&
            &  0.0000000000,0.6283185307,0.0000000000,&
            &  0.2520539002,1.004583161,0.6283185307, &
            &  0.3293628477,0.9272742138,0.0000000000,&
            &  0.3141592654,0.9424777961,0.6283185307,&
            &  0.2989556830,0.9576813784,0.0000000000,&
            &  0.3762646305,0.8803724309,0.6283188307,&
            &  0.0000000000/)

    FIR   = 1.256637061

    IDUM = (/1, 6, 2, 32, 36, 37, 1, 2, 3, 33, 32, 38, 1, 3, 4, 34,  &
     &   33, 39, 1, 4, 5, 35, 34, 40, 1, 5, 6, 36, 35, 41, 7, 2, 6, 51,&
     &   42, 37, 8, 3, 2, 47, 43, 38, 9, 4, 3, 48, 44, 39, 10, 5, 4,   &
     &   49, 45, 40, 11, 6, 5, 50, 46, 41, 8, 2, 12, 62, 47, 52, 9,    &
     &   3, 13, 63, 48, 53, 10, 4, 14, 64, 49, 54, 11, 5, 15, 65, 50,  &
     &   55, 7, 6, 16, 66, 51, 56, 7, 12, 2, 42, 57, 52, 8, 13, 3,     &
     &   43, 58, 53, 9, 14, 4, 44, 59, 54, 10, 15, 5, 45, 60, 55, 11,  &
     &   16, 6, 46, 61, 56, 8, 12, 18, 68, 62, 77, 9, 13, 19, 69, 63,  &
     &   78, 10, 14, 20, 70, 64, 79, 11, 15, 21, 71, 65, 80, 7, 16,    &
     &   17, 67, 66, 81, 7, 17, 12, 57, 67, 72, 8, 18, 13, 58, 68, 73, &
     &   9, 19, 14, 59, 69, 74, 10, 20, 15, 60, 70, 75, 11, 21, 16,    &
     &   61, 71, 76, 22, 12, 17, 87, 82, 72, 23, 13, 18, 88, 83, 73,   &
     &   24, 14, 19, 89, 84, 74, 25, 15, 20, 90, 85, 75, 26, 16, 21,   &
     &   91, 86, 76, 22, 18, 12, 82, 92, 77, 23, 19, 13, 83, 93, 78,   &
     &   24, 20, 14, 84, 94, 79, 25, 21, 15, 85, 95, 80, 26, 17, 16,   &
     &   86, 96, 81, 22, 17, 27, 102, 87, 97, 23, 18, 28, 103, 88, 98, &
     &   24, 19, 29, 104, 89, 99, 25, 20, 30, 105, 90, 100, 26, 21,    &
     &   31, 106, 91, 101, 22, 28, 18, 92, 107, 98, 23, 29, 19, 93,    &
     &   108, 99, 24, 30, 20, 94, 109, 100, 25, 31, 21, 95, 110, 101,  &
     &   26, 27, 17, 96, 111, 97, 22, 27, 28, 107, 102, 112, 23, 28,   &
     &   29, 108, 103, 113, 24, 29, 30, 109, 104, 114, 25, 30, 31,     &
     &   110, 105, 115, 26, 31, 27, 111, 106, 116, 122, 28, 27, 117,   &
     &   118, 112, 122, 29, 28, 118, 119, 113, 122, 30, 29, 119, 120,  &
     &   114, 122, 31, 30, 120, 121, 115, 122, 27, 31, 121, 117, 116 /)

!
!     ADAPTED FROM PCM CODE (IN GAMESS)
!     HUI LI, NOV 26, 2011
!
!     COORDINATES OF VERTICES OF TESSERAE IN A SPHERE WITH UNIT RADIUS.
!
!                                    1
!
!                                 4     5
!
!                              3     6     2
!
      CV(  1,1) =  0.0
      CV(  1,2) =  0.0
      CV(  1,3) =  1.0
      CV(122,1) =  0.0
      CV(122,2) =  0.0
      CV(122,3) = -1.0
      II=1
      DO I=1,24
         TH=THEV(I)
         FI=FIV(I)
         CTH=COS(TH)
         STH=SIN(TH)
         DO J=1,5
            FI=FI+FIR
            IF(J.EQ.1) FI=FIV(I)
            II=II+1
            CV(II,1)=STH*COS(FI)
            CV(II,2)=STH*SIN(FI)
            CV(II,3)=CTH
         END DO
      END DO
!
      IF(NTSATM.EQ.4) THEN
         VT20(1,1)  = 1.0
         VT20(2,1)  = 1.0
         VT20(3,1)  = 1.0
         VT20(1,2)  =-1.0
         VT20(2,2)  =-1.0
         VT20(3,2)  = 1.0
         VT20(1,3)  =-1.0
         VT20(2,3)  = 1.0
         VT20(3,3)  =-1.0
         VT20(1,4)  = 1.0
         VT20(2,4)  =-1.0
         VT20(3,4)  =-1.0
         CALL DSCAL(12,SQRT13,VT20,1)
      END IF
!
      IF(NTSATM.EQ.6) THEN
         VT20(1,1)  = 1.0
         VT20(2,1)  = 0.0
         VT20(3,1)  = 0.0
         VT20(1,2)  =-1.0
         VT20(2,2)  = 0.0
         VT20(3,2)  = 0.0
         VT20(1,3)  = 0.0
         VT20(2,3)  = 1.0
         VT20(3,3)  = 0.0
         VT20(1,4)  = 0.0
         VT20(2,4)  =-1.0
         VT20(3,4)  = 0.0
         VT20(1,5)  = 0.0
         VT20(2,5)  = 0.0
         VT20(3,5)  = 1.0
         VT20(1,6)  = 0.0
         VT20(2,6)  = 0.0
         VT20(3,6)  =-1.0
      END IF
!
      IF(NTSATM.EQ.8) THEN
         VT20(1,1)  = 1.0
         VT20(2,1)  = 1.0
         VT20(3,1)  = 1.0
         VT20(1,2)  =-1.0
         VT20(2,2)  = 1.0
         VT20(3,2)  = 1.0
         VT20(1,3)  = 1.0
         VT20(2,3)  =-1.0
         VT20(3,3)  = 1.0
         VT20(1,4)  = 1.0
         VT20(2,4)  = 1.0
         VT20(3,4)  =-1.0
         VT20(1,5)  =-1.0
         VT20(2,5)  =-1.0
         VT20(3,5)  = 1.0
         VT20(1,6)  =-1.0
         VT20(2,6)  = 1.0
         VT20(3,6)  =-1.0
         VT20(1,7)  = 1.0
         VT20(2,7)  =-1.0
         VT20(3,7)  =-1.0
         VT20(1,8)  =-1.0
         VT20(2,8)  =-1.0
         VT20(3,8)  =-1.0
         CALL DSCAL(24,SQRT13,VT20,1)
      END IF
!
      IF(NTSATM.EQ.12) THEN
         VT20(1,1)  = 0.0
         VT20(2,1)  = 1.0
         VT20(3,1)  =    GOLD
         VT20(1,2)  = 0.0
         VT20(2,2)  = 1.0
         VT20(3,2)  =   -GOLD
         VT20(1,3)  = 0.0
         VT20(2,3)  =-1.0
         VT20(3,3)  =    GOLD
         VT20(1,4)  = 0.0
         VT20(2,4)  =-1.0
         VT20(3,4)  =   -GOLD
         VT20(3,5)  = 0.0
         VT20(1,5)  = 1.0
         VT20(2,5)  =    GOLD
         VT20(3,6)  = 0.0
         VT20(1,6)  = 1.0
         VT20(2,6)  =   -GOLD
         VT20(3,7)  = 0.0
         VT20(1,7)  =-1.0
         VT20(2,7)  =    GOLD
         VT20(3,8)  = 0.0
         VT20(1,8)  =-1.0
         VT20(2,8)  =   -GOLD
         VT20(2,9)  = 0.0
         VT20(3,9)  = 1.0
         VT20(1,9)  =    GOLD
         VT20(2,10) = 0.0
         VT20(3,10) = 1.0
         VT20(1,10) =   -GOLD
         VT20(2,11) = 0.0
         VT20(3,11) =-1.0
         VT20(1,11) =    GOLD
         VT20(2,12) = 0.0
         VT20(3,12) =-1.0
         VT20(1,12) =   -GOLD
         CALL DSCAL(36,0.525731112119134,VT20,1)
      END IF
!
      IF(NTSATM.EQ.20) THEN
         VT20(1,1)  = 1.0
         VT20(2,1)  = 1.0
         VT20(3,1)  = 1.0
         VT20(1,2)  = 1.0
         VT20(2,2)  = 1.0
         VT20(3,2)  =-1.0
         VT20(1,3)  = 1.0
         VT20(2,3)  =-1.0
         VT20(3,3)  = 1.0
         VT20(1,4)  =-1.0
         VT20(2,4)  = 1.0
         VT20(3,4)  = 1.0
         VT20(1,5)  = 1.0
         VT20(2,5)  =-1.0
         VT20(3,5)  =-1.0
         VT20(1,6)  =-1.0
         VT20(2,6)  = 1.0
         VT20(3,6)  =-1.0
         VT20(1,7)  =-1.0
         VT20(2,7)  =-1.0
         VT20(3,7)  = 1.0
         VT20(1,8)  =-1.0
         VT20(2,8)  =-1.0
         VT20(3,8)  =-1.0
         VT20(1,9)  = 0.0
         VT20(2,9)  = ONEGOLD
         VT20(3,9)  =    GOLD
         VT20(1,10) = 0.0
         VT20(2,10) = ONEGOLD
         VT20(3,10) =   -GOLD
         VT20(1,11) = 0.0
         VT20(2,11) =-ONEGOLD
         VT20(3,11) =    GOLD
         VT20(1,12) = 0.0
         VT20(2,12) =-ONEGOLD
         VT20(3,12) =   -GOLD
         VT20(3,13) = 0.0
         VT20(1,13) = ONEGOLD
         VT20(2,13) =    GOLD
         VT20(3,14) = 0.0
         VT20(1,14) = ONEGOLD
         VT20(2,14) =   -GOLD
         VT20(3,15) = 0.0
         VT20(1,15) =-ONEGOLD
         VT20(2,15) =    GOLD
         VT20(3,16) = 0.0
         VT20(1,16) =-ONEGOLD
         VT20(2,16) =   -GOLD
         VT20(2,17) = 0.0
         VT20(3,17) = ONEGOLD
         VT20(1,17) =    GOLD
         VT20(2,18) = 0.0
         VT20(3,18) = ONEGOLD
         VT20(1,18) =   -GOLD
         VT20(2,19) = 0.0
         VT20(3,19) =-ONEGOLD
         VT20(1,19) =    GOLD
         VT20(2,20) = 0.0
         VT20(3,20) =-ONEGOLD
         VT20(1,20) =   -GOLD
         CALL DSCAL(60,SQRT13,VT20,1)
      END IF
!
      DO I = 1, NCENTS
         INUC = INT(all_charges(I))
         RFIX(I) = 2.400*aa2bohr

         IF(INUC.EQ. 0) RFIX(I) = 0.001*aa2bohr
         IF(INUC.EQ. 1) RFIX(I) = 1.400*aa2bohr
         IF(INUC.EQ. 3) RFIX(I) = 1.400*aa2bohr
         IF(INUC.EQ. 4) RFIX(I) = 1.400*aa2bohr
         IF(INUC.EQ. 5) RFIX(I) = 1.400*aa2bohr
         IF(INUC.EQ. 6) RFIX(I) = 2.100*aa2bohr
         IF(INUC.EQ. 7) RFIX(I) = 2.000*aa2bohr
         IF(INUC.EQ. 8) RFIX(I) = 1.900*aa2bohr
         IF(INUC.EQ. 9) RFIX(I) = 1.800*aa2bohr
         IF(INUC.EQ.10) RFIX(I) = 1.800*aa2bohr
         IF(INUC.EQ.11) RFIX(I) = 1.800*aa2bohr
         IF(INUC.EQ.12) RFIX(I) = 1.800*aa2bohr
         IF(INUC.EQ.13) RFIX(I) = 1.800*aa2bohr
         IF(INUC.EQ.14) RFIX(I) = 2.000*aa2bohr
         IF(INUC.EQ.15) RFIX(I) = 2.200*aa2bohr
         IF(INUC.EQ.16) RFIX(I) = 2.400*aa2bohr
         IF(INUC.EQ.17) RFIX(I) = 2.760*aa2bohr
         IF(INUC.EQ.18) RFIX(I) = 3.000*aa2bohr
      END DO
!
!     -- NOTE: 'ME' STARTS AT 0 --
!
      xts = 0.0
      yts = 0.0
      zts = 0.0
      afix = 0.0
      idatom = 0
      dai = 0.0
      iddai = 0
!      NFFTS   = ME*(MXFFTS/NPROC - 1)
!      IPCOUNT = ME - 1
      DO IFFAT = 1, NCENTS
         IF(RFIX(IFFAT).LE.0.1) cycle
!         IF(GOPARR) THEN
!           IPCOUNT = IPCOUNT + 1
!           IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 500
!         END IF
         CALL FIXPVA2(IFFAT,all_centers,JVT1,CV,CDTST,AST,  &
                    & XTS,YTS,ZTS,AFIX,RFIX,IDATOM,  &
                    & DAI,IDDAI,DAIT,IDTMP,TMPTS,VT20)
      end do
!     Parallel code not available in dalton yet
!     IF(GOPARR) THEN
!        CALL DDI_GSUMF(2418,XTS   ,      MXFFTS)
!        CALL DDI_GSUMF(2419,YTS   ,      MXFFTS)
!        CALL DDI_GSUMF(2420,ZTS   ,      MXFFTS)
!        CALL DDI_GSUMF(2421,AFIX  ,      MXFFTS)
!        CALL DDI_GSUMI(2422,IDATOM,      MXFFTS)
!        CALL DDI_GSUMF(2423,DAI   , 3*40*MXFFTS)
!        CALL DDI_GSUMI(2424,IDDAI ,   41*MXFFTS)
!     END IF
!

      KFFTS = 0
      DO ITS=1,MXFFTS
         IF(AFIX(ITS).LT.1.0e-4) cycle  ! 1.0e-04
         KFFTS = KFFTS + 1
         IF(KFFTS.GT.MXFFTS) THEN
            ERROR STOP 'FIXTES: PLEASE INCREASE MXFFTS.'
         END IF
         XTS(KFFTS)      = XTS(ITS)
         YTS(KFFTS)      = YTS(ITS)
         ZTS(KFFTS)      = ZTS(ITS)
         AFIX(KFFTS)     = AFIX(ITS)
         IDATOM(KFFTS)   = IDATOM(ITS)
         IDDAI(41,KFFTS) = IDDAI(41,ITS)
         DO JJJ = 1, IDDAI(41,ITS)
            IDDAI(JJJ,KFFTS) = IDDAI(JJJ,ITS)
            DAI(1,JJJ,KFFTS) = DAI(1,JJJ,ITS)
            DAI(2,JJJ,KFFTS) = DAI(2,JJJ,ITS)
            DAI(3,JJJ,KFFTS) = DAI(3,JJJ,ITS)
         END DO
      end do
      NFFTS = KFFTS   ! NFFTS IS GLOBAL

!      FIXA = 0.0
!      DO IFFTS=1,NFFTS
!         FIXA = FIXA + AFIX(IFFTS)
!      END DO
!      write(luout,*) 'Surface point coordinates'
!      do iffts= 1, nffts
!         write(luout,*) XTS(IFFTS), YTS(IFFTS), ZTS(IFFTS)
!      end do
!      FIXA = FIXA*bohr2aa**2
!
!      write(luout,*) 'Surface area FIXA', FIXA, '/A**2'
!      write(luout,*) 'Number of surface charges', NFFTS
!      call openfile(surfile, surf, 'new', 'formatted')
!      rewind(lu)
!      write(surf,'(i6)') nffts
!      write(surf,'(a)') 'AA'
!      do i = 1, nffts
!          write(surf,'(4f12.6)') XTS(I), YTS(I), ZTS(I), AFIX(I)
!      end do
      nsurp = nffts
      allocate(Rsp(3,nsurp))
      allocate(Sa(nsurp))
      allocate(idatm(nsurp))
      allocate(Radsp(nsurp))
      do i = 1, nsurp
         Sa(i) = AFIX(i)
         Rsp(1,i) = XTS(i)
         Rsp(2,i) = YTS(i)
         Rsp(3,i) = ZTS(i)
         idatm(i) = idatom(i)
         Radsp(i) = rfix(idatom(i))
      end do
!      close(lu)

end subroutine fixtes

!-----------------------------------------------------------------------------

subroutine fixpva2(iffat,cord,jvt1,cv,cdtst,ast,&
                & xts,yts,zts,afix,rfix,idatom,&
                & dai, iddai, dait,idtmp,tmpts,vt20)

    use pelib_constants, only: pi

    real(rp), dimension(:), intent(out) :: XTS, YTS, ZTS, AFIX, RFIX

    real(rp), dimension(3,NCENTS) :: CORD
    real(rp), dimension(960) :: AST
    real(rp), dimension(3,960) :: CDTST
    real(rp), dimension(3,40,960) :: DAIT
    integer(ip), dimension(41,960) :: IDTMP
    integer(ip), dimension(41) :: IDTMPTS
    real(rp), dimension(6,60) :: JVT1
    real(rp), dimension(122,3) :: CV
    real(rp), dimension(3,20) :: VT20
    real(rp), dimension(3,MXFFTS) :: TMPTS
    real(rp) :: XII, YII, ZII, RII, AREA0, DUMMY
    real(rp) :: P1X, P1Y, P1Z, P2X, P2Y, P2Z, P3X, P3Y, P3Z, P4X, P4Y, P4Z
    real(rp) :: P5X, P5Y, P5Z, P6X, P6Y, P6Z, P7X, P7Y, P7Z
    real(rp) :: PTS11, PTS21, PTS31, PTS12, PTS22, PTS32, PTS13, PTS23, PTS33
    real(rp) :: DNORM4, DNORM5, DNORM6, DNORM7, FACTOR, SCALE4, SCALE5, SCALE6, SCALE7
    integer(ip) :: IFFAT, KFFAT, ITS, KKK, MFFAT, III
    integer(ip) :: N1, N2, N3
    integer(ip) :: JJJ, JTS, KTS, LTS
    real(rp) :: FOURPI = 4.0 * pi
    real(rp) :: ONE3RD=1.0/3.0
    real(rp), dimension(3,40,MXFFTS) :: DAI
    integer(ip), dimension(MXFFTS) :: IDATOM
    integer(ip), dimension(41,MXFFTS) :: IDDAI

!     PARTITION OF THE CAVITY SURFACE INTO TESSERAE
!     HUI LI, NOV 26, 2011

     XII    = CORD(1,IFFAT)
     YII    = CORD(2,IFFAT)
     ZII    = CORD(3,IFFAT)
     RII    = RFIX(IFFAT)
     AREA0  = FOURPI*RII*RII/NTSATM
     cdtst  = 0.0
     ast    = 0.0
     dait   = 0.0
     idtmp  = 0



!     --- USE 20 TESSERAE FOR EACH SPHERE ---
     IF(NTSATM.EQ.20) THEN
     DO ITS = 1, NTSATM

        CDTST(1,ITS) = XII + VT20(1,ITS)*RII
        CDTST(2,ITS) = YII + VT20(2,ITS)*RII
        CDTST(3,ITS) = ZII + VT20(3,ITS)*RII

        tmpts = 0.0
        idtmpts = 0
        KFFAT=IFFAT
        AST(ITS) = AREA0
        CALL FIXPVASWF(KFFAT,CORD,CDTST,ITS,AST(ITS),RFIX,&
                     & TMPTS,IDTMPTS)
        IF(AST(ITS).GT.0.9e-4) THEN
           DO KKK = 1, IDTMPTS(41)
              MFFAT = IDTMPTS(KKK)
              DUMMY=ABS(TMPTS(1,MFFAT))+&
                  & ABS(TMPTS(2,MFFAT))+&
                  & ABS(TMPTS(3,MFFAT))
              IF(DUMMY.GT.1.0e-5) THEN
                 IDTMP(41,ITS) = IDTMP(41,ITS) + 1
                 III = IDTMP(41,ITS)
                 IDTMP(III,ITS)  = MFFAT
                 DAIT(1,III,ITS) = TMPTS(1,MFFAT)
                 DAIT(2,III,ITS) = TMPTS(2,MFFAT)
                 DAIT(3,III,ITS) = TMPTS(3,MFFAT)
              END IF
           ENDDO
        END IF
     END DO
     END IF



!     --- USE 60 TESSERAE FOR EACH SPHERE ---
     IF(NTSATM.EQ.60) THEN
     DO ITS = 1, NTSATM

        N1 = JVT1(1,ITS)
        N2 = JVT1(2,ITS)
        N3 = JVT1(3,ITS)
        P4X = (CV(N1,1)+CV(N2,1)+CV(N3,1))*ONE3RD
        P4Y = (CV(N1,2)+CV(N2,2)+CV(N3,2))*ONE3RD
        P4Z = (CV(N1,3)+CV(N2,3)+CV(N3,3))*ONE3RD
        DNORM4 = P4X**2+P4Y**2+P4Z**2
        SCALE4 = 1.0/SQRT(DNORM4)
        CDTST(1,ITS) = XII + P4X*SCALE4*RII
        CDTST(2,ITS) = YII + P4Y*SCALE4*RII
        CDTST(3,ITS) = ZII + P4Z*SCALE4*RII

        tmpts = 0.0
        idtmpts = 0
        KFFAT=IFFAT
        AST(ITS) = AREA0
        CALL FIXPVASWF(KFFAT,CORD,CDTST,ITS,AST(ITS),RFIX,&
                     & TMPTS,IDTMPTS)
        IF(AST(ITS).GT.0.9D-04) THEN
           DO KKK = 1, IDTMPTS(41)
              MFFAT = IDTMPTS(KKK)
              DUMMY=ABS(TMPTS(1,MFFAT))+&
                  & ABS(TMPTS(2,MFFAT))+&
                  & ABS(TMPTS(3,MFFAT))
              IF(DUMMY.GT.1.0e-05) THEN
                 IDTMP(41,ITS) = IDTMP(41,ITS) + 1
                 III = IDTMP(41,ITS)
                 IDTMP(III,ITS)  = MFFAT
                 DAIT(1,III,ITS) = TMPTS(1,MFFAT)
                 DAIT(2,III,ITS) = TMPTS(2,MFFAT)
                 DAIT(3,III,ITS) = TMPTS(3,MFFAT)
              END IF
           ENDDO
        END IF
     end do
     END IF



!     --- USE 240 TESSERAE FOR EACH SPHERE ---
     IF(NTSATM.EQ.240) THEN
     DO KTS = 1, 60

     DO JTS = 1, 4
        ITS = (KTS-1)*4 + JTS
        IF(JTS.EQ.1) THEN
           N1 = JVT1(1,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(5,KTS)
        END IF
        IF(JTS.EQ.2) THEN
           N1 = JVT1(6,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(5,KTS)
        END IF
        IF(JTS.EQ.3) THEN
           N1 = JVT1(3,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(6,KTS)
        END IF
        IF(JTS.EQ.4) THEN
           N1 = JVT1(2,KTS)
           N2 = JVT1(6,KTS)
           N3 = JVT1(5,KTS)
        END IF
        P4X = (CV(N1,1)+CV(N2,1)+CV(N3,1))*ONE3RD
        P4Y = (CV(N1,2)+CV(N2,2)+CV(N3,2))*ONE3RD
        P4Z = (CV(N1,3)+CV(N2,3)+CV(N3,3))*ONE3RD
        DNORM4 = P4X**2+P4Y**2+P4Z**2
        SCALE4 = 1.0/SQRT(DNORM4)
        CDTST(1,ITS) = XII + P4X*SCALE4*RII
        CDTST(2,ITS) = YII + P4Y*SCALE4*RII
        CDTST(3,ITS) = ZII + P4Z*SCALE4*RII
        KKK = MOD(ITS,4)
        IF(KKK.EQ.1) FACTOR = 0.97527227808
        IF(KKK.EQ.2) FACTOR = 1.04680294088
        IF(KKK.EQ.3) FACTOR = 0.98896281888
        IF(KKK.EQ.0) FACTOR = 0.98896281888

        tmpts = 0.0
        idtmpts = 0
        KFFAT=IFFAT
        AST(ITS) = AREA0*FACTOR
        CALL FIXPVASWF(KFFAT,CORD,CDTST,ITS,AST(ITS),RFIX,&
                     & TMPTS,IDTMPTS)
        IF(AST(ITS).GT.0.9D-04) THEN
           DO KKK = 1, IDTMPTS(41)
              MFFAT = IDTMPTS(KKK)
              DUMMY=ABS(TMPTS(1,MFFAT))+&
                  & ABS(TMPTS(2,MFFAT))+&
                  & ABS(TMPTS(3,MFFAT))
              IF(DUMMY.GT.1.0e-05) THEN
                 IDTMP(41,ITS) = IDTMP(41,ITS) + 1
                 III = IDTMP(41,ITS)
                 IDTMP(III,ITS)  = MFFAT
                 DAIT(1,III,ITS) = TMPTS(1,MFFAT)
                 DAIT(2,III,ITS) = TMPTS(2,MFFAT)
                 DAIT(3,III,ITS) = TMPTS(3,MFFAT)
              END IF
           ENDDO
        END IF
     END DO
     END DO
     END IF



!     --- USE 960 TESSERAE FOR EACH SPHERE ---
     IF(NTSATM.EQ.960) THEN
     DO KTS = 1, 60

     DO JTS = 1, 4
        IF(JTS.EQ.1) THEN
           N1 = JVT1(1,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(5,KTS)
        END IF
        IF(JTS.EQ.2) THEN
           N1 = JVT1(6,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(5,KTS)
        END IF
        IF(JTS.EQ.3) THEN
           N1 = JVT1(3,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(6,KTS)
        END IF
        IF(JTS.EQ.4) THEN
           N1 = JVT1(2,KTS)
           N2 = JVT1(6,KTS)
           N3 = JVT1(5,KTS)
        END IF
        P1X = CV(N1,1)
        P1Y = CV(N1,2)
        P1Z = CV(N1,3)
        P2X = CV(N2,1)
        P2Y = CV(N2,2)
        P2Z = CV(N2,3)
        P3X = CV(N3,1)
        P3Y = CV(N3,2)
        P3Z = CV(N3,3)
!          COMPUTE THE COORDINATES OF POINTS 4, 5, 6
!
!                         1
!
!                      4     5
!
!                   3     6     2


        P4X = (P1X+P3X)/2.0
        P4Y = (P1Y+P3Y)/2.0
        P4Z = (P1Z+P3Z)/2.0
        P5X = (P1X+P2X)/2.0
        P5Y = (P1Y+P2Y)/2.0
        P5Z = (P1Z+P2Z)/2.0
        P6X = (P2X+P3X)/2.0
        P6Y = (P2Y+P3Y)/2.0
        P6Z = (P2Z+P3Z)/2.0

        DNORM4 = P4X**2+P4Y**2+P4Z**2
        DNORM5 = P5X**2+P5Y**2+P5Z**2
        DNORM6 = P6X**2+P6Y**2+P6Z**2
        SCALE4 = 1.0/SQRT(DNORM4)
        SCALE5 = 1.0/SQRT(DNORM5)
        SCALE6 = 1.0/SQRT(DNORM6)
        P4X = P4X*SCALE4
        P4Y = P4Y*SCALE4
        P4Z = P4Z*SCALE4
        P5X = P5X*SCALE5
        P5Y = P5Y*SCALE5
        P5Z = P5Z*SCALE5
        P6X = P6X*SCALE6
        P6Y = P6Y*SCALE6
        P6Z = P6Z*SCALE6

    DO  LTS = 1, 4
        ITS = ((KTS-1)*4 + JTS-1)*4 + LTS
        IF(LTS.EQ.1) THEN
          PTS11=P1X
          PTS21=P1Y
          PTS31=P1Z
          PTS12=P4X
          PTS22=P4Y
          PTS32=P4Z
          PTS13=P5X
          PTS23=P5Y
          PTS33=P5Z
        ELSE IF(LTS.EQ.2) THEN
          PTS11=P6X
          PTS21=P6Y
          PTS31=P6Z
          PTS12=P4X
          PTS22=P4Y
          PTS32=P4Z
          PTS13=P5X
          PTS23=P5Y
          PTS33=P5Z
        ELSE IF(LTS.EQ.3) THEN
          PTS11=P3X
          PTS21=P3Y
          PTS31=P3Z
          PTS12=P4X
          PTS22=P4Y
          PTS32=P4Z
          PTS13=P6X
          PTS23=P6Y
          PTS33=P6Z
        ELSE IF(LTS.EQ.4) THEN
          PTS11=P2X
          PTS21=P2Y
          PTS31=P2Z
          PTS12=P6X
          PTS22=P6Y
          PTS32=P6Z
          PTS13=P5X
          PTS23=P5Y
          PTS33=P5Z
        END IF

        P7X = (PTS11+PTS12+PTS13)*ONE3RD
        P7Y = (PTS21+PTS22+PTS23)*ONE3RD
        P7Z = (PTS31+PTS32+PTS33)*ONE3RD
        DNORM7 = P7X**2+P7Y**2+P7Z**2
        SCALE7 = 1.0/SQRT(DNORM7)
        CDTST(1,ITS) = XII + P7X*SCALE7*RII
        CDTST(2,ITS) = YII + P7Y*SCALE7*RII
        CDTST(3,ITS) = ZII + P7Z*SCALE7*RII
        KKK = MOD(ITS,16)
        IF(KKK.EQ. 1) FACTOR = 0.96853843384
        IF(KKK.EQ. 2) FACTOR = 0.98634758299
        IF(KKK.EQ. 3) FACTOR = 0.97310155328
        IF(KKK.EQ. 4) FACTOR = 0.97310155328
        IF(KKK.EQ. 5) FACTOR = 1.04026074020
        IF(KKK.EQ. 6) FACTOR = 1.05941468448
        IF(KKK.EQ. 7) FACTOR = 1.04376817374
        IF(KKK.EQ. 8) FACTOR = 1.04376817369
        IF(KKK.EQ. 9) FACTOR = 0.98544774054
        IF(KKK.EQ.10) FACTOR = 1.00020007354
        IF(KKK.EQ.11) FACTOR = 0.98676349755
        IF(KKK.EQ.12) FACTOR = 0.98343824680
        IF(KKK.EQ.13) FACTOR = 0.98544774307
        IF(KKK.EQ.14) FACTOR = 1.00020007615
        IF(KKK.EQ.15) FACTOR = 0.98343824928
        IF(KKK.EQ. 0) FACTOR = 0.98676350013

        tmpts = 0.0
        idtmpts = 0
        KFFAT=IFFAT
        AST(ITS) = AREA0*FACTOR
        CALL FIXPVASWF(KFFAT,CORD,CDTST,ITS,AST(ITS),RFIX,&
                     & TMPTS,IDTMPTS)
       IF(AST(ITS).GT.0.9D-04) THEN
           DO KKK = 1, IDTMPTS(41)
              MFFAT = IDTMPTS(KKK)
              DUMMY=ABS(TMPTS(1,MFFAT))+&
                  & ABS(TMPTS(2,MFFAT))+&
                  & ABS(TMPTS(3,MFFAT))
              IF(DUMMY.GT.1.0e-05) THEN
                 IDTMP(41,ITS) = IDTMP(41,ITS) + 1
                 III = IDTMP(41,ITS)
                 IDTMP(III,ITS)  = MFFAT
                 DAIT(1,III,ITS) = TMPTS(1,MFFAT)
                 DAIT(2,III,ITS) = TMPTS(2,MFFAT)
                 DAIT(3,III,ITS) = TMPTS(3,MFFAT)
              END IF
           ENDDO
        END IF
     end do  ! lts = 1, 4
     end do  ! jts = 1, 4
     end do  ! kts = 1, 60
     END IF


     DO ITS=1,NTSATM
        IF(AST(ITS).LT.0.95D-04) cycle      !  0.95D-04
        NFFTS = NFFTS + 1
        if (NFFTS .GT. MXFFTS) then
           ERROR STOP 'FIXPVA2: PLEASE INCREASE MXFFTS.'
        end if
        XTS(NFFTS)      = CDTST(1,ITS)
        YTS(NFFTS)      = CDTST(2,ITS)
        ZTS(NFFTS)      = CDTST(3,ITS)
        AFIX(NFFTS)     = AST(ITS)
        IDATOM(NFFTS)   = IFFAT
        IDDAI(41,NFFTS) = IDTMP(41,ITS)
        DO JJJ = 1, IDTMP(41,ITS)
           IDDAI(JJJ,NFFTS) = IDTMP(JJJ,ITS)
           DAI(1,JJJ,NFFTS) = DAIT(1,JJJ,ITS)
           DAI(2,JJJ,NFFTS) = DAIT(2,JJJ,ITS)
           DAI(3,JJJ,NFFTS) = DAIT(3,JJJ,ITS)
        ENDDO
     END DO

     RETURN

end subroutine fixpva2

!-----------------------------------------------------------------------------

subroutine FIXPVASWF(IFFAT,CORD,CDTST,ITS,AREA,RFIX,TMP,IDTMPTS)

    use pelib_constants

   real(rp), intent(out) :: AREA
   integer(ip), intent(in) :: iffat
   real(rp), dimension(3,NCENTS) :: cord, tmp
   real(rp), dimension(3,960) :: CDTST
   real(rp), dimension(MXFFTS) :: RFIX
   integer(ip), dimension(41) :: IDTMPTS
   real(rp) :: ZERO=0.0
   real(rp) :: ONE=1.0
   real(rp) :: TWO=2.0
   real(rp) :: PT5=0.5
   real(rp) :: TEN=10.0
   real(rp) :: DISM0, ONEDISM0, DISN1, DISN2, RA, RB
   real(rp) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X5, Y5, Z5
   real(rp) :: DISC2, ONEDISC2, DISC, ONEDISC, DISC3, DIS132, DIS13, DUM
   real(rp) :: DISN, SWF2, DSWF2X, DSWF2Y, DSWF2Z
   real(rp) :: DUW, DX5DX3, DY5DX3, DZ5DX3, DX5DY3, DY5DY3, DZ5DY3
   real(rp) :: DX5DZ3, DY5DZ3, DZ5DZ3, DDISN2X, DDISN2Y, DDISN2Z
   real(rp) :: SWF1, DSWF1X, DSWF1Y, DSWF1Z, DISD
   real(rp) :: FB, FA, R1, R2, R3, R4, R5, R6, R7, R8, FAB, ONEFAB, DISM
   real(rp) :: ONEDISM, ONEDISD, DDDIS
   real(rp) :: CX3, CY3, CZ3, DX3, DY3, DZ3, R1X, R1Y, R1Z
   real(rp) :: R2X, R2Y, R2Z, R3X, R3Y, R3Z, R4X, R4Y, R4Z, R5X, R5Y, R5Z
   real(rp) :: R6X, R6Y, R6Z, R7X, R7Y, R7Z, R8X, R8Y, R8Z
   real(rp) :: FABX, FABY, FABZ, TEMP, DDISM2X, DDISM2Y, DDISM2Z, DUMMY
   integer(ip) :: JFFAT, KKK, ITS, KFFAT

!
!     DETERMINE THE AREA AND AREA DERIVATIVES FOR A TESSERA
!     USING FIXPVA2
!     HUI LI, NOV 27, 2011, LINCOLN
!
!     P1 = X1, Y1, Z1 (THE CURRENT TESSERA)
!     P2 = X2, Y2, Z2 (THE CENTER OF THE SPHERE OF THE TESSERA)
!     P3 = X3, Y3, Z3 (THE CENTER OF A NEIGHBOR SPHERE)
!
!
      DISM0 = 0.30*aa2bohr
      ONEDISM0 = ONE/DISM0
      DISN1 = 0.70*aa2bohr
      DISN2 = 1.20*aa2bohr    ! DISN2 MUST BE < RA
!
      X1 = CDTST(1,ITS)
      Y1 = CDTST(2,ITS)
      Z1 = CDTST(3,ITS)
      X2 = CORD(1,IFFAT)
      Y2 = CORD(2,IFFAT)
      Z2 = CORD(3,IFFAT)
      RA = RFIX(IFFAT)
!
      DO JFFAT=1,NCENTS
!
!        IF AREA=0, RETURN
!        SKIP DISTANT SPHERES
!
         IF(JFFAT.EQ.IFFAT .OR. RFIX(JFFAT).LT.0.10) cycle
         X3 = CORD(1,JFFAT)
         Y3 = CORD(2,JFFAT)
         Z3 = CORD(3,JFFAT)
         RB = RFIX(JFFAT)
         IF(ABS(X2-X3).GT.(RA+RB+DISN2)) cycle
         IF(ABS(Y2-Y3).GT.(RA+RB+DISN2)) cycle
         IF(ABS(Z2-Z3).GT.(RA+RB+DISN2)) cycle
         DISC2 = (X2-X3)**2 + (Y2-Y3)**2 + (Z2-Z3)**2
         ONEDISC2 = ONE/DISC2
         IF(DISC2.GE.(RA+RB+DISN2)**2) cycle
         IF(DISC2.LE.(RA-RB)**2 .AND. RB.LT.RA) cycle
         IF(DISC2.LE.(RA-RB)**2 .AND. RB.GT.RA) THEN
            AREA = ZERO
            RETURN
         END IF
         DISC   = SQRT(DISC2)
         ONEDISC= ONE/DISC
         DISC3  = DISC*DISC2
         DIS132 = (X1-X3)**2+(Y1-Y3)**2+(Z1-Z3)**2
         DIS13  = SQRT(DIS132)
!
         DUM = RB*ONEDISC
         X5  = X3 + (X2-X3)*DUM
         Y5  = Y3 + (Y2-Y3)*DUM
         Z5  = Z3 + (Z2-Z3)*DUM
         DISN = SQRT((X1-X5)**2+(Y1-Y5)**2+(Z1-Z5)**2)
         IF(DISN.GE.DISN2 .OR. DISC.LE.RB) THEN
            SWF2   = ONE
            DSWF2X = ZERO
            DSWF2Y = ZERO
            DSWF2Z = ZERO
         ELSE IF (DISN.LE.DISN1) THEN
            AREA   = ZERO
            RETURN
         ELSE
            DUW  = (DISN**2 - DISN1**2)/(DISN2**2-DISN1**2)
            SWF2 = 10.0*DUW**3 - 15.0*DUW**4 + 6.0*DUW**5
            DUM    = RB/DISC3
            DX5DX3 = ONE - RB*ONEDISC + DUM*(X2-X3)*(X2-X3)
            DY5DX3 =                    DUM*(Y2-Y3)*(X2-X3)
            DZ5DX3 =                    DUM*(Z2-Z3)*(X2-X3)
            DX5DY3 =                    DUM*(X2-X3)*(Y2-Y3)
            DY5DY3 = ONE - RB*ONEDISC + DUM*(Y2-Y3)*(Y2-Y3)
            DZ5DY3 =                    DUM*(Z2-Z3)*(Y2-Y3)
            DX5DZ3 =                    DUM*(X2-X3)*(Z2-Z3)
            DY5DZ3 =                    DUM*(Y2-Y3)*(Z2-Z3)
            DZ5DZ3 = ONE - RB*ONEDISC + DUM*(Z2-Z3)*(Z2-Z3)
            DDISN2X=-TWO*((X1-X5)*DX5DX3+(Y1-Y5)*DY5DX3+(Z1-Z5)*DZ5DX3)
            DDISN2Y=-TWO*((X1-X5)*DX5DY3+(Y1-Y5)*DY5DY3+(Z1-Z5)*DZ5DY3)
            DDISN2Z=-TWO*((X1-X5)*DX5DZ3+(Y1-Y5)*DY5DZ3+(Z1-Z5)*DZ5DZ3)
            DUM = (30.0*DUW**2-60.0*DUW**3+30.0*DUW**4)/&
               &  (DISN2**2-DISN1**2)
            DSWF2X = DUM*DDISN2X
            DSWF2Y = DUM*DDISN2Y
            DSWF2Z = DUM*DDISN2Z
         END IF
!
!
         IF(DISC.GE.(RA+RB)) THEN
            SWF1   = ONE
            DSWF1X = ZERO
            DSWF1Y = ZERO
            DSWF1Z = ZERO
         ELSE
            DISD = DIS13
            FB =  RA**2 + DISC2 - DISD**2
            FA =  RA**2 + DISC2 - RB**2
            R1 =  RA    + DISC  + DISD
            R2 =  RA    + DISC  - DISD
            R3 =  RA    - DISC  + DISD
            R4 = -RA    + DISC  + DISD
            R5 =  RA    + DISC  + RB
            R6 =  RA    + DISC  - RB
            R7 =  RA    - DISC  + RB
            R8 = -RA    + DISC  + RB
            FAB = SQRT(ABS(R1*R2*R3*R4*R5*R6*R7*R8))
            ONEFAB = ONE/FAB
            DISM = SQRT(ABS(TWO*RA*RA - (FA*FB + FAB)*PT5*ONEDISC2))
            IF(DISM.GE.DISM0) THEN
               IF(DIS13.LT.RB) THEN
                  AREA   = ZERO
                  RETURN
               END IF
               SWF1   = ONE
               DSWF1X = ZERO
               DSWF1Y = ZERO
               DSWF1Z = ZERO
            ELSE
               ONEDISM  = 1.0e+04
               IF(DISM.GT.1.0e-04) ONEDISM  = ONE/DISM
               ONEDISD  = ONE/DISD
               DDDIS    = (DISM-DISM0)*ONEDISM0
               SWF1     = PT5*DDDIS*DDDIS
               IF(DIS13.GE.RB) SWF1 = ONE - SWF1
               CX3 =  (X3-X2)*ONEDISC
               CY3 =  (Y3-Y2)*ONEDISC
               CZ3 =  (Z3-Z2)*ONEDISC
               DX3 =  (X3-X1)*ONEDISD
               DY3 =  (Y3-Y1)*ONEDISD
               DZ3 =  (Z3-Z1)*ONEDISD
               R1X =  CX3 + DX3
               R1Y =  CY3 + DY3
               R1Z =  CZ3 + DZ3
               R2X =  CX3 - DX3
               R2Y =  CY3 - DY3
               R2Z =  CZ3 - DZ3
               R3X = -CX3 + DX3
               R3Y = -CY3 + DY3
               R3Z = -CZ3 + DZ3
               R4X =  CX3 + DX3
               R4Y =  CY3 + DY3
               R4Z =  CZ3 + DZ3
               R5X =  CX3
               R5Y =  CY3
               R5Z =  CZ3
               R6X =  CX3
               R6Y =  CY3
               R6Z =  CZ3
               R7X = -CX3
               R7Y = -CY3
               R7Z = -CZ3
               R8X =  CX3
               R8Y =  CY3
               R8Z =  CZ3
               FABX = (R1X*R2*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2X*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3X*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4X*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5X*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6X*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7X*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7*R8X)*ONEFAB*PT5
               FABY = (R1Y*R2*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2Y*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3Y*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4Y*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5Y*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6Y*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7Y*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7*R8Y)*ONEFAB*PT5
               FABZ = (R1Z*R2*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2Z*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3Z*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4Z*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5Z*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6Z*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7Z*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7*R8Z)*ONEFAB*PT5
               TEMP = (FA*FB + FAB)*ONEDISC2*ONEDISC2
               DDISM2X= (X3-X2)*TEMP&
                    &  -((X1-X2)*FA+(X3-X2)*FB)*ONEDISC2&
                    &  -FABX*PT5*ONEDISC2
               DDISM2Y= (Y3-Y2)*TEMP&
                    &  -((Y1-Y2)*FA+(Y3-Y2)*FB)*ONEDISC2&
                    &  -FABY*PT5*ONEDISC2
               DDISM2Z= (Z3-Z2)*TEMP&
                    &  -((Z1-Z2)*FA+(Z3-Z2)*FB)*ONEDISC2&
                    &  -FABZ*PT5*ONEDISC2
               DUM =    (DISM-DISM0)*ONEDISM0*ONEDISM0&
                    &   *PT5*ONEDISM
               IF(DIS13.GE.RB) DUM = -DUM
               DSWF1X = DUM*DDISM2X
               DSWF1Y = DUM*DDISM2Y
               DSWF1Z = DUM*DDISM2Z
               IF(DSWF1X.GT. TEN) DSWF1X =  TEN  ! AVOID SINGULARITY WHEN DISM=0
               IF(DSWF1X.LT.-TEN) DSWF1X = -TEN
               IF(DSWF1Y.GT. TEN) DSWF1Y =  TEN
               IF(DSWF1Y.LT.-TEN) DSWF1Y = -TEN
               IF(DSWF1Z.GT. TEN) DSWF1Z =  TEN
               IF(DSWF1Z.LT.-TEN) DSWF1Z = -TEN
            END IF
         END IF
!
!        - NOTE: LOOP 150 MEANS MULTI-SPHERE SCALING
!                AREA MUST BE SCALED AFTER AREA DERIVATIVES
!
         TMP(1,JFFAT) = AREA*(SWF1*DSWF2X + SWF2*DSWF1X)
         TMP(2,JFFAT) = AREA*(SWF1*DSWF2Y + SWF2*DSWF1Y)
         TMP(3,JFFAT) = AREA*(SWF1*DSWF2Z + SWF2*DSWF1Z)
         IF(JFFAT.GE.2) THEN
            DO KFFAT = 1, JFFAT-1
               TMP(1,KFFAT) = TMP(1,KFFAT)*SWF1*SWF2
               TMP(2,KFFAT) = TMP(2,KFFAT)*SWF1*SWF2
               TMP(3,KFFAT) = TMP(3,KFFAT)*SWF1*SWF2
            ENDDO
         END IF
         AREA = AREA*SWF1*SWF2
         IF(AREA.LT.0.9D-04) THEN   ! 0.9D-04
            AREA = ZERO
            RETURN
         END IF
         DUMMY=ABS(TMP(1,JFFAT))+ABS(TMP(2,JFFAT))+ABS(TMP(3,JFFAT))
         IF(DUMMY.LT.0.9D-05) cycle
         IDTMPTS(41)  = IDTMPTS(41) + 1
         KKK          = IDTMPTS(41)
         IF(KKK.EQ.39) THEN
            ERROR STOP 'FIXPVASWF: IDTMPTS EXCEEDED 40.'
         END IF
         IDTMPTS(KKK) = JFFAT
      end do ! jffat =1,NCENTS
!
!
      IDTMPTS(41)  = IDTMPTS(41) + 1
      KKK          = IDTMPTS(41)
      IDTMPTS(KKK) = IFFAT
      DO KFFAT = 1, NCENTS
         IF(KFFAT.NE.IFFAT) THEN
            TMP(1,IFFAT) = TMP(1,IFFAT) - TMP(1,KFFAT)
            TMP(2,IFFAT) = TMP(2,IFFAT) - TMP(2,KFFAT)
            TMP(3,IFFAT) = TMP(3,IFFAT) - TMP(3,KFFAT)
         END IF
      ENDDO
!
      RETURN
end subroutine FIXPVASWF

end module pelib_cavity_generators
