      PARAMETER ( MXTPOP = MAXLBL , MBTPFR = MAXLBL)
      LOGICAL TPAMP,TPALP,ATPOP,BTPOP
      CHARACTER*8 ATPLB,BTPLB
      COMMON /INFTPA/ BTPFR(MBTPFR), NATPOP(8),
     *                NBTPOP(8), NTPCN1(8), NTPCN2(8),
     *                ATPOP(MXTPOP), BTPOP(MXTPOP),
     *                TPAMP,TPALP,IPRTPA,NBTPFR
      COMMON /CHRTPA/ ATPLB(8,MXTPOP),BTPLB(8,MXTPOP)