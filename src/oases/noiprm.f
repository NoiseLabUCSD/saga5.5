      DIMENSION SLS(4),NWS(3),SLD(4),NWD(3)
      DIMENSION CONMAX(3)
      DIMENSION SLEVF(NRD),SLEVDB(NRD)
      DIMENSION ZDN(NSMAX),XDN(NSMAX),YDN(NSMAX),DNLEVDB(NSMAX),
     -          DNLEV(NSMAX)
      COMMON /NOIPRM/ SNLEVDB,WNLEVDB,DPLEVDB,NDNS,
     &                SLEVEL ,WNLEV  ,DPLEV  ,
     &                CMINN  ,CMAXN  ,NWVNON , NWS , SLS ,
     &                CMINP  ,CMAXP  ,NWVNOP , NWD , SLD , DPSD,
     &                CMIND  ,CMAXD  ,NWVNOD , ICUT1D , ICUT2D ,
     &                SLOW1D ,SLOW2D ,CONMAX , SLEVF  ,SLEVDB  ,
     &                ZDN    ,XDN    ,YDN    ,DNLEVDB ,DNLEV
