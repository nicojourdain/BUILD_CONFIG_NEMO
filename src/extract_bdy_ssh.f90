program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to build netcdf coordinate file for BDY
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read input file dimensions in first existing file for specified time window
! 3- Process all gridT files over specified period
!
! history : - Feb. 2017: version with namelist (N. Jourdain)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

use gsw_mod_kinds
use gsw_mod_netcdf
use gsw_mod_toolbox

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /bdy_data/ nn_yeari, nn_yearf, data_dir, data_prefix, nn_bdy_eosmatch, &
&                   data_suffix_T, data_suffix_S, data_suffix_U, data_suffix_V, &
&                   data_suffix_ssh, data_suffix_ice, file_data_mask, file_data_zgr
CHARACTER(LEN=50)                    :: config
CHARACTER(LEN=150)                   :: config_dir, data_dir, data_prefix, data_suffix_T, data_suffix_S, &
&                                       data_suffix_U, data_suffix_V, data_suffix_ssh, data_suffix_ice,  &
&                                       file_data_mask, file_data_zgr
INTEGER                              :: nn_yeari, nn_yearf, nn_bdy_eosmatch

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbt, myb, mxbt, glamt_ID, gphit_ID, &
&                                       e1t_ID, e2t_ID, nbit_ID, nbjt_ID, nbrt_ID, mtime, dimID_x, dimID_y,   &
&                                       mlon, mlat, kday, kmonth, kyear, kbdy, nfmt,                          &
&                                       kt, kz, lon_ID, lat_ID, fidM,                                         &
&                                       sossheig_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidT,           &
&                                       dimID_time_counter, time_ID, dimID_time, fidS,                        &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO,               &
&                                       ai, aj, bi, bj, kfmt
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_SSH, file_bdy_gridSSH,   &
&                                       file_in_coord_REG, file_in_gridS, command_str
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: list_fmt
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbit, nbjt, nbrt
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: glamt, gphit, e1t, e2t, nav_lon, nav_lat, nav_lon_bdy, nav_lat_bdy
REAL*4,ALLOCATABLE,DIMENSION(:,:,:)  :: sossheig, sossheig_bdy
REAL*8,ALLOCATABLE,DIMENSION(:)      :: time
LOGICAL                              :: existfile

!=================================================================================
!- 0- Initialiartions
!=================================================================================

call gsw_saar_init (.true.)

write(*,*) 'Reading namelist parameters'

! Default values (replaced with namelist values if specified):
config_dir        = '.'
nn_bdy_eosmatch   =   1

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=bdy_data)
CLOSE(1)

!- bdy coordinates :
write(file_coord,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/coordinates_bdy_',a,'.nc')

!- name of regional coordinates (input) :
write(file_in_coord_REG,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coordinates_',a,'.nc')

!=================================================================================
! 1- Read information on grids
!=================================================================================

!- Read global attributes of coordinate file to get grid correspondance :
!       i_ORCA12 = ai * i_ORCA025 + bi
!       j_ORCA12 = aj * j_ORCA025 + bj

write(*,*) 'Reading parameters for grid correspondence in the attributes of ', TRIM(file_in_coord_REG)
status = NF90_OPEN(TRIM(file_in_coord_REG),0,fidC); call erreur(status,.TRUE.,"read global grid coord input")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "ai", ai); call erreur(status,.TRUE.,"read att1")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "bi", bi); call erreur(status,.TRUE.,"read att2")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "aj", aj); call erreur(status,.TRUE.,"read att3")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "bj", bj); call erreur(status,.TRUE.,"read att4")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "imin_extraction", imin_ORCA12); call erreur(status,.TRUE.,"read att5")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "jmin_extraction", jmin_ORCA12); call erreur(status,.TRUE.,"read att6")
status = NF90_CLOSE(fidC)                         ; call erreur(status,.TRUE.,"end read fidC")

!- Read BDY coordinates 

write(*,*) 'Reading BDY coordinates in ', TRIM(file_coord)
status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read bdy coordinate file") 

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbt",dimID_xbt) ; call erreur(status,.TRUE.,"inq_dimID_xbt")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbt,len=mxbt) ; call erreur(status,.TRUE.,"inq_dim_xbt")

ALLOCATE(  glamt(mxbt,myb)  ) 
ALLOCATE(  gphit(mxbt,myb)  ) 
ALLOCATE(  e1t(mxbt,myb)  ) 
ALLOCATE(  e2t(mxbt,myb)  ) 
ALLOCATE(  nbit(mxbt,myb)  ) 
ALLOCATE(  nbjt(mxbt,myb)  ) 
ALLOCATE(  nbrt(mxbt,myb)  ) 

status = NF90_INQ_VARID(fidCOORD,"glamt",glamt_ID) ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidCOORD,"gphit",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidCOORD,"e1t",e1t_ID)     ; call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidCOORD,"e2t",e2t_ID)     ; call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidCOORD,"nbit",nbit_ID)   ; call erreur(status,.TRUE.,"inq_nbit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjt",nbjt_ID)   ; call erreur(status,.TRUE.,"inq_nbjt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrt",nbrt_ID)   ; call erreur(status,.TRUE.,"inq_nbrt_ID")

status = NF90_GET_VAR(fidCOORD,glamt_ID,glamt) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidCOORD,gphit_ID,gphit) ; call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidCOORD,e1t_ID,e1t)     ; call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidCOORD,e2t_ID,e2t)     ; call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidCOORD,nbit_ID,nbit)   ; call erreur(status,.TRUE.,"getvar_nbit")
status = NF90_GET_VAR(fidCOORD,nbjt_ID,nbjt)   ; call erreur(status,.TRUE.,"getvar_nbjt")
status = NF90_GET_VAR(fidCOORD,nbrt_ID,nbrt)   ; call erreur(status,.TRUE.,"getvar_nbrt")

status = NF90_CLOSE(fidCOORD) ; call erreur(status,.TRUE.,"close coordinate file")     

!=================================================================================
! 2- Read input file dimensions in first existing file for specified time window
!=================================================================================

!- accepted input format :
191 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')  ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD_<data_suffix_ssh>.nc
192 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')           ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_<data_suffix_ssh>.nc
193 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')        ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD.nc
194 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'.nc')                 ! <data_dir>/YYYY/<data_prefix>_YYYY_MM.nc
195 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')           ! <data_dir>/<data_prefix>_YYYY_MM_DD_<data_suffix_ssh>.nc
196 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')                    ! <data_dir>/<data_prefix>_YYYY_MM_<data_suffix_ssh>.nc
197 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')                 ! <data_dir>/<data_prefix>_YYYY_MM_DD.nc
198 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'.nc')                          ! <data_dir>/<data_prefix>_YYYY_MM.nc
291 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')          ! <data_dir>/YYYY/<data_prefix>_YYYYMMDD_<data_suffix_ssh>.nc
292 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'_',a,'.nc')               ! <data_dir>/YYYY/<data_prefix>_YYYYMM_<data_suffix_ssh>.nc
293 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'.nc')                ! <data_dir>/YYYY/<data_prefix>_YYYYMMDD.nc
294 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'.nc')                     ! <data_dir>/YYYY/<data_prefix>_YYYYMM.nc
295 FORMAT(a,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')                   ! <data_dir>/<data_prefix>_YYYYMMDD_<data_suffix_ssh>.nc
296 FORMAT(a,'/',a,'_',i4.4,i2.2,'_',a,'.nc')                        ! <data_dir>/<data_prefix>_YYYYMM_<data_suffix_ssh>.nc
297 FORMAT(a,'/',a,'_',i4.4,i2.2,i2.2,'.nc')                         ! <data_dir>/<data_prefix>_YYYYMMDD.nc
298 FORMAT(a,'/',a,'_',i4.4,i2.2,'.nc')                              ! <data_dir>/<data_prefix>_YYYYMM.nc

kyear=nn_yeari
kmonth=1
DO kday=1,31

  ALLOCATE(list_fmt(16))
  list_fmt=(/191,192,193,194,195,196,197,198,291,292,293,294,295,296,297,298/)

  do kfmt=1,size(list_fmt)
     nfmt=list_fmt(kfmt)
     SELECT CASE(nfmt)
        CASE(191)
          write(file_in_SSH,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_ssh)
        CASE(192)
          write(file_in_SSH,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_ssh) 
        CASE(193)
          write(file_in_SSH,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_SSH,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_SSH,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_ssh)
        CASE(196) 
          write(file_in_SSH,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_ssh)
        CASE(197) 
          write(file_in_SSH,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_SSH,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE(291)
          write(file_in_SSH,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_ssh)
        CASE(292)
          write(file_in_SSH,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_ssh) 
        CASE(293)
          write(file_in_SSH,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_SSH,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_SSH,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_ssh)
        CASE(296) 
          write(file_in_SSH,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_ssh)
        CASE(297) 
          write(file_in_SSH,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_SSH,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     inquire(file=file_in_SSH, exist=existfile)
     if ( existfile ) exit
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading T,S input dimensions in ', TRIM(file_in_SSH)
    status = NF90_OPEN(TRIM(file_in_SSH),0,fidT)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidT,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidT,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidT,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")

    status = NF90_INQUIRE_DIMENSION(fidT,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidT,dimID_x,len=mlon)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidT,dimID_y,len=mlat)         ; call erreur(status,.TRUE.,"inq_dim_y")

    ALLOCATE( nav_lon(mlon,mlat), nav_lat(mlon,mlat) )

    status = NF90_INQ_VARID(fidT,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidT,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
        
    status = NF90_GET_VAR(fidT,lon_ID,nav_lon)                    ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidT,lat_ID,nav_lat)                    ; call erreur(status,.TRUE.,"getvar_lat")

    status = NF90_CLOSE(fidT)                                     ; call erreur(status,.TRUE.,"fin_lecture")

    exit

  ELSEIF ( kday .eq. 31 ) THEN

    write(*,*) 'No SSH file found for first month of ', kyear
    write(*,*) '        >>>>>>>>>>>>>>>>>> stop !!'
    stop

  ENDIF

ENDDO

!--

ALLOCATE( nav_lon_bdy(mxbt,1), nav_lat_bdy(mxbt,1) )

do kbdy=1,mxbt
  iGLO=NINT(FLOAT(nbit(kbdy,1)+imin_ORCA12-1-bi)/ai)
  jGLO=NINT(FLOAT(nbjt(kbdy,1)+jmin_ORCA12-1-bj)/aj)
  nav_lon_bdy(kbdy,1) = nav_lon(iGLO,jGLO)
  nav_lat_bdy(kbdy,1) = nav_lat(iGLO,jGLO)
enddo

!--

write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir ',a,'/BDY')
CALL system(TRIM(command_str))

!=================================================================================
! 3- Process all gridT files over specified period
!=================================================================================

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_SSH,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_ssh)
        CASE(192)
          write(file_in_SSH,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_ssh) 
        CASE(193)
          write(file_in_SSH,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_SSH,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_SSH,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_ssh)
        CASE(196) 
          write(file_in_SSH,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_ssh)
        CASE(197) 
          write(file_in_SSH,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_SSH,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE(291)
          write(file_in_SSH,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_ssh)
        CASE(292)
          write(file_in_SSH,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_ssh) 
        CASE(293)
          write(file_in_SSH,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_SSH,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_SSH,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_ssh)
        CASE(296) 
          write(file_in_SSH,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_ssh)
        CASE(297) 
          write(file_in_SSH,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_SSH,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_SSH, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 .or. nfmt .eq. 195 .or. nfmt .eq. 197 &
        &   .or. nfmt .eq. 291 .or. nfmt .eq. 293 .or. nfmt .eq. 295 .or. nfmt .eq. 297 ) then
          401 FORMAT(a,'/BDY/bdyT_ssh_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridSSH,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 .or. nfmt .eq. 196 .or. nfmt .eq. 198 &
        &   .or. nfmt .eq. 292 .or. nfmt .eq. 294 .or. nfmt .eq. 296 .or. nfmt .eq. 298 ) then
          402 FORMAT(a,'/BDY/bdyT_ssh_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridSSH,402) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_bdy_gridSSH  >>>> stop'
          stop
        endif

        ALLOCATE( sossheig(mlon,mlat,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read input SSH :

        write(*,*) 'Reading SSH in ', TRIM(file_in_SSH)
        
        status = NF90_OPEN(TRIM(file_in_SSH),0,fidT)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidT,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidT,"sossheig",sossheig_ID)          ; call erreur(status,.TRUE.,"inq_sossheig_ID")
        
        status = NF90_GET_VAR(fidT,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_VAR(fidT,sossheig_ID,sossheig)              ; call erreur(status,.TRUE.,"getvar_sossheig")

        status = NF90_GET_ATT(fidT,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidT,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidT)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Remove possible NaNs :
 
        do i=1,mlon
        do j=1,mlat
        do l=1,mtime
          if ( .not. sossheig(i,j,l) .lt. 50.0 ) then
            sossheig(i,j,l) = 0.0
          endif
        enddo
        enddo
        enddo

        !---------------------------------------
        ! Fill values on bdyT :
      
        ALLOCATE( sossheig_bdy(mxbt,1,mtime)  )
 
        do kbdy=1,mxbt
          iGLO=NINT(FLOAT(nbit(kbdy,1)+imin_ORCA12-1-bi)/ai)
          jGLO=NINT(FLOAT(nbjt(kbdy,1)+jmin_ORCA12-1-bj)/aj)
          do kt=1,mtime
            sossheig_bdy(kbdy,1,kt) = sossheig( iGLO, jGLO, kt )
          enddo
        enddo

        !--------------------------------------
        ! Write BDY netcdf file

        write(*,*) 'Creating ', TRIM(file_bdy_gridSSH)
        status = NF90_CREATE(TRIM(file_bdy_gridSSH),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidM,"xbT",mxbT,dimID_xbT)                             ; call erreur(status,.TRUE.,"def_dimID_xbT")

        status = NF90_DEF_VAR(fidM,"sossheig",NF90_FLOAT,(/dimID_xbT,dimID_yb,dimID_time_counter/),sossheig_ID)
        call erreur(status,.TRUE.,"def_var_sossheig_ID")
        status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_xbT,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_xbT,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"nbrdta",NF90_INT,(/dimID_xbT,dimID_yb/),nbrt_ID)
        call erreur(status,.TRUE.,"def_var_nbrt_ID")
        status = NF90_DEF_VAR(fidM,"nbjdta",NF90_INT,(/dimID_xbT,dimID_yb/),nbjt_ID)
        call erreur(status,.TRUE.,"def_var_nbjt_ID")
        status = NF90_DEF_VAR(fidM,"nbidta",NF90_INT,(/dimID_xbT,dimID_yb/),nbit_ID)
        call erreur(status,.TRUE.,"def_var_nbit_ID")
       
        status = NF90_PUT_ATT(fidM,sossheig_ID,"long_name","Sea Surface Height")    ; call erreur(status,.TRUE.,"put_att_sossheig_ID")
        status = NF90_PUT_ATT(fidM,sossheig_ID,"units","m")                         ; call erreur(status,.TRUE.,"put_att_sossheig_ID")
        status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,nbrt_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
        status = NF90_PUT_ATT(fidM,nbrt_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
        status = NF90_PUT_ATT(fidM,nbjt_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
        status = NF90_PUT_ATT(fidM,nbjt_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
        status = NF90_PUT_ATT(fidM,nbit_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbit_ID")
        status = NF90_PUT_ATT(fidM,nbit_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbit_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bdy_gridT.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidM,sossheig_ID,sossheig_bdy) ; call erreur(status,.TRUE.,"var_sossheig_ID")
        status = NF90_PUT_VAR(fidM,nav_lat_ID,nav_lat_bdy)   ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidM,nav_lon_ID,nav_lon_bdy)   ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidM,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidM,nbrt_ID,nbrt)             ; call erreur(status,.TRUE.,"var_nbrt_ID")
        status = NF90_PUT_VAR(fidM,nbjt_ID,nbjt)             ; call erreur(status,.TRUE.,"var_nbjt_ID")
        status = NF90_PUT_VAR(fidM,nbit_ID,nbit)             ; call erreur(status,.TRUE.,"var_nbit_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close BDY file")
        
        !--       
        DEALLOCATE( sossheig, time )
        DEALLOCATE( sossheig_bdy )

        !--
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 .or. nfmt .eq. 195 .or. nfmt .eq. 197 &
        &   .or. nfmt .eq. 291 .or. nfmt .eq. 293 .or. nfmt .eq. 295 .or. nfmt .eq. 297 ) then
          write(*,*) 'Looking for next existing day in this month/year'
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 .or. nfmt .eq. 196 .or. nfmt .eq. 198 &
        &   .or. nfmt .eq. 292 .or. nfmt .eq. 294 .or. nfmt .eq. 296 .or. nfmt .eq. 298 ) then
          write(*,*) 'Only one file per month => switching to next month'
          exit
        else
          write(*,*) 'Keep in mind to decide how to end the main loop according to new file format  >>>> stop'
          stop
        endif

      ENDIF

    ENDDO !-kday
  ENDDO !- kmonth
ENDDO !- kyear

end program modif

!=============================================

SUBROUTINE erreur(iret, lstop, chaine)
  ! pour les messages d'erreur
  USE netcdf
  INTEGER, INTENT(in)                     :: iret
  LOGICAL, INTENT(in)                     :: lstop
  CHARACTER(LEN=*), INTENT(in)            :: chaine
  !
  CHARACTER(LEN=80)                       :: message
  !
  IF ( iret .NE. 0 ) THEN
    WRITE(*,*) 'ROUTINE: ', TRIM(chaine)
    WRITE(*,*) 'ERREUR: ', iret
    message=NF90_STRERROR(iret)
    WRITE(*,*) 'CA VEUT DIRE:',TRIM(message)
    IF ( lstop ) STOP
  ENDIF
  !
END SUBROUTINE erreur
