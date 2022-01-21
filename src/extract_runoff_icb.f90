program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Used to build netcdf files with both prescribed runoff and prescribed iceberg melt
! (if one of these is not provided, the script sets the variable to zero everywhere).
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read child mask (CHLD)
! 3- Read input file dimensions in first existing file for specified time window
!    (can be 2d or 3d salinity)
! 4- Process all gridT files over specified period
!
! history : - Mar. 2015: initial version (N. Jourdain, CNRS-LGGE)
!           - Feb. 2017: version with namelist (N. Jourdain, CNRS-IGE)
!           - Jun. 2018: enable time-dependent files (N. Jourdain)
!           - Nov. 2021: version for iceberg runoff separated from standard runoff (N. Jourdain)
!                        (i.e. works with ln_rnf_icb = .true. or .false. in NEMO4.2)
!           - Jan. 2022: cleaning and new naming convention (EXT/PAR/CHLD)
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
namelist /runoff/ nn_rrr_yeari, nn_rrr_yearf, rrr_dir, rrr_prefix, nn_band, &
&                   rrr_suffix, file_mask_runoff, rrr_sep1, rrr_sep2
CHARACTER(LEN=50)                    :: config, rrr_sep1, rrr_sep2
CHARACTER(LEN=150)                   :: config_dir, rrr_dir, rrr_prefix, rrr_suffix, file_mask_runoff, inputdir, &
&                                       file_in_coord_extract, file_in_bathy_extract, file_in_bathy_bdy, file_in_coord_bdy
INTEGER                              :: nn_rrr_yeari, nn_rrr_yearf, nn_band, nn_imin_extract, nn_imax_extract,   &
&                                       nn_jmin_extract, nn_jmax_extract, nn_perio, nn_isfcav
LOGICAL                              :: ln_dateline

INTEGER                              :: status, dimID_x, dimID_y, mx_PAR, my_PAR, kday, kmonth, kyear, i0, j0, &
&                                       runoff_ID, icb_melt_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidRNF,&
&                                       dimID_time_counter, time_ID, dimID_time, fidS, kfmt, mx_tmp, mtime,    &
&                                       i, j, k, l, fidC, imin_EXT, jmin_EXT, iPAR, jPAR, dimID_z, mb,   &
&                                       nfmt, lon_ID, lat_ID, iCHLD, jCHLD, tmask_PAR_ID, ntest, ntest2, mz_PAR, &
&                                       fidMSKCHLD, mx_CHLD, my_CHLD, mz_CHLD, fidMSKIN, tmask_CHLD_ID, iii, jjj,   &
&                                       kiter, rs, dij, im1, ip1, jm1, jp1, fidM, my_tmp, socoefr_ID, statusb, &
&                                       ai, aj, bi, bj
CHARACTER(LEN=100)                           :: calendar, time_units
CHARACTER(LEN=150)                           :: file_coord, file_in_RNF, file_CHLD_RNF, file_in_mask_CHLD, &
&                                               file_in_coord_CHLD, file_in_gridS, command_str
INTEGER(KIND=4),ALLOCATABLE,DIMENSION(:)     :: list_fmt
INTEGER(KIND=1),ALLOCATABLE,DIMENSION(:,:)   :: tmask_PAR, tmask_CHLD, missing, tmp_missing
INTEGER(KIND=1),ALLOCATABLE,DIMENSION(:,:,:) :: tmask
REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)      :: nav_lon, nav_lat, nav_lon_CHLD, nav_lat_CHLD, socoefr
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)    :: runoff_PAR, runoff_CHLD, tmp_runoff_CHLD, icb_melt_PAR, icb_melt_CHLD, tmp_icb_melt_CHLD
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)        :: time
REAL(KIND=8)                                 :: chkland, eps
LOGICAL                                      :: existfile, iout !, ll_climato

!=================================================================================
!- 0- Initialiartions
!=================================================================================

call gsw_saar_init (.true.)

write(*,*) 'Reading namelist parameters'

! Default values (replaced with namelist values if specified):
config_dir        = '.'

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=runoff)
CLOSE(1)

write(*,*) 'nn_rrr_yeari = ', nn_rrr_yeari

!- name of regional coordinates (input) :
write(file_in_coord_CHLD,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coordinates_',a,'.nc')

!- name of regional mesh_mask (input) :
write(file_in_mask_CHLD,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

eps=1.d-9

!=================================================================================
! 1- Read information on grids
!=================================================================================

!- Read global attributes of coordinate file to get grid correspondance :
!       i_EXT = ai * i_PAR + bi
!       j_EXT = aj * j_PAR + bj

write(*,*) 'Reading parameters for grid correspondence in the attributes of ', TRIM(file_in_coord_CHLD)
status = NF90_OPEN(TRIM(file_in_coord_CHLD),0,fidC); call erreur(status,.TRUE.,"read global grid coord input")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "ai", ai); call erreur(status,.TRUE.,"read att1")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "bi", bi); call erreur(status,.TRUE.,"read att2")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "aj", aj); call erreur(status,.TRUE.,"read att3")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "bj", bj); call erreur(status,.TRUE.,"read att4")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "imin_extraction", imin_EXT); call erreur(status,.TRUE.,"read att5")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "jmin_extraction", jmin_EXT); call erreur(status,.TRUE.,"read att6")
status = NF90_CLOSE(fidC)                         ; call erreur(status,.TRUE.,"end read fidC")

!=================================================================================
! 2- Read child mask (CHLD) :
!=================================================================================

write(*,*) 'Reading regional mask in ', TRIM(file_in_mask_CHLD)
status = NF90_OPEN(TRIM(file_in_mask_CHLD),0,fidMSKCHLD); call erreur(status,.TRUE.,"read regional mask") 
!-
status = NF90_INQ_DIMID(fidMSKCHLD,"z",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKCHLD,"nav_lev",dimID_z)
call erreur(status,.TRUE.,"inq_dimID_z_CHLD")
status = NF90_INQ_DIMID(fidMSKCHLD,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_CHLD")
status = NF90_INQ_DIMID(fidMSKCHLD,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_CHLD")
!-
status = NF90_INQUIRE_DIMENSION(fidMSKCHLD,dimID_z,len=mz_CHLD); call erreur(status,.TRUE.,"inq_dim_z_CHLD")
status = NF90_INQUIRE_DIMENSION(fidMSKCHLD,dimID_y,len=my_CHLD); call erreur(status,.TRUE.,"inq_dim_y_CHLD")
status = NF90_INQUIRE_DIMENSION(fidMSKCHLD,dimID_x,len=mx_CHLD); call erreur(status,.TRUE.,"inq_dim_x_CHLD")
!-
ALLOCATE(  tmask(mx_CHLD,my_CHLD,mz_CHLD), tmask_CHLD(mx_CHLD,my_CHLD)  ) 
status = NF90_INQ_VARID(fidMSKCHLD,"tmask",tmask_CHLD_ID) ; call erreur(status,.TRUE.,"inq_tmask_CHLD_ID")
status = NF90_GET_VAR(fidMSKCHLD,tmask_CHLD_ID,tmask)     ; call erreur(status,.TRUE.,"getvar_tmask_CHLD")
tmask_CHLD(:,:)=tmask(:,:,1)
DEALLOCATE(tmask)
!-
status = NF90_CLOSE(fidMSKCHLD); call erreur(status,.TRUE.,"end read fidMSKCHLD")

!=================================================================================
! 3- Read input file dimensions in first existing file for specified time window
!=================================================================================

!- accepted input format :
191 FORMAT(a,'/',i4.4,'/',a,i4.4,a,i2.2,a,i2.2,a,'.nc')  ! <rrr_dir>/YYYY/<rrr_prefix>YYYY<rrr_sep1>MM<rrr_sep2>DD<rrr_suffix>.nc  
192 FORMAT(a,'/',i4.4,'/',a,i4.4,a,i2.2,a,'.nc')         ! <rrr_dir>/YYYY/<rrr_prefix>YYYY<rrr_sep1>MM<rrr_suffix>.nc  
193 FORMAT(a,'/',a,i4.4,a,i2.2,a,i2.2,a,'.nc')           ! <rrr_dir>/<rrr_prefix>YYYY<rrr_sep1>MM<rrr_sep2>DD<rrr_suffix>.nc  
194 FORMAT(a,'/',a,i4.4,a,i2.2,a,'.nc')                  ! <rrr_dir>/<rrr_prefix>YYYY<rrr_sep1>MM<rrr_suffix>.nc 
195 FORMAT(a,'/',a,'.nc')                                ! <rrr_dir>/<rrr_prefix>.nc

ALLOCATE(list_fmt(5))
list_fmt=(/191,192,193,194,195/)

!ll_climato = .false.

kyear=nn_rrr_yeari
kmonth=1
DO kday=1,31

  do kfmt=1,size(list_fmt)
     nfmt=list_fmt(kfmt)
     SELECT CASE(nfmt)
        CASE(191)
          write(file_in_RNF,191) TRIM(rrr_dir), kyear, TRIM(rrr_prefix), kyear, TRIM(rrr_sep1), kmonth, TRIM(rrr_sep2), kday, TRIM(rrr_suffix)
        CASE(192)
          write(file_in_RNF,192) TRIM(rrr_dir), kyear, TRIM(rrr_prefix), kyear, TRIM(rrr_sep1), kmonth, TRIM(rrr_suffix) 
        CASE(193)
          write(file_in_RNF,193) TRIM(rrr_dir), TRIM(rrr_prefix), kyear, TRIM(rrr_sep1), kmonth, TRIM(rrr_sep2), kday, TRIM(rrr_suffix)
        CASE(194)
          write(file_in_RNF,194) TRIM(rrr_dir), TRIM(rrr_prefix), kyear, TRIM(rrr_sep1), kmonth, TRIM(rrr_suffix)
        CASE(195)
          write(file_in_RNF,195) TRIM(rrr_dir), TRIM(rrr_prefix)
          !ll_climato = .true.
          !write(*,*) 'CLIMATOLOGICAL FILE IS BEING USED'
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     write(*,*) 'Looking for existence of ', TRIM(file_in_RNF)
     inquire(file=file_in_RNF, exist=existfile)
     if ( existfile ) then; write(*,*) 'BINGO !'; exit ; endif
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading runoff input dimensions in ', TRIM(file_in_RNF)
    status = NF90_OPEN(TRIM(file_in_RNF),0,fidRNF)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidRNF,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidRNF,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidRNF,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")

    status = NF90_INQUIRE_DIMENSION(fidRNF,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidRNF,dimID_x,len=mx_PAR)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidRNF,dimID_y,len=my_PAR)         ; call erreur(status,.TRUE.,"inq_dim_y")

    write(*,*) 'dimensions :', mx_PAR, my_PAR

    ALLOCATE( nav_lon(mx_PAR,my_PAR), nav_lat(mx_PAR,my_PAR) )

    status = NF90_INQ_VARID(fidRNF,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidRNF,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
        
    status = NF90_GET_VAR(fidRNF,lon_ID,nav_lon)                    ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidRNF,lat_ID,nav_lat)                    ; call erreur(status,.TRUE.,"getvar_lat")

    status = NF90_CLOSE(fidRNF)                                     ; call erreur(status,.TRUE.,"fin_lecture")

    exit

  ELSEIF ( kday .eq. 31 ) THEN

    write(*,*) 'No RNF file found for first month of ', kyear
    write(*,*) '        >>>>>>>>>>>>>>>>>> stop !!'
    stop

  ENDIF

ENDDO

!- Read tmask in parent file:
status = NF90_OPEN(TRIM(file_mask_runoff),0,fidMSKIN);     call erreur(status,.TRUE.,"read mask_PAR") 
status = NF90_INQ_DIMID(fidMSKIN,"z",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"depth",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"deptht",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"nav_lev",dimID_z)
status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_z,len=mz_PAR)
ALLOCATE(  tmask(mx_PAR,my_PAR,mz_PAR), tmask_PAR(mx_PAR,my_PAR)  ) 
status = NF90_INQ_VARID(fidMSKIN,"tmask",tmask_PAR_ID); call erreur(status,.TRUE.,"inq_tmask_PAR_ID")
status = NF90_GET_VAR(fidMSKIN,tmask_PAR_ID,tmask);     call erreur(status,.TRUE.,"getvar_tmask_PAR")
tmask_PAR(:,:)=tmask(:,:,1)
DEALLOCATE(tmask)
status = NF90_CLOSE(fidMSKIN);                          call erreur(status,.TRUE.,"end read mask_PAR")

!- create RNF directory :
write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir -pv ',a,'/RNF')
CALL system(TRIM(command_str))

!=================================================================================
! 4- Process all gridT files over specified period
!=================================================================================

DO kyear=nn_rrr_yeari,nn_rrr_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_RNF,191) TRIM(rrr_dir), kyear, TRIM(rrr_prefix), kyear, TRIM(rrr_sep1), kmonth, TRIM(rrr_sep2), kday, TRIM(rrr_suffix)
        CASE(192)
          write(file_in_RNF,192) TRIM(rrr_dir), kyear, TRIM(rrr_prefix), kyear, TRIM(rrr_sep1), kmonth, TRIM(rrr_suffix) 
        CASE(193)
          write(file_in_RNF,193) TRIM(rrr_dir), TRIM(rrr_prefix), kyear, TRIM(rrr_sep1), kmonth, TRIM(rrr_sep2), kday, TRIM(rrr_suffix)
        CASE(194)
          write(file_in_RNF,194) TRIM(rrr_dir), TRIM(rrr_prefix), kyear, TRIM(rrr_sep1), kmonth, TRIM(rrr_suffix)
        CASE(195)
          write(file_in_RNF,195) TRIM(rrr_dir), TRIM(rrr_prefix)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_RNF, exist=existfile)
      write(*,*) 'toto ', TRIM(file_in_RNF), existfile

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          401 FORMAT(a,'/RNF/runoff_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_CHLD_RNF,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 ) then
          402 FORMAT(a,'/RNF/runoff_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_CHLD_RNF,402) TRIM(config_dir), kyear, kmonth, TRIM(config)
        elseif ( nfmt .eq. 195 ) then
          403 FORMAT(a,'/runoff_climato_',a,'.nc')
          write(file_CHLD_RNF,403) TRIM(config_dir), TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_CHLD_RNF  >>>> stop'
          stop
        endif

        ALLOCATE( runoff_PAR(mx_PAR,my_PAR,mtime)  )
        ALLOCATE( icb_melt_PAR(mx_PAR,my_PAR,mtime)  )
        ALLOCATE( time(mtime) )
        
        !--------------------------------------------
        ! Read GLOBAL liquid runoff and iceberg melt:

        write(*,*) 'Reading runoff in ', TRIM(file_in_RNF)
        
        status = NF90_OPEN(TRIM(file_in_RNF),0,fidRNF)              ; call erreur(status,.TRUE.,"read EXT TS") 
        
        status = NF90_INQ_VARID(fidRNF,"time_counter",time_ID)      ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_GET_VAR(fidRNF,time_ID,time)                  ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_ATT(fidRNF,time_ID,"calendar",calendar)   ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidRNF,time_ID,"units",time_units)    ; call erreur(status,.TRUE.,"getatt_units")

        status = NF90_INQ_VARID(fidRNF,"runoff",runoff_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"RUNOFF",runoff_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"rnf",runoff_ID)
        if ( status .ne. 0 ) then
          statusb = 1
          write(*,*) '~!@#$%^*  WARNING : liquid runoff is assumed to be zero everywhere  ~!@#$%^*'
          runoff_PAR(:,:,:) = 0.e0
        else
          statusb = 0
          status = NF90_GET_VAR(fidRNF,runoff_ID,runoff_PAR) ; call erreur(status,.TRUE.,"getvar_runoff")
        endif

        status = NF90_INQ_VARID(fidRNF,"berg_melt",icb_melt_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"Melt",icb_melt_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"melt",icb_melt_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"icb_melt",icb_melt_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"iceberg_cea",icb_melt_ID)
        if ( status .ne. 0 ) then
          statusb = statusb + 1
          write(*,*) '~!@#$%^*  WARNING : iceberg melt is assumed to be zero everywhere  ~!@#$%^*'
          icb_melt_PAR(:,:,:) = 0.e0
        else
          status = NF90_GET_VAR(fidRNF,icb_melt_ID,icb_melt_PAR) ; call erreur(status,.TRUE.,"getvar_icb_melt")
        endif

        if ( statusb .eq. 2 ) call erreur(status,.TRUE.,"Found no liquid runoff and no iceberg melt in input file")
        
        status = NF90_CLOSE(fidRNF) ; call erreur(status,.TRUE.,"fin_lecture")     

        !----------------------------------------------------------------------
        ! Projection onto regional grid :

        write(*,*) 'start extraction...'
        ALLOCATE( runoff_CHLD(mx_CHLD,my_CHLD,mtime), icb_melt_CHLD(mx_CHLD,my_CHLD,mtime) )
        ALLOCATE( tmp_runoff_CHLD(mx_CHLD,my_CHLD,mtime), tmp_icb_melt_CHLD(mx_CHLD,my_CHLD,mtime) )
        missing(:,:)=0
        runoff_CHLD(:,:,:)=0.d0
        icb_melt_CHLD(:,:,:)=0.d0

        do iCHLD=1,mx_CHLD
        do jCHLD=1,my_CHLD
          iPAR=NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai)
          jPAR=NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
          if ( iPAR .ge. 1 .and. jPAR .ge. 1 ) then
            if ( tmask_PAR(iPAR,jPAR) .eq. 1 ) then
              do l=1,mtime
                runoff_CHLD(iCHLD,jCHLD,l) = runoff_PAR(iPAR,jPAR,l) * tmask_CHLD(iCHLD,jCHLD)
                icb_melt_CHLD(iCHLD,jCHLD,l) = icb_melt_PAR(iPAR,jPAR,l) * tmask_CHLD(iCHLD,jCHLD)
              enddo
            else
              runoff_CHLD(iCHLD,jCHLD,:) = 0.e0 
              icb_melt_CHLD(iCHLD,jCHLD,:) = 0.e0 
            endif
          else
            runoff_CHLD(iCHLD,jCHLD,:) = 0.e0
            icb_melt_CHLD(iCHLD,jCHLD,:) = 0.e0
          endif
        enddo
        enddo
        
        !- Smoothing (i.e. bi-linear interpolation) :
        dij=INT(MAX(ai,aj)*0.5)
        write(*,*) 'Smoothing over plus and minus :', dij, ' pts'
        tmp_runoff_CHLD(:,:,:)=runoff_CHLD(:,:,:)
        tmp_icb_melt_CHLD(:,:,:)=icb_melt_CHLD(:,:,:)
        do iCHLD=1,mx_CHLD
        do jCHLD=1,my_CHLD
            im1=MAX(iCHLD-dij,1) ; ip1=MIN(iCHLD+dij,mx_CHLD) 
            jm1=MAX(jCHLD-dij,1) ; jp1=MIN(jCHLD+dij,my_CHLD)
            if ( tmask_CHLD(iCHLD,jCHLD) .eq. 1 ) then 
              do l=1,mtime
                tmp_runoff_CHLD(iCHLD,jCHLD,l) =   SUM( SUM( runoff_CHLD(im1:ip1,jm1:jp1,l) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
                &                             / SUM( SUM(                          1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
                tmp_icb_melt_CHLD(iCHLD,jCHLD,l) =   SUM( SUM( icb_melt_CHLD(im1:ip1,jm1:jp1,l) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
                &                               / SUM( SUM(                            1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
              enddo
            else
              tmp_runoff_CHLD(iCHLD,jCHLD,:) = 0.e0
              tmp_icb_melt_CHLD(iCHLD,jCHLD,:) = 0.e0
            endif
        enddo
        enddo
        runoff_CHLD(:,:,:)=tmp_runoff_CHLD(:,:,:)
        icb_melt_CHLD(:,:,:)=tmp_icb_melt_CHLD(:,:,:)
        
        !------------------------------------------------------------------------------
        ! Define socoefr to avoid SSS restoring within some distance from the coast:
        ! socoefr=0.5 where no SSS relaxation is applied, and = 0.0 everywhere else.
        ! (and between 0 and 0.5 for partial SSS relaxation).
        
        mb=nn_band
        
        ALLOCATE( socoefr(mx_CHLD,my_CHLD) )
        socoefr(:,:) = 0.0
        
        !! manual corrections :
        if (      TRIM(config) == 'AMU12'  &
        &    .or. TRIM(config) == 'AMU12y' &
        &    .or. TRIM(config) == 'AMU12r' ) then
          i0 = imin_EXT - 1811  ! To keep the boxes at the same position even if
          j0 = jmin_EXT -  571  ! nn_jmin_extract & nn_jmax_extract are changed.
          tmask_CHLD(i0+240:i0+320,j0+1:j0+80) = 0.0 ! to remove SSS restoring in the TG & PIG Bay
        elseif ( TRIM(config) == 'WED12' ) then
          i0 = imin_EXT - 2464  ! To keep the boxes at the same position even if
          j0 = jmin_EXT -  151  ! nn_jmin_extract & nn_jmax_extract are changed.
          tmask_CHLD(:,j0+1:j0+415) = 0.0  ! region where the old ORCA025 grid was masked
        endif
        
        DO i=1,mx_CHLD
        DO j=1,my_CHLD
          ! mask surrounding domain excluded :
          chkland=SUM(SUM(1.00000000000000*(1-tmask_CHLD(MAX(2,i-mb):MIN(i+mb,mx_CHLD-1),MAX(2,j-mb):MIN(j+mb,my_CHLD-1))),2),1)
          if ( chkland .gt. 5.5 ) socoefr(i,j) = 0.5 * MIN(chkland,0.15*4*mb**2) / (0.15*4*mb**2)  ! 0.15 is for 15% of the square with land
        ENDDO
        ENDDO
        
        DO i=1,mx_CHLD
        DO j=1,my_CHLD
          if ( tmask_CHLD(i,j) .eq. 0 ) socoefr(i,j) = 0.5  ! Important to mask where there are ice shelf cavities !
        ENDDO
        ENDDO

        !--------------------------------------
        ! Write RNF netcdf file

        write(*,*) 'Creating ', TRIM(file_CHLD_RNF)
        status = NF90_CREATE(TRIM(file_CHLD_RNF),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create RNF file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"x",mx_CHLD,dimID_x)                       ; call erreur(status,.TRUE.,"def_dimID_x")
        status = NF90_DEF_DIM(fidM,"y",my_CHLD,dimID_y)                       ; call erreur(status,.TRUE.,"def_dimID_y")
        
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_DOUBLE,(/dimID_time/),time_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"runoff",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time/),runoff_ID)
        call erreur(status,.TRUE.,"def_var_runoff_ID")
        status = NF90_DEF_VAR(fidM,"icb_melt",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time/),icb_melt_ID)
        call erreur(status,.TRUE.,"def_var_icb_melt_ID")
        status = NF90_DEF_VAR(fidM,"socoefr",NF90_FLOAT,(/dimID_x,dimID_y/),socoefr_ID)
        call erreur(status,.TRUE.,"def_var_socoefr_ID")
        
        status = NF90_PUT_ATT(fidM,runoff_ID,"associate","time_counter, y, x")  ; call erreur(status,.TRUE.,"put_att_runoff_ID")
        status = NF90_PUT_ATT(fidM,runoff_ID,"units","kg/m2/s")                 ; call erreur(status,.TRUE.,"put_att_runoff_ID")
        status = NF90_PUT_ATT(fidM,runoff_ID,"long_name","liquid runoff")       ; call erreur(status,.TRUE.,"put_att_runoff_ID")
        status = NF90_PUT_ATT(fidM,icb_melt_ID,"associate","time_counter, y, x"); call erreur(status,.TRUE.,"put_att_icb_melt_ID")
        status = NF90_PUT_ATT(fidM,icb_melt_ID,"units","kg/m2/s")               ; call erreur(status,.TRUE.,"put_att_icb_melt_ID")
        status = NF90_PUT_ATT(fidM,icb_melt_ID,"long_name","iceberg melt")      ; call erreur(status,.TRUE.,"put_att_icb_melt_ID")
        status = NF90_PUT_ATT(fidM,socoefr_ID,"units","-")                      ; call erreur(status,.TRUE.,"put_att_socoefr_ID")
        status = NF90_PUT_ATT(fidM,socoefr_ID,"short_name","socoefr")           ; call erreur(status,.TRUE.,"put_att_socoefr_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"units",TRIM(time_units))            ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"calendar",TRIM(calendar))           ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"title","Time")                      ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"long_name","Time axis")             ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"standard_name","time")              ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"axis","T")                          ; call erreur(status,.TRUE.,"put_att_time_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_runoff_icb.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"fin_definition") 
        
        status = NF90_PUT_VAR(fidM,time_ID,time)             ; call erreur(status,.TRUE.,"var_time_ID")
        status = NF90_PUT_VAR(fidM,runoff_ID,runoff_CHLD)     ; call erreur(status,.TRUE.,"var_runoff_ID")
        status = NF90_PUT_VAR(fidM,icb_melt_ID,icb_melt_CHLD) ; call erreur(status,.TRUE.,"var_icb_melt_ID")
        status = NF90_PUT_VAR(fidM,socoefr_ID,socoefr)       ; call erreur(status,.TRUE.,"var_socoefr_ID")       
 
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close new RNF file")

        !--       
        DEALLOCATE( runoff_PAR, icb_melt_PAR, time )
        DEALLOCATE( runoff_CHLD, icb_melt_CHLD, socoefr )
        DEALLOCATE( tmp_runoff_CHLD, tmp_icb_melt_CHLD )
        
        !--
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          write(*,*) 'Looking for next existing day in this month/year'
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 .or. nfmt .eq. 195 ) then
          write(*,*) 'Only one file per month => switching to next month'
          exit
        else
          write(*,*) 'Keep in mind to decide how to end the main loop according to new file format  >>>> stop'
          stop
        endif

      ENDIF

    ENDDO !-kday
    !if (ll_climato) exit
  ENDDO !- kmonth
  !if (ll_climato) exit
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
