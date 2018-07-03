program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to build netcdf files with prescribed runoff
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read REGIONAL mask
! 3- Read input file dimensions in first existing file for specified time window
!    (can be 2d or 3d salinity)
! 4- Process all gridT files over specified period
!
! history : - Feb. 2017: version with namelist (N. Jourdain)
!           - Jun. 2018: enable time-dependent files (N. Jourdain)
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

INTEGER                              :: status, dimID_x, dimID_y, mx_GLO, my_GLO, kday, kmonth, kyear, i0, j0, &
&                                       runoff_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidRNF, mz_GLO,    &
&                                       dimID_time_counter, time_ID, dimID_time, fidS, kfmt, mx_tmp, mtime,    &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO, dimID_z, mb,   &
&                                       nfmt, lon_ID, lat_ID, iREG, jREG, tmask_GLO_ID, ntest, ntest2,         &
&                                       fidMSKREG, mx_REG, my_REG, mz_REG, fidMSKIN, tmask_REG_ID, iii, jjj,   &
&                                       kiter, rs, dij, im1, ip1, jm1, jp1, fidM, my_tmp, socoefr_ID,          &
&                                       ai, aj, bi, bj
CHARACTER(LEN=100)                           :: calendar, time_units
CHARACTER(LEN=150)                           :: file_coord, file_in_RNF, file_REG_RNF, file_in_mask_REG, &
&                                               file_in_coord_REG, file_in_gridS, command_str
INTEGER(KIND=4),ALLOCATABLE,DIMENSION(:)     :: list_fmt
INTEGER(KIND=1),ALLOCATABLE,DIMENSION(:,:)   :: tmask_GLO, tmask_REG, missing, tmp_missing
INTEGER(KIND=1),ALLOCATABLE,DIMENSION(:,:,:) :: tmask
REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)      :: nav_lon, nav_lat, nav_lon_REG, nav_lat_REG, socoefr
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)    :: runoff_GLO, runoff_REG, tmp_runoff_REG
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)        :: time
REAL(KIND=8)                                 :: chkland, eps
LOGICAL                                      :: existfile, iout, ll_climato

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
write(file_in_coord_REG,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coordinates_',a,'.nc')

!- name of regional mesh_mask (input) :
write(file_in_mask_REG,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

eps=1.d-9

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

!=================================================================================
! 2- Read REGIONAL mask :
!=================================================================================

write(*,*) 'Reading regional mask in ', TRIM(file_in_mask_REG)
status = NF90_OPEN(TRIM(file_in_mask_REG),0,fidMSKREG); call erreur(status,.TRUE.,"read regional mask") 
!-
status = NF90_INQ_DIMID(fidMSKREG,"z",dimID_z); call erreur(status,.TRUE.,"inq_dimID_z_REG")
status = NF90_INQ_DIMID(fidMSKREG,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_REG")
status = NF90_INQ_DIMID(fidMSKREG,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_REG")
!-
status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_z,len=mz_REG); call erreur(status,.TRUE.,"inq_dim_z_REG")
status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_y,len=my_REG); call erreur(status,.TRUE.,"inq_dim_y_REG")
status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_x,len=mx_REG); call erreur(status,.TRUE.,"inq_dim_x_REG")
!-
ALLOCATE(  tmask(mx_REG,my_REG,mz_REG), tmask_REG(mx_REG,my_REG)  ) 
status = NF90_INQ_VARID(fidMSKREG,"tmask",tmask_REG_ID) ; call erreur(status,.TRUE.,"inq_tmask_REG_ID")
status = NF90_GET_VAR(fidMSKREG,tmask_REG_ID,tmask)     ; call erreur(status,.TRUE.,"getvar_tmask_REG")
tmask_REG(:,:)=tmask(:,:,1)
DEALLOCATE(tmask)
!-
status = NF90_CLOSE(fidMSKREG); call erreur(status,.TRUE.,"end read fidMSKREG")

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

ll_climato = .false.

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
          ll_climato = .true.
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     write(*,*) 'Looking for existence of ', TRIM(file_in_RNF)
     inquire(file=file_in_RNF, exist=existfile)
     if ( existfile ) then; write(*,*) 'BINGO !'; exit ; endif
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading T,S input dimensions in ', TRIM(file_in_RNF)
    status = NF90_OPEN(TRIM(file_in_RNF),0,fidRNF)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidRNF,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidRNF,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidRNF,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")

    status = NF90_INQUIRE_DIMENSION(fidRNF,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidRNF,dimID_x,len=mx_GLO)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidRNF,dimID_y,len=my_GLO)         ; call erreur(status,.TRUE.,"inq_dim_y")

    write(*,*) 'dimensions :', mx_GLO, my_GLO

    ALLOCATE( nav_lon(mx_GLO,my_GLO), nav_lat(mx_GLO,my_GLO) )

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

!- Read tmask in large-scale/global file:
status = NF90_OPEN(TRIM(file_mask_runoff),0,fidMSKIN);     call erreur(status,.TRUE.,"read mask_GLO") 
status = NF90_INQ_DIMID(fidMSKIN,"z",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"depth",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"deptht",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"nav_lev",dimID_z)
status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_z,len=mz_GLO)
ALLOCATE(  tmask(mx_GLO,my_GLO,mz_GLO), tmask_GLO(mx_GLO,my_GLO)  ) 
status = NF90_INQ_VARID(fidMSKIN,"tmask",tmask_GLO_ID); call erreur(status,.TRUE.,"inq_tmask_GLO_ID")
status = NF90_GET_VAR(fidMSKIN,tmask_GLO_ID,tmask);     call erreur(status,.TRUE.,"getvar_tmask_GLO")
tmask_GLO(:,:)=tmask(:,:,1)
DEALLOCATE(tmask)
status = NF90_CLOSE(fidMSKIN);                          call erreur(status,.TRUE.,"end read mask_GLO")

!- create RNF directory :
write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir ',a,'/RNF')
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

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          401 FORMAT(a,'/RNF/runoff_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_REG_RNF,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 ) then
          402 FORMAT(a,'/RNF/runoff_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_REG_RNF,402) TRIM(config_dir), kyear, kmonth, TRIM(config)
        elseif ( nfmt .eq. 195 ) then
          403 FORMAT(a,'/runoff_climato_',a,'.nc')
          write(file_REG_RNF,403) TRIM(config_dir), TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_REG_RNF  >>>> stop'
          stop
        endif

        ALLOCATE( runoff_GLO(mx_GLO,my_GLO,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read GLOBAL runoff :

        write(*,*) 'Reading runoff in ', TRIM(file_in_RNF)
        
        status = NF90_OPEN(TRIM(file_in_RNF),0,fidRNF)              ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidRNF,"time_counter",time_ID)      ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidRNF,"runoff",runoff_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"RUNOFF",runoff_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"berg_melt",runoff_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"Melt",runoff_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRNF,"melt",runoff_ID)
        call erreur(status,.TRUE.,"inq_runoff_ID")
        
        status = NF90_GET_VAR(fidRNF,time_ID,time)                  ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_VAR(fidRNF,runoff_ID,runoff_GLO)          ; call erreur(status,.TRUE.,"getvar_runoff")

        status = NF90_GET_ATT(fidRNF,time_ID,"calendar",calendar)   ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidRNF,time_ID,"units",time_units)    ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidRNF)                                 ; call erreur(status,.TRUE.,"fin_lecture")     

        !----------------------------------------------------------------------
        ! Projection onto regional grid :

        write(*,*) 'start extraction...'
        !ALLOCATE( missing(mx_REG,my_REG) )
        !ALLOCATE( tmp_missing(mx_REG,my_REG) )
        ALLOCATE( runoff_REG(mx_REG,my_REG,mtime) )
        ALLOCATE( tmp_runoff_REG(mx_REG,my_REG,mtime) )
        missing(:,:)=0
        runoff_REG(:,:,:)=0.d0

        do iREG=1,mx_REG
        do jREG=1,my_REG
          iGLO=NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai)
          jGLO=NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
          if ( iGLO .ge. 1 .and. jGLO .ge. 1 ) then
            if ( tmask_GLO(iGLO,jGLO) .eq. 1 ) then
              do l=1,mtime
                runoff_REG(iREG,jREG,l) = runoff_GLO(iGLO,jGLO,l) * tmask_REG(iREG,jREG)
              enddo
            else
              runoff_REG(iREG,jREG,:) = 0.0 
            endif
          else
            runoff_REG(iREG,jREG,:) = 0.0
          endif
        enddo
        enddo
        
        !- Smoothing (i.e. bi-linear interpolation) :
        dij=INT(MAX(ai,aj)*0.5)
        write(*,*) 'Smoothing over plus and minus :', dij, ' pts'
        tmp_runoff_REG(:,:,:)=runoff_REG(:,:,:)
        do iREG=1,mx_REG
        do jREG=1,my_REG
            im1=MAX(iREG-dij,1) ; ip1=MIN(iREG+dij,mx_REG) 
            jm1=MAX(jREG-dij,1) ; jp1=MIN(jREG+dij,my_REG)
            if ( tmask_REG(iREG,jREG) .eq. 1 ) then 
              do l=1,mtime
                tmp_runoff_REG(iREG,jREG,l) =   SUM( SUM( runoff_REG(im1:ip1,jm1:jp1,l) * tmask_REG(im1:ip1,jm1:jp1), 2), 1) &
                &                             / SUM( SUM(                          1.0  * tmask_REG(im1:ip1,jm1:jp1), 2), 1)
              enddo
            else
              tmp_runoff_REG(iREG,jREG,:) = 0.e0
            endif
        enddo
        enddo
        runoff_REG(:,:,:)=tmp_runoff_REG(:,:,:)
        
        !------------------------------------------------------------------------------
        ! Define socoefr to avoid SSS restoring within some distance from the coast:
        ! socoefr=0.5 where no SSS relaxation is applied, and = 0.0 everywhere else.
        ! (and between 0 and 0.5 for partial SSS relaxation).
        
        mb=nn_band
        
        ALLOCATE( socoefr(mx_REG,my_REG) )
        socoefr(:,:) = 0.0
        
        !! manual corrections :
        if (      TRIM(config) == 'AMU12'  &
        &    .or. TRIM(config) == 'AMU12y' &
        &    .or. TRIM(config) == 'AMU12r' ) then
          i0 = imin_ORCA12 - 1811  ! To keep the boxes at the same position even if
          j0 = jmin_ORCA12 -  571  ! nn_jmin_extract & nn_jmax_extract are changed.
          tmask_REG(i0+240:i0+320,j0+1:j0+80) = 0.0 ! to remove SSS restoring in the TG & PIG Bay
        elseif ( TRIM(config) == 'WED12' ) then
          i0 = imin_ORCA12 - 2464  ! To keep the boxes at the same position even if
          j0 = jmin_ORCA12 -  151  ! nn_jmin_extract & nn_jmax_extract are changed.
          tmask_REG(:,j0+1:j0+415) = 0.0  ! region where the old ORCA025 grid was masked
        endif
        
        DO i=1,mx_REG
        DO j=1,my_REG
          ! mask surrounding domain excluded :
          chkland=SUM(SUM(1.00000000000000*(1-tmask_REG(MAX(2,i-mb):MIN(i+mb,mx_REG-1),MAX(2,j-mb):MIN(j+mb,my_REG-1))),2),1)
          if ( chkland .gt. 5.5 ) socoefr(i,j) = 0.5 * MIN(chkland,0.15*4*mb**2) / (0.15*4*mb**2)  ! 0.15 is for 15% of the square with land
        ENDDO
        ENDDO
        
        DO i=1,mx_REG
        DO j=1,my_REG
          if ( tmask_REG(i,j) .eq. 0 ) socoefr(i,j) = 0.5  ! Important to mask where there are ice shelf cavities !
        ENDDO
        ENDDO

        !--------------------------------------
        ! Write RNF netcdf file

        write(*,*) 'Creating ', TRIM(file_REG_RNF)
        status = NF90_CREATE(TRIM(file_REG_RNF),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create RNF file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"x",mx_REG,dimID_x)                       ; call erreur(status,.TRUE.,"def_dimID_x")
        status = NF90_DEF_DIM(fidM,"y",my_REG,dimID_y)                       ; call erreur(status,.TRUE.,"def_dimID_y")
        
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_DOUBLE,(/dimID_time/),time_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"runoff",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time/),runoff_ID)
        call erreur(status,.TRUE.,"def_var_runoff_ID")
        status = NF90_DEF_VAR(fidM,"socoefr",NF90_FLOAT,(/dimID_x,dimID_y/),socoefr_ID)
        call erreur(status,.TRUE.,"def_var_socoefr_ID")
        
        status = NF90_PUT_ATT(fidM,runoff_ID,"associate","time_counter, y, x")  ; call erreur(status,.TRUE.,"put_att_runoff_ID")
        status = NF90_PUT_ATT(fidM,runoff_ID,"missing_value",0.)                ; call erreur(status,.TRUE.,"put_att_runoff_ID")
        status = NF90_PUT_ATT(fidM,runoff_ID,"_FillValue",0.)                   ; call erreur(status,.TRUE.,"put_att_runoff_ID")
        status = NF90_PUT_ATT(fidM,runoff_ID,"units","kg/m2/s")                 ; call erreur(status,.TRUE.,"put_att_runoff_ID")
        status = NF90_PUT_ATT(fidM,runoff_ID,"long_name","runoff")              ; call erreur(status,.TRUE.,"put_att_runoff_ID")
        status = NF90_PUT_ATT(fidM,socoefr_ID,"units","-")                      ; call erreur(status,.TRUE.,"put_att_socoefr_ID")
        status = NF90_PUT_ATT(fidM,socoefr_ID,"short_name","socoefr")           ; call erreur(status,.TRUE.,"put_att_socoefr_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"units",TRIM(time_units))            ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"calendar",TRIM(calendar))           ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"title","Time")                      ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"long_name","Time axis")             ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"standard_name","time")              ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"axis","T")                          ; call erreur(status,.TRUE.,"put_att_time_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_runoff.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"fin_definition") 
        
        status = NF90_PUT_VAR(fidM,time_ID,time)         ; call erreur(status,.TRUE.,"var_time_ID")
        status = NF90_PUT_VAR(fidM,runoff_ID,runoff_REG) ; call erreur(status,.TRUE.,"var_runoff_ID")
        status = NF90_PUT_VAR(fidM,socoefr_ID,socoefr)   ; call erreur(status,.TRUE.,"var_socoefr_ID")       
 
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close new RNF file")

        !--       
        DEALLOCATE( runoff_GLO, time )
        DEALLOCATE( runoff_REG, socoefr )
        DEALLOCATE( tmp_runoff_REG )
        
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
    if (ll_climato) exit
  ENDDO !- kmonth
  if (ll_climato) exit
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
