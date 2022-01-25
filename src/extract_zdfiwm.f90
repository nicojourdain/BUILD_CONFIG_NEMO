program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Extracts internal wave mixing parameters from parent grid to child grid
!
! 0- Initializations
! 1- Read information on grids
! 2- Read zdf_iwm coefficients
! 3- Read CHILD mesh/mask file
! 4- Projection onto child grid
! 5- Writing child zdf_iwm coefficients file 
!
! history: - May  2021: Creation (N. Jourdain, CNRS-IGE)
!          - Jan. 2022: New naming convention (pAR/CHLD/EXT)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /zdfiwm/ file_zdfiwm_in
CHARACTER(LEN=50)                        :: config
CHARACTER(LEN=150)                       :: file_zdfiwm_in, config_dir

INTEGER                                  :: fidA, status, dimID_x, dimID_y, dimID_time_counter, dimID_z, dimID_t, mx, my, mtime_counter, &
&                                           mzreg, mtreg, time_counter_ID, nav_lon_ID, nav_lat_ID, fidM, fidglo, jmin_EXT,            &
&                                           i, j, iPAR, jPAR, kk, rr, rs, mxglo, myglo, mzglo, mtglo, fidMSH, l, tmask_ID, mb, kki, kkj, &
&                                           fidC, ai, aj, bi, bj, iCHLD, jCHLD, mx_CHLD, my_CHLD, imin_EXT, dij, im1, ip1,                &
&                                           jm1, jp1, scale_bot_ID, scale_cri_ID, mixing_cri_ID, mixing_pyc_ID, mixing_bot_ID
CHARACTER(LEN=150)                       :: file_zdfiwm_out, file_in_mask_CHLD, file_in_coord_CHLD
REAL*4,ALLOCATABLE,DIMENSION(:)          :: time_counter           
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: nav_lon, nav_lat, lonreg, latreg
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: scale_bot, scale_cri, mixing_cri, mixing_pyc, mixing_bot
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: scale_bot_CHLD, scale_cri_CHLD, mixing_cri_CHLD, mixing_pyc_CHLD, mixing_bot_CHLD
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: tmp_scale_bot_CHLD, tmp_scale_cri_CHLD, tmp_mixing_cri_CHLD, tmp_mixing_pyc_CHLD, tmp_mixing_bot_CHLD
INTEGER*1,ALLOCATABLE,DIMENSION(:,:)     :: tmask_CHLD
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:,:) :: tmask

!=================================================================================
!- 0- Initialiartions
!=================================================================================

! Default values (replaced with namelist values if specified):
config_dir        = '.'

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=zdfiwm)
CLOSE(1)

!- name of child mesh_mask (input) :
write(file_in_mask_CHLD,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

!- name of child coordinates (input) :
write(file_in_coord_CHLD,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

!- output file names :
write(file_zdfiwm_out,201)  TRIM(config_dir), TRIM(config)
201 FORMAT(a,'/zdfiwm_',a,'.nc')

!=================================================================================
! 1- Read information on grids
!=================================================================================

!- Read EXT attributes of coordinate file to get grid correspondance :
!       i_EXT = ai * i_PAR + bi
!       j_EXT = aj * j_PAR + bj

write(*,*) 'Reading parameters for grid correspondence in the attributes of ', TRIM(file_in_coord_CHLD)
status = NF90_OPEN(TRIM(file_in_coord_CHLD),0,fidC); call erreur(status,.TRUE.,"read EXT grid coord input")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "ai", ai); call erreur(status,.TRUE.,"read att1")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "bi", bi); call erreur(status,.TRUE.,"read att2")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "aj", aj); call erreur(status,.TRUE.,"read att3")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "bj", bj); call erreur(status,.TRUE.,"read att4")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "imin_extraction", imin_EXT); call erreur(status,.TRUE.,"read att5")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "jmin_extraction", jmin_EXT); call erreur(status,.TRUE.,"read att6")
status = NF90_CLOSE(fidC)                         ; call erreur(status,.TRUE.,"end read fidC")

!=================================================================================
! 2- Read fileds for vertical diffusion due to mixing by internal waves :
!=================================================================================

write(*,*) 'Reading ', TRIM(file_zdfiwm_in)
 
status = NF90_OPEN(TRIM(file_zdfiwm_in),0,fidA)  ; call erreur(status,.TRUE.,"read input zdfiwm") 
                                   
status = NF90_INQ_DIMID(fidA,"x",dimID_x)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"XAXIS",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidA,"y",dimID_y)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"YAXIS",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidA,"time_counter",dimID_time_counter)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"MONTH_CHLD",dimID_time_counter)
call erreur(status,.TRUE.,"inq_dimID_time_counter")

status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx)  ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my)  ; call erreur(status,.TRUE.,"inq_dim_y")

ALLOCATE(  nav_lon(mx,my)  ) 
ALLOCATE(  nav_lat(mx,my)  ) 
ALLOCATE(  scale_bot(mx,my)  )
ALLOCATE(  scale_cri(mx,my)  )
ALLOCATE(  mixing_cri(mx,my)  )
ALLOCATE(  mixing_pyc(mx,my)  )
ALLOCATE(  mixing_bot(mx,my)  )

status = NF90_INQ_VARID(fidA,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"glamt",nav_lon_ID) 
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"lon",nav_lon_ID) 
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"XAXIS",nav_lon_ID) 
call erreur(status,.TRUE.,"inq_nav_lon_ID")
status = NF90_INQ_VARID(fidA,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"gphit",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"YAXIS",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID")
!-
status = NF90_INQ_VARID(fidA,"scale_bot",scale_bot_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"decay_scale_bot",scale_bot_ID)
call erreur(status,.TRUE.,"inq_scale_bot_ID")
!-
status = NF90_INQ_VARID(fidA,"scale_cri",scale_cri_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"decay_scale_cri",scale_cri_ID)
call erreur(status,.TRUE.,"inq_scale_cri_ID")
!-
status = NF90_INQ_VARID(fidA,"mixing_cri",mixing_cri_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"mixing_power_cri",mixing_cri_ID)
call erreur(status,.TRUE.,"inq_mixing_cri_ID")
!-
status = NF90_INQ_VARID(fidA,"mixing_pyc",mixing_pyc_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"mixing_power_pyc",mixing_pyc_ID)
call erreur(status,.TRUE.,"inq_mixing_pyc_ID")
!-
status = NF90_INQ_VARID(fidA,"mixing_bot",mixing_bot_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"mixing_power_bot",mixing_bot_ID)
call erreur(status,.TRUE.,"inq_mixing_bot_ID")

status = NF90_GET_VAR(fidA,scale_bot_ID,scale_bot)  ; call erreur(status,.TRUE.,"getvar_scale_bot")
status = NF90_GET_VAR(fidA,scale_cri_ID,scale_cri)  ; call erreur(status,.TRUE.,"getvar_scale_cri")
status = NF90_GET_VAR(fidA,mixing_cri_ID,mixing_cri); call erreur(status,.TRUE.,"getvar_mixing_cri")
status = NF90_GET_VAR(fidA,mixing_pyc_ID,mixing_pyc); call erreur(status,.TRUE.,"getvar_mixing_pyc")
status = NF90_GET_VAR(fidA,mixing_bot_ID,mixing_bot); call erreur(status,.TRUE.,"getvar_mixing_bot")
status = NF90_GET_VAR(fidA,nav_lon_ID,nav_lon)      ; call erreur(status,.TRUE.,"getvar_nav_lon")
status = NF90_GET_VAR(fidA,nav_lat_ID,nav_lat)      ; call erreur(status,.TRUE.,"getvar_nav_lat")

status = NF90_CLOSE(fidA) ; call erreur(status,.TRUE.,"end_reading")     
                                      
!=================================================================================
! 3- Read CHLDIONAL mesh/mask file
!=================================================================================

write(*,*) 'Reading ', TRIM(file_in_mask_CHLD)

status = NF90_OPEN(TRIM(file_in_mask_CHLD),0,fidMSH) ; call erreur(status,.TRUE.,"read AMU12 mesh mask") 

status = NF90_INQ_DIMID(fidMSH,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidMSH,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidMSH,"z",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSH,"nav_lev",dimID_z)
call erreur(status,.TRUE.,"inq_dimID_z")
status = NF90_INQ_DIMID(fidMSH,"t",dimID_t)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSH,"time_counter",dimID_t)
call erreur(status,.TRUE.,"inq_dimID_t")

status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_x,len=mx_CHLD) ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_y,len=my_CHLD) ; call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_z,len=mzreg) ; call erreur(status,.TRUE.,"inq_dim_z")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_t,len=mtreg) ; call erreur(status,.TRUE.,"inq_dim_t")

ALLOCATE(  tmask(mx_CHLD,my_CHLD,mzreg,mtreg)  )
ALLOCATE(  tmask_CHLD(mx_CHLD,my_CHLD)  )
ALLOCATE(  lonreg(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  latreg(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  scale_bot_CHLD(mx_CHLD,my_CHLD)  )
ALLOCATE(  scale_cri_CHLD(mx_CHLD,my_CHLD)  )
ALLOCATE(  mixing_cri_CHLD(mx_CHLD,my_CHLD)  )
ALLOCATE(  mixing_pyc_CHLD(mx_CHLD,my_CHLD)  )
ALLOCATE(  mixing_bot_CHLD(mx_CHLD,my_CHLD)  )

status = NF90_INQ_VARID(fidMSH,"tmask",tmask_ID)     ; call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_INQ_VARID(fidMSH,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidMSH,"glamt",nav_lon_ID)
call erreur(status,.TRUE.,"inq_nav_lon_ID")
status = NF90_INQ_VARID(fidMSH,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidMSH,"gphit",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID")

status = NF90_GET_VAR(fidMSH,tmask_ID,tmask)         ; call erreur(status,.TRUE.,"getvar_tmask")
status = NF90_GET_VAR(fidMSH,nav_lon_ID,lonreg)      ; call erreur(status,.TRUE.,"getvar_nav_lon")
status = NF90_GET_VAR(fidMSH,nav_lat_ID,latreg)      ; call erreur(status,.TRUE.,"getvar_nav_lat")

tmask_CHLD(:,:)=tmask(:,:,1,1)
DEALLOCATE(tmask)

status = NF90_CLOSE(fidMSH) ; call erreur(status,.TRUE.,"fin_lecture")     

!=================================================================================
! 4- Projection onto child grid :
!=================================================================================

do iCHLD=1,mx_CHLD
do jCHLD=1,my_CHLD
  iPAR=NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai)
  jPAR=NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
  if ( iPAR .ge. 1 .and. jPAR .ge. 1 ) then
    scale_bot_CHLD (iCHLD,jCHLD) = scale_bot (iPAR,jPAR) * tmask_CHLD(iCHLD,jCHLD)
    scale_cri_CHLD (iCHLD,jCHLD) = scale_cri (iPAR,jPAR) * tmask_CHLD(iCHLD,jCHLD)
    mixing_cri_CHLD(iCHLD,jCHLD) = mixing_cri(iPAR,jPAR) * tmask_CHLD(iCHLD,jCHLD)
    mixing_pyc_CHLD(iCHLD,jCHLD) = mixing_pyc(iPAR,jPAR) * tmask_CHLD(iCHLD,jCHLD)
    mixing_bot_CHLD(iCHLD,jCHLD) = mixing_bot(iPAR,jPAR) * tmask_CHLD(iCHLD,jCHLD)
  endif
enddo
enddo

write(*,*) lonreg(1,1), latreg(1,1)
iPAR=NINT(FLOAT(1+imin_EXT-1-bi)/ai)
jPAR=NINT(FLOAT(1+jmin_EXT-1-bj)/aj)
write(*,*) nav_lon(iPAR,jPAR), nav_lat(iPAR,jPAR)

write(*,*) 'max scale_bot : ', MAXVAL(scale_bot), MAXVAL(scale_bot_CHLD)
write(*,*) 'max scale_cri : ', MAXVAL(scale_cri), MAXVAL(scale_cri_CHLD)
write(*,*) 'max mixing_pyc : ', MAXVAL(mixing_pyc), MAXVAL(mixing_pyc_CHLD)

!- Smoothing (i.e. bi-linear interpolation) :
dij=INT(MAX(ai,aj)*0.5)
write(*,*) 'Smoothing over plus and minus :', dij, ' pts'
if ( dij .gt. 0 ) then
  ALLOCATE(  tmp_scale_bot_CHLD(mx_CHLD,my_CHLD)  )
  ALLOCATE(  tmp_scale_cri_CHLD(mx_CHLD,my_CHLD)  )
  ALLOCATE(  tmp_mixing_cri_CHLD(mx_CHLD,my_CHLD)  )
  ALLOCATE(  tmp_mixing_pyc_CHLD(mx_CHLD,my_CHLD)  )
  ALLOCATE(  tmp_mixing_bot_CHLD(mx_CHLD,my_CHLD)  )
  tmp_scale_bot_CHLD (:,:)=scale_bot_CHLD (:,:)
  tmp_scale_cri_CHLD (:,:)=scale_cri_CHLD (:,:)
  tmp_mixing_cri_CHLD(:,:)=mixing_cri_CHLD(:,:)
  tmp_mixing_pyc_CHLD(:,:)=mixing_pyc_CHLD(:,:)
  tmp_mixing_bot_CHLD(:,:)=mixing_bot_CHLD(:,:)
  do iCHLD=1,mx_CHLD
  do jCHLD=1,my_CHLD
      im1=MAX(iCHLD-dij,1) ; ip1=MIN(iCHLD+dij,mx_CHLD) 
      jm1=MAX(jCHLD-dij,1) ; jp1=MIN(jCHLD+dij,my_CHLD)
      if ( tmask_CHLD(iCHLD,jCHLD) .eq. 1 ) then 
          tmp_scale_bot_CHLD (iCHLD,jCHLD) =   SUM( SUM( scale_bot_CHLD (im1:ip1,jm1:jp1) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
          &                                  / SUM( SUM(                            1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
          tmp_scale_cri_CHLD (iCHLD,jCHLD) =   SUM( SUM( scale_cri_CHLD (im1:ip1,jm1:jp1) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
          &                                  / SUM( SUM(                            1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
          tmp_mixing_cri_CHLD(iCHLD,jCHLD) =   SUM( SUM( mixing_cri_CHLD(im1:ip1,jm1:jp1) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
          &                                  / SUM( SUM(                            1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
          tmp_mixing_pyc_CHLD(iCHLD,jCHLD) =   SUM( SUM( mixing_pyc_CHLD(im1:ip1,jm1:jp1) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
          &                                  / SUM( SUM(                            1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
          tmp_mixing_bot_CHLD(iCHLD,jCHLD) =   SUM( SUM( mixing_bot_CHLD(im1:ip1,jm1:jp1) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
          &                                  / SUM( SUM(                            1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
      else
        tmp_scale_bot_CHLD (iCHLD,jCHLD) = 0.e0
        tmp_scale_cri_CHLD (iCHLD,jCHLD) = 0.e0
        tmp_mixing_cri_CHLD(iCHLD,jCHLD) = 0.e0
        tmp_mixing_pyc_CHLD(iCHLD,jCHLD) = 0.e0
        tmp_mixing_bot_CHLD(iCHLD,jCHLD) = 0.e0
      endif
  enddo
  enddo
  scale_bot_CHLD (:,:)=tmp_scale_bot_CHLD (:,:)
  scale_cri_CHLD (:,:)=tmp_scale_cri_CHLD (:,:)
  mixing_cri_CHLD(:,:)=tmp_mixing_cri_CHLD(:,:)
  mixing_pyc_CHLD(:,:)=tmp_mixing_pyc_CHLD(:,:)
  mixing_bot_CHLD(:,:)=tmp_mixing_bot_CHLD(:,:)
endif

!=================================================================================
! 5- Writing new zdfiwm coeffs file on CHLDIONAL grid
!=================================================================================
                                      
status = NF90_CREATE(TRIM(file_zdfiwm_out),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create new zdfiwm file')
                                        
status = NF90_DEF_DIM(fidM,"x",mx_CHLD,dimID_x)                  ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"y",my_CHLD,dimID_y)                  ; call erreur(status,.TRUE.,"def_dimID_y")
                                      
status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lon_ID)      ; call erreur(status,.TRUE.,"def_var_nav_lon_ID")
status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lat_ID)      ; call erreur(status,.TRUE.,"def_var_nav_lat_ID")
status = NF90_DEF_VAR(fidM,"decay_scale_bot",NF90_FLOAT,(/dimID_x,dimID_y/),scale_bot_ID)  ; call erreur(status,.TRUE.,"def_var_scale_bot_ID")
status = NF90_DEF_VAR(fidM,"decay_scale_cri",NF90_FLOAT,(/dimID_x,dimID_y/),scale_cri_ID)  ; call erreur(status,.TRUE.,"def_var_scale_cri_ID")
status = NF90_DEF_VAR(fidM,"mixing_power_cri",NF90_FLOAT,(/dimID_x,dimID_y/),mixing_cri_ID); call erreur(status,.TRUE.,"def_var_mixing_cri_ID")
status = NF90_DEF_VAR(fidM,"mixing_power_pyc",NF90_FLOAT,(/dimID_x,dimID_y/),mixing_pyc_ID); call erreur(status,.TRUE.,"def_var_mixing_pyc_ID")
status = NF90_DEF_VAR(fidM,"mixing_power_bot",NF90_FLOAT,(/dimID_x,dimID_y/),mixing_bot_ID); call erreur(status,.TRUE.,"def_var_mixing_bot_ID")
           
status = NF90_PUT_ATT(fidM,nav_lon_ID,"nav_model","Default grid") ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"long_name","Longitude")    ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"valid_max",180.)           ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"valid_min",-180.)          ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")     ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"nav_model","Default grid") ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"long_name","Latitude")     ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"valid_max",90.)            ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"valid_min",-90.)           ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")    ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,scale_bot_ID,"valid_max",6100.)        ; call erreur(status,.TRUE.,"put_att_scale_bot_ID")
status = NF90_PUT_ATT(fidM,scale_bot_ID,"valid_min",100.)         ; call erreur(status,.TRUE.,"put_att_scale_bot_ID")
status = NF90_PUT_ATT(fidM,scale_cri_ID,"valid_max",6100.)        ; call erreur(status,.TRUE.,"put_att_scale_cri_ID")
status = NF90_PUT_ATT(fidM,scale_cri_ID,"valid_min",100.)         ; call erreur(status,.TRUE.,"put_att_scale_cri_ID")
status = NF90_PUT_ATT(fidM,mixing_cri_ID,"valid_max",1.)          ; call erreur(status,.TRUE.,"put_att_mixing_cri_ID")
status = NF90_PUT_ATT(fidM,mixing_cri_ID,"valid_min",1.e-10)      ; call erreur(status,.TRUE.,"put_att_mixing_cri_ID")
status = NF90_PUT_ATT(fidM,mixing_pyc_ID,"valid_max",1.)          ; call erreur(status,.TRUE.,"put_att_mixing_pyc_ID")
status = NF90_PUT_ATT(fidM,mixing_pyc_ID,"valid_min",1.e-10)      ; call erreur(status,.TRUE.,"put_att_mixing_pyc_ID")
status = NF90_PUT_ATT(fidM,mixing_bot_ID,"valid_max",1.)          ; call erreur(status,.TRUE.,"put_att_mixing_bot_ID")
status = NF90_PUT_ATT(fidM,mixing_bot_ID,"valid_min",1.e-10)      ; call erreur(status,.TRUE.,"put_att_mixing_bot_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_zdfiwm.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
                              
status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"fin_definition") 
                              
status = NF90_PUT_VAR(fidM,nav_lon_ID,lonreg)        ; call erreur(status,.TRUE.,"var_nav_lon_ID")
status = NF90_PUT_VAR(fidM,nav_lat_ID,latreg)        ; call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,scale_bot_ID,scale_bot_CHLD)   ; call erreur(status,.TRUE.,"var_scale_bot_ID")
status = NF90_PUT_VAR(fidM,scale_cri_ID,scale_cri_CHLD)   ; call erreur(status,.TRUE.,"var_scale_cri_ID")
status = NF90_PUT_VAR(fidM,mixing_cri_ID,mixing_cri_CHLD) ; call erreur(status,.TRUE.,"var_mixing_cri_ID")
status = NF90_PUT_VAR(fidM,mixing_pyc_ID,mixing_pyc_CHLD) ; call erreur(status,.TRUE.,"var_mixing_pyc_ID")
status = NF90_PUT_VAR(fidM,mixing_bot_ID,mixing_bot_CHLD) ; call erreur(status,.TRUE.,"var_mixing_bot_ID")
                              
status = NF90_CLOSE(fidM) ;call erreur(status,.TRUE.,"final")         

end program modif

!----------------------------

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
