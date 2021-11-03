program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, IGE-CNRS, May. 2021
!
! purpose: extract internal wave mixing parameters from large-scale grid to regional grid
!
! 0- Initializations
! 1- Read information on grids
! 2- Read zdf_iwm coefficients
! 3- Read REGIONAL mesh/mask file
! 4- Projection onto regional grid
! 5- Writing regional zdf_iwm coefficients file 
!
! history: - May 2021: Creation (N. Jourdain, IGE)
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
&                                           mzreg, mtreg, time_counter_ID, nav_lon_ID, nav_lat_ID, fidM, fidglo, jmin_ORCA12,            &
&                                           i, j, iGLO, jGLO, kk, rr, rs, mxglo, myglo, mzglo, mtglo, fidMSH, l, tmask_ID, mb, kki, kkj, &
&                                           fidC, ai, aj, bi, bj, iREG, jREG, mx_REG, my_REG, imin_ORCA12, dij, im1, ip1,                &
&                                           jm1, jp1, scale_bot_ID, scale_cri_ID, mixing_cri_ID, mixing_pyc_ID, mixing_bot_ID
CHARACTER(LEN=150)                       :: file_zdfiwm_out, file_in_mask_REG, file_in_coord_REG
REAL*4,ALLOCATABLE,DIMENSION(:)          :: time_counter           
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: nav_lon, nav_lat, lonreg, latreg
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: scale_bot, scale_cri, mixing_cri, mixing_pyc, mixing_bot
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: scale_bot_REG, scale_cri_REG, mixing_cri_REG, mixing_pyc_REG, mixing_bot_REG
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: tmp_scale_bot_REG, tmp_scale_cri_REG, tmp_mixing_cri_REG, tmp_mixing_pyc_REG, tmp_mixing_bot_REG
INTEGER*1,ALLOCATABLE,DIMENSION(:,:)     :: tmask_REG
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

!- name of regional mesh_mask (input) :
write(file_in_mask_REG,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

!- name of regional coordinates (input) :
write(file_in_coord_REG,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

!- output file names :
write(file_zdfiwm_out,201)  TRIM(config_dir), TRIM(config)
201 FORMAT(a,'/zdfiwm_',a,'.nc')

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
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"MONTH_REG",dimID_time_counter)
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
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"XAXIS",nav_lon_ID) 
call erreur(status,.TRUE.,"inq_nav_lon_ID")
status = NF90_INQ_VARID(fidA,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"YAXIS",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID")
status = NF90_INQ_VARID(fidA,"scale_bot",scale_bot_ID)  ; call erreur(status,.TRUE.,"inq_scale_bot_ID")
status = NF90_INQ_VARID(fidA,"scale_cri",scale_cri_ID)  ; call erreur(status,.TRUE.,"inq_scale_cri_ID")
status = NF90_INQ_VARID(fidA,"mixing_cri",mixing_cri_ID); call erreur(status,.TRUE.,"inq_mixing_cri_ID")
status = NF90_INQ_VARID(fidA,"mixing_pyc",mixing_pyc_ID); call erreur(status,.TRUE.,"inq_mixing_pyc_ID")
status = NF90_INQ_VARID(fidA,"mixing_bot",mixing_bot_ID); call erreur(status,.TRUE.,"inq_mixing_bot_ID")

status = NF90_GET_VAR(fidA,scale_bot_ID,scale_bot)  ; call erreur(status,.TRUE.,"getvar_scale_bot")
status = NF90_GET_VAR(fidA,scale_cri_ID,scale_cri)  ; call erreur(status,.TRUE.,"getvar_scale_cri")
status = NF90_GET_VAR(fidA,mixing_cri_ID,mixing_cri); call erreur(status,.TRUE.,"getvar_mixing_cri")
status = NF90_GET_VAR(fidA,mixing_pyc_ID,mixing_pyc); call erreur(status,.TRUE.,"getvar_mixing_pyc")
status = NF90_GET_VAR(fidA,mixing_bot_ID,mixing_bot); call erreur(status,.TRUE.,"getvar_mixing_bot")
status = NF90_GET_VAR(fidA,nav_lon_ID,nav_lon)      ; call erreur(status,.TRUE.,"getvar_nav_lon")
status = NF90_GET_VAR(fidA,nav_lat_ID,nav_lat)      ; call erreur(status,.TRUE.,"getvar_nav_lat")

status = NF90_CLOSE(fidA) ; call erreur(status,.TRUE.,"end_reading")     
                                      
!=================================================================================
! 3- Read REGIONAL mesh/mask file
!=================================================================================

write(*,*) 'Reading ', TRIM(file_in_mask_REG)

status = NF90_OPEN(TRIM(file_in_mask_REG),0,fidMSH) ; call erreur(status,.TRUE.,"read AMU12 mesh mask") 

status = NF90_INQ_DIMID(fidMSH,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidMSH,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidMSH,"z",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSH,"nav_lev",dimID_z)
call erreur(status,.TRUE.,"inq_dimID_z")
status = NF90_INQ_DIMID(fidMSH,"t",dimID_t)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSH,"time_counter",dimID_t)
call erreur(status,.TRUE.,"inq_dimID_t")

status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_x,len=mx_REG) ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_y,len=my_REG) ; call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_z,len=mzreg) ; call erreur(status,.TRUE.,"inq_dim_z")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_t,len=mtreg) ; call erreur(status,.TRUE.,"inq_dim_t")

ALLOCATE(  tmask(mx_REG,my_REG,mzreg,mtreg)  )
ALLOCATE(  tmask_REG(mx_REG,my_REG)  )
ALLOCATE(  lonreg(mx_REG,my_REG)  ) 
ALLOCATE(  latreg(mx_REG,my_REG)  ) 
ALLOCATE(  scale_bot_REG(mx_REG,my_REG)  )
ALLOCATE(  scale_cri_REG(mx_REG,my_REG)  )
ALLOCATE(  mixing_cri_REG(mx_REG,my_REG)  )
ALLOCATE(  mixing_pyc_REG(mx_REG,my_REG)  )
ALLOCATE(  mixing_bot_REG(mx_REG,my_REG)  )

status = NF90_INQ_VARID(fidMSH,"tmask",tmask_ID)     ; call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_INQ_VARID(fidMSH,"nav_lon",nav_lon_ID) ; call erreur(status,.TRUE.,"inq_nav_lon_ID")
status = NF90_INQ_VARID(fidMSH,"nav_lat",nav_lat_ID) ; call erreur(status,.TRUE.,"inq_nav_lat_ID")

status = NF90_GET_VAR(fidMSH,tmask_ID,tmask)         ; call erreur(status,.TRUE.,"getvar_tmask")
status = NF90_GET_VAR(fidMSH,nav_lon_ID,lonreg)      ; call erreur(status,.TRUE.,"getvar_nav_lon")
status = NF90_GET_VAR(fidMSH,nav_lat_ID,latreg)      ; call erreur(status,.TRUE.,"getvar_nav_lat")

tmask_REG(:,:)=tmask(:,:,1,1)
DEALLOCATE(tmask)

status = NF90_CLOSE(fidMSH) ; call erreur(status,.TRUE.,"fin_lecture")     

!=================================================================================
! 4- Projection onto regional grid :
!=================================================================================

do iREG=1,mx_REG
do jREG=1,my_REG
  iGLO=NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai)
  jGLO=NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
  if ( iGLO .ge. 1 .and. jGLO .ge. 1 ) then
    scale_bot_REG (iREG,jREG) = scale_bot (iGLO,jGLO) * tmask_REG(iREG,jREG)
    scale_cri_REG (iREG,jREG) = scale_cri (iGLO,jGLO) * tmask_REG(iREG,jREG)
    mixing_cri_REG(iREG,jREG) = mixing_cri(iGLO,jGLO) * tmask_REG(iREG,jREG)
    mixing_pyc_REG(iREG,jREG) = mixing_pyc(iGLO,jGLO) * tmask_REG(iREG,jREG)
    mixing_bot_REG(iREG,jREG) = mixing_bot(iGLO,jGLO) * tmask_REG(iREG,jREG)
  endif
enddo
enddo

write(*,*) lonreg(1,1), latreg(1,1)
iGLO=NINT(FLOAT(1+imin_ORCA12-1-bi)/ai)
jGLO=NINT(FLOAT(1+jmin_ORCA12-1-bj)/aj)
write(*,*) nav_lon(iGLO,jGLO), nav_lat(iGLO,jGLO)

write(*,*) 'max scale_bot : ', MAXVAL(scale_bot), MAXVAL(scale_bot_REG)
write(*,*) 'max scale_cri : ', MAXVAL(scale_cri), MAXVAL(scale_cri_REG)
write(*,*) 'max mixing_pyc : ', MAXVAL(mixing_pyc), MAXVAL(mixing_pyc_REG)

!- Smoothing (i.e. bi-linear interpolation) :
dij=INT(MAX(ai,aj)*0.5)
write(*,*) 'Smoothing over plus and minus :', dij, ' pts'
if ( dij .gt. 0 ) then
  ALLOCATE(  tmp_scale_bot_REG(mx_REG,my_REG)  )
  ALLOCATE(  tmp_scale_cri_REG(mx_REG,my_REG)  )
  ALLOCATE(  tmp_mixing_cri_REG(mx_REG,my_REG)  )
  ALLOCATE(  tmp_mixing_pyc_REG(mx_REG,my_REG)  )
  ALLOCATE(  tmp_mixing_bot_REG(mx_REG,my_REG)  )
  tmp_scale_bot_REG (:,:)=scale_bot_REG (:,:)
  tmp_scale_cri_REG (:,:)=scale_cri_REG (:,:)
  tmp_mixing_cri_REG(:,:)=mixing_cri_REG(:,:)
  tmp_mixing_pyc_REG(:,:)=mixing_pyc_REG(:,:)
  tmp_mixing_bot_REG(:,:)=mixing_bot_REG(:,:)
  do iREG=1,mx_REG
  do jREG=1,my_REG
      im1=MAX(iREG-dij,1) ; ip1=MIN(iREG+dij,mx_REG) 
      jm1=MAX(jREG-dij,1) ; jp1=MIN(jREG+dij,my_REG)
      if ( tmask_REG(iREG,jREG) .eq. 1 ) then 
          tmp_scale_bot_REG (iREG,jREG) =   SUM( SUM( scale_bot_REG (im1:ip1,jm1:jp1) * tmask_REG(im1:ip1,jm1:jp1), 2), 1) &
          &                               / SUM( SUM(                            1.0  * tmask_REG(im1:ip1,jm1:jp1), 2), 1)
          tmp_scale_cri_REG (iREG,jREG) =   SUM( SUM( scale_cri_REG (im1:ip1,jm1:jp1) * tmask_REG(im1:ip1,jm1:jp1), 2), 1) &
          &                               / SUM( SUM(                            1.0  * tmask_REG(im1:ip1,jm1:jp1), 2), 1)
          tmp_mixing_cri_REG(iREG,jREG) =   SUM( SUM( mixing_cri_REG(im1:ip1,jm1:jp1) * tmask_REG(im1:ip1,jm1:jp1), 2), 1) &
          &                               / SUM( SUM(                            1.0  * tmask_REG(im1:ip1,jm1:jp1), 2), 1)
          tmp_mixing_pyc_REG(iREG,jREG) =   SUM( SUM( mixing_pyc_REG(im1:ip1,jm1:jp1) * tmask_REG(im1:ip1,jm1:jp1), 2), 1) &
          &                               / SUM( SUM(                            1.0  * tmask_REG(im1:ip1,jm1:jp1), 2), 1)
          tmp_mixing_bot_REG(iREG,jREG) =   SUM( SUM( mixing_bot_REG(im1:ip1,jm1:jp1) * tmask_REG(im1:ip1,jm1:jp1), 2), 1) &
          &                               / SUM( SUM(                            1.0  * tmask_REG(im1:ip1,jm1:jp1), 2), 1)
      else
        tmp_scale_bot_REG (iREG,jREG) = 0.e0
        tmp_scale_cri_REG (iREG,jREG) = 0.e0
        tmp_mixing_cri_REG(iREG,jREG) = 0.e0
        tmp_mixing_pyc_REG(iREG,jREG) = 0.e0
        tmp_mixing_bot_REG(iREG,jREG) = 0.e0
      endif
  enddo
  enddo
  scale_bot_REG (:,:)=tmp_scale_bot_REG (:,:)
  scale_cri_REG (:,:)=tmp_scale_cri_REG (:,:)
  mixing_cri_REG(:,:)=tmp_mixing_cri_REG(:,:)
  mixing_pyc_REG(:,:)=tmp_mixing_pyc_REG(:,:)
  mixing_bot_REG(:,:)=tmp_mixing_bot_REG(:,:)
endif

!=================================================================================
! 5- Writing new zdfiwm coeffs file on REGIONAL grid
!=================================================================================
                                      
status = NF90_CREATE(TRIM(file_zdfiwm_out),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create new zdfiwm file')
                                        
status = NF90_DEF_DIM(fidM,"x",mx_REG,dimID_x)                  ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"y",my_REG,dimID_y)                  ; call erreur(status,.TRUE.,"def_dimID_y")
                                      
status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lon_ID)      ; call erreur(status,.TRUE.,"def_var_nav_lon_ID")
status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lat_ID)      ; call erreur(status,.TRUE.,"def_var_nav_lat_ID")
status = NF90_DEF_VAR(fidM,"scale_bot",NF90_FLOAT,(/dimID_x,dimID_y/),scale_bot_ID)  ; call erreur(status,.TRUE.,"def_var_scale_bot_ID")
status = NF90_DEF_VAR(fidM,"scale_cri",NF90_FLOAT,(/dimID_x,dimID_y/),scale_cri_ID)  ; call erreur(status,.TRUE.,"def_var_scale_cri_ID")
status = NF90_DEF_VAR(fidM,"mixing_cri",NF90_FLOAT,(/dimID_x,dimID_y/),mixing_cri_ID); call erreur(status,.TRUE.,"def_var_mixing_cri_ID")
status = NF90_DEF_VAR(fidM,"mixing_pyc",NF90_FLOAT,(/dimID_x,dimID_y/),mixing_pyc_ID); call erreur(status,.TRUE.,"def_var_mixing_pyc_ID")
status = NF90_DEF_VAR(fidM,"mixing_bot",NF90_FLOAT,(/dimID_x,dimID_y/),mixing_bot_ID); call erreur(status,.TRUE.,"def_var_mixing_bot_ID")
           
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
status = NF90_PUT_VAR(fidM,scale_bot_ID,scale_bot_REG)   ; call erreur(status,.TRUE.,"var_scale_bot_ID")
status = NF90_PUT_VAR(fidM,scale_cri_ID,scale_cri_REG)   ; call erreur(status,.TRUE.,"var_scale_cri_ID")
status = NF90_PUT_VAR(fidM,mixing_cri_ID,mixing_cri_REG) ; call erreur(status,.TRUE.,"var_mixing_cri_ID")
status = NF90_PUT_VAR(fidM,mixing_pyc_ID,mixing_pyc_REG) ; call erreur(status,.TRUE.,"var_mixing_pyc_ID")
status = NF90_PUT_VAR(fidM,mixing_bot_ID,mixing_bot_REG) ; call erreur(status,.TRUE.,"var_mixing_bot_ID")
                              
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
