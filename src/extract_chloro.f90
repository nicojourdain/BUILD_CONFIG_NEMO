program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, Feb. 2015
!
! purpose: extract chlororophyll from large-scale grid to regional grid
!
! 0- Initializations
! 1- Read information on grids
! 2- Read chlorophyll
! 3- Read REGIONAL mesh/mask file
! 4- Projection onto regional grid
! 5- Writing regional chlorophyll file 
!
! history: - March 2017: GitHub version with namelist (N. Jourdain, IGE)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /chloro/ file_chloro_in, rn_chla
CHARACTER(LEN=50)                        :: config
CHARACTER(LEN=150)                       :: file_chloro_in, config_dir
REAL*4                                   :: rn_chla

INTEGER                                  :: fidA, status, dimID_x, dimID_y, dimID_time_counter, dimID_z, dimID_t, mx, my, mtime_counter, &
&                                           mzreg, mtreg, time_counter_ID, nav_lon_ID, nav_lat_ID, CHLA_ID, fidM, fidglo, jmin_ORCA12,   &
&                                           i, j, iGLO, jGLO, kk, rr, rs, mxglo, myglo, mzglo, mtglo, fidMSH, l, tmask_ID, mb, kki, kkj, &
&                                           fidC, ai, aj, bi, bj, iREG, jREG, mx_REG, my_REG, imin_ORCA12, dij, im1, ip1,    &
&                                           jm1, jp1
CHARACTER(LEN=150)                       :: file_chloro_out, file_in_mask_REG, file_in_coord_REG
REAL*4,ALLOCATABLE,DIMENSION(:)          :: time_counter           
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: nav_lon, nav_lat, lonreg, latreg
REAL*4,ALLOCATABLE,DIMENSION(:,:,:)      :: CHLA, CHLAreg, tmp_CHLAreg
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
READ (UNIT=1, NML=chloro)
CLOSE(1)

!- name of regional mesh_mask (input) :
write(file_in_mask_REG,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

!- name of regional coordinates (input) :
write(file_in_coord_REG,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

!- output file names :
write(file_chloro_out,201)  TRIM(config_dir), TRIM(config)
201 FORMAT(a,'/chlorophyll_',a,'.nc')

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
! 2- Read iceberg chloro :
!=================================================================================

write(*,*) 'Reading ', TRIM(file_chloro_in)
 
status = NF90_OPEN(TRIM(file_chloro_in),0,fidA)  ; call erreur(status,.TRUE.,"read input chloro") 
                                   
status = NF90_INQ_DIMID(fidA,"x",dimID_x)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"XAXIS",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidA,"y",dimID_y)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"YAXIS",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidA,"time_counter",dimID_time_counter)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"MONTH_REG",dimID_time_counter)
call erreur(status,.TRUE.,"inq_dimID_time_counter")
                                       
status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx)                       ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my)                       ; call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_time_counter,len=mtime_counter) ; call erreur(status,.TRUE.,"inq_dim_time_counter")
       
ALLOCATE(  time_counter(mtime_counter)  ) 
ALLOCATE(  nav_lon(mx,my)  ) 
ALLOCATE(  nav_lat(mx,my)  ) 
ALLOCATE(  CHLA(mx,my,mtime_counter)  ) 
         
status = NF90_INQ_VARID(fidA,"time_counter",time_counter_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"MONTH_REG",time_counter_ID)
call erreur(status,.TRUE.,"inq_time_counter_ID")
status = NF90_INQ_VARID(fidA,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"XAXIS",nav_lon_ID) 
call erreur(status,.TRUE.,"inq_nav_lon_ID")
status = NF90_INQ_VARID(fidA,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"YAXIS",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID")
status = NF90_INQ_VARID(fidA,"CHLA",CHLA_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"chlorophyll",CHLA_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"Chlorophyll",CHLA_ID)
call erreur(status,.TRUE.,"inq_CHLA_ID")
                                      
status = NF90_GET_VAR(fidA,time_counter_ID,time_counter) ; call erreur(status,.TRUE.,"getvar_time_counter")
status = NF90_GET_VAR(fidA,nav_lon_ID,nav_lon)           ; call erreur(status,.TRUE.,"getvar_nav_lon")
status = NF90_GET_VAR(fidA,nav_lat_ID,nav_lat)           ; call erreur(status,.TRUE.,"getvar_nav_lat")
status = NF90_GET_VAR(fidA,CHLA_ID,CHLA)                 ; call erreur(status,.TRUE.,"getvar_CHLA")
                              
status = NF90_CLOSE(fidA) ; call erreur(status,.TRUE.,"fin_lecture")     
                                      
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
ALLOCATE(  CHLAreg(mx_REG,my_REG,12) )
ALLOCATE(  tmp_CHLAreg(mx_REG,my_REG,12) )
ALLOCATE(  lonreg(mx_REG,my_REG)  ) 
ALLOCATE(  latreg(mx_REG,my_REG)  ) 

status = NF90_INQ_VARID(fidMSH,"tmask",tmask_ID)     ; call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_INQ_VARID(fidMSH,"glamt",nav_lon_ID)   ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidMSH,"gphit",nav_lat_ID)   ; call erreur(status,.TRUE.,"inq_gphit_ID")

status = NF90_GET_VAR(fidMSH,tmask_ID,tmask)         ; call erreur(status,.TRUE.,"getvar_tmask")
status = NF90_GET_VAR(fidMSH,nav_lon_ID,lonreg)      ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidMSH,nav_lat_ID,latreg)      ; call erreur(status,.TRUE.,"getvar_gphit")

tmask_REG(:,:)=tmask(:,:,1,1)
DEALLOCATE(tmask)

status = NF90_CLOSE(fidMSH) ; call erreur(status,.TRUE.,"fin_lecture")     

!=================================================================================
! 4- Projection onto regional grid :
!=================================================================================

! Default value :
CHLAreg(:,:,:) = rn_chla

do iREG=1,mx_REG
do jREG=1,my_REG
  iGLO=NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai)
  jGLO=NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
  if ( iGLO .ge. 1 .and. jGLO .ge. 1 ) then
    do l=1,mtime_counter
      CHLAreg(iREG,jREG,l) = CHLA(iGLO,jGLO,l) * tmask_REG(iREG,jREG)
    enddo
  endif
enddo
enddo

!- Smoothing (i.e. bi-linear interpolation) :
dij=INT(MAX(ai,aj)*0.5)
write(*,*) 'Smoothing over plus and minus :', dij, ' pts'
tmp_CHLAreg(:,:,:)=CHLAreg(:,:,:)
do iREG=1,mx_REG
do jREG=1,my_REG
    im1=MAX(iREG-dij,1) ; ip1=MIN(iREG+dij,mx_REG) 
    jm1=MAX(jREG-dij,1) ; jp1=MIN(jREG+dij,my_REG)
    if ( tmask_REG(iREG,jREG) .eq. 1 ) then 
      do l=1,mtime_counter
        tmp_CHLAreg(iREG,jREG,l) =   SUM( SUM( CHLAreg(im1:ip1,jm1:jp1,l) * tmask_REG(im1:ip1,jm1:jp1), 2), 1) &
        &                          / SUM( SUM(                       1.0  * tmask_REG(im1:ip1,jm1:jp1), 2), 1)
      enddo
    else
      tmp_CHLAreg(iREG,jREG,:) = 0.e0
    endif
enddo
enddo
CHLAreg(:,:,:)=tmp_CHLAreg(:,:,:)

!=================================================================================
! 5- Writing new chlorophyll file on REGIONAL grid
!=================================================================================
                                      
status = NF90_CREATE(TRIM(file_chloro_out),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create new chloro file')                     
                                        
status = NF90_DEF_DIM(fidM,"x",mx_REG,dimID_x)                               ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"y",my_REG,dimID_y)                               ; call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
                                      
status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID) ; call erreur(status,.TRUE.,"def_var_time_counter_ID")
status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lon_ID)              ; call erreur(status,.TRUE.,"def_var_nav_lon_ID")
status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lat_ID)              ; call erreur(status,.TRUE.,"def_var_nav_lat_ID")
status = NF90_DEF_VAR(fidM,"CHLA",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),CHLA_ID)
call erreur(status,.TRUE.,"def_var_CHLA_ID")
           
status = NF90_PUT_ATT(fidM,time_counter_ID,"long_name","Time axis")                  ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"title","Time")                           ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"time_origin","01-JAN-0000 00:00:00")     ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"units","hour since 0000-01-01 00:00:00") ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar","noleap")                      ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
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
status = NF90_PUT_ATT(fidM,CHLA_ID,"axis","T")                    ; call erreur(status,.TRUE.,"put_att_CHLA_ID")
status = NF90_PUT_ATT(fidM,CHLA_ID,"online_operation","N/A")      ; call erreur(status,.TRUE.,"put_att_CHLA_ID")
status = NF90_PUT_ATT(fidM,CHLA_ID,"short_name","CHLA")           ; call erreur(status,.TRUE.,"put_att_CHLA_ID")
status = NF90_PUT_ATT(fidM,CHLA_ID,"long_name","seawifs Chlorophyll-A")   ; call erreur(status,.TRUE.,"put_att_CHLA_ID")
status = NF90_PUT_ATT(fidM,CHLA_ID,"units","-")                   ; call erreur(status,.TRUE.,"put_att_CHLA_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_chloro.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
                              
status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"fin_definition") 
                              
status = NF90_PUT_VAR(fidM,time_counter_ID,time_counter) ; call erreur(status,.TRUE.,"var_time_counter_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,lonreg)            ; call erreur(status,.TRUE.,"var_nav_lon_ID")
status = NF90_PUT_VAR(fidM,nav_lat_ID,latreg)            ; call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,CHLA_ID,CHLAreg)              ; call erreur(status,.TRUE.,"var_CHLA_ID")
                              
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
