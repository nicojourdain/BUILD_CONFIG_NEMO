program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate chlororophyll from parent grid onto child grid
!
! 0- Initializations
! 1- Read information on grids
! 2- Read chlorophyll
! 3- Read child mesh/mask file (CHLD)
! 4- Projection onto child grid
! 5- Writing child chlorophyll file 
!
! history : - Jan. 2015: initial version (N. Jourdain, CNRS-LGGE)
!           - Mar. 2017: GitHub version with namelist (N. Jourdain, CNRS-IGE)
!           - Jan. 2022: cleaning and new naming convention (PAR/CHLD)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /chloro/ file_chloro_in, rn_chla, ln_datelinex, nn_periox
CHARACTER(LEN=50)                        :: config
CHARACTER(LEN=150)                       :: file_chloro_in, config_dir
REAL*4                                   :: rn_chla
LOGICAL                                  :: ln_datelinex
INTEGER                                  :: nn_periox

REAL*4                                   :: minlon, maxlon, minlat, maxlat
INTEGER                                  :: fidA, status, dimID_x, dimID_y, dimID_time_counter, dimID_z, dimID_t, mlon, mlat, mt,        &
&                                           mzreg, mtreg, time_counter_ID, nav_lon_ID, nav_lat_ID, CHLA_ID, fidM, fidglo, ii, jj, kt,    &
&                                           i, j, iPAR, jPAR, kk, rr, rs, mxglo, myglo, mzglo, mtglo, fidMSH, l, tmask_ID, mb, kki, kkj, &
&                                           iCHLD, jCHLD, mx_CHLD, my_CHLD, dij, im1, ip1, jm1, jp1, imin, imax, jmin, jmax, iip1,       &
&                                           iinf, jinf, isup, jsup
CHARACTER(LEN=150)                       :: file_chloro_out, file_in_mask_CHLD, file_in_coord_CHLD
REAL*4,ALLOCATABLE,DIMENSION(:)          :: time_counter           
REAL*4,ALLOCATABLE,DIMENSION(:)          :: lon, lat, zlon, tmp1, tmp2, tmp3, tmp4
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: glamt_CHLD, gphit_CHLD, zglamt_CHLD
REAL*4,ALLOCATABLE,DIMENSION(:,:,:)      :: CHLA, CHLA_CHLD
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
READ (UNIT=1, NML=chloro)
CLOSE(1)

!- name of child mesh_mask (input) :
write(file_in_mask_CHLD,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

!- name of child coordinates (input) :
write(file_in_coord_CHLD,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

!- output file names :
write(file_chloro_out,201)  TRIM(config_dir), TRIM(config)
201 FORMAT(a,'/chlorophyll_',a,'.nc')

!=================================================================================
! 1- Read chlorophyll in (lon,lat) coordinates :
!=================================================================================

write(*,*) 'Reading ', TRIM(file_chloro_in)
 
status = NF90_OPEN(TRIM(file_chloro_in),0,fidA)  ; call erreur(status,.TRUE.,"read input chloro") 
                                   
status = NF90_INQ_DIMID(fidA,"lon",dimID_x)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"x",dimID_x)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"XAXIS",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidA,"lat",dimID_y)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"y",dimID_y)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"YAXIS",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidA,"time_counter",dimID_time_counter)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"time",dimID_time_counter)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"MONTH_CHLD",dimID_time_counter)
call erreur(status,.TRUE.,"inq_dimID_time_counter")
                                       
status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mlon)          ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=mlat)          ; call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_time_counter,len=mt) ; call erreur(status,.TRUE.,"inq_dim_time_counter")
       
ALLOCATE(  time_counter(mt)  ) 
ALLOCATE(  lon(mlon), zlon(mlon)  ) 
ALLOCATE(  lat(mlat)  ) 
ALLOCATE(  CHLA(mlon,mlat,mt)  ) 
         
status = NF90_INQ_VARID(fidA,"time_counter",time_counter_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"MONTH_CHLD",time_counter_ID)
call erreur(status,.TRUE.,"inq_time_counter_ID")
status = NF90_INQ_VARID(fidA,"lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"nav_lon",nav_lon_ID) 
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"longitude",nav_lon_ID) 
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"Longitude",nav_lon_ID) 
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"Lon",nav_lon_ID) 
call erreur(status,.TRUE.,"inq_nav_lon_ID")
status = NF90_INQ_VARID(fidA,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"latitude",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"Latitude",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"Lat",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID")
status = NF90_INQ_VARID(fidA,"CHLA",CHLA_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"chlorophyll",CHLA_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"Chlorophyll",CHLA_ID)
call erreur(status,.TRUE.,"inq_CHLA_ID")
                                      
status = NF90_GET_VAR(fidA,time_counter_ID,time_counter) ; call erreur(status,.TRUE.,"getvar_time_counter")
status = NF90_GET_VAR(fidA,nav_lon_ID,lon)               ; call erreur(status,.TRUE.,"getvar_nav_lon")
status = NF90_GET_VAR(fidA,nav_lat_ID,lat)               ; call erreur(status,.TRUE.,"getvar_nav_lat")
status = NF90_GET_VAR(fidA,CHLA_ID,CHLA)                 ; call erreur(status,.TRUE.,"getvar_CHLA")
                              
status = NF90_CLOSE(fidA) ; call erreur(status,.TRUE.,"fin_lecture")     
                            
zlon(:) = lon(:)
if ( ln_datelinex ) then
  where( lon(:) .lt. 0.e0 )
    zlon(:) = 360.e0 + lon(:)
  endwhere
else
  where( lon(:) .gt. 180.e0 )
    zlon(:) = lon(:) - 360.e0
  endwhere
endif
          
!=================================================================================
! 2- Read child mesh/mask file (CHLD)
!=================================================================================

write(*,*) 'Reading ', TRIM(file_in_mask_CHLD)

status = NF90_OPEN(TRIM(file_in_mask_CHLD),0,fidMSH) ; call erreur(status,.TRUE.,"read child mesh mask") 

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
ALLOCATE(  CHLA_CHLD(mx_CHLD,my_CHLD,12) )
ALLOCATE(  glamt_CHLD(mx_CHLD,my_CHLD), zglamt_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  gphit_CHLD(mx_CHLD,my_CHLD)  ) 

status = NF90_INQ_VARID(fidMSH,"tmask",tmask_ID)     ; call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_INQ_VARID(fidMSH,"glamt",nav_lon_ID)   ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidMSH,"gphit",nav_lat_ID)   ; call erreur(status,.TRUE.,"inq_gphit_ID")

status = NF90_GET_VAR(fidMSH,tmask_ID,tmask)         ; call erreur(status,.TRUE.,"getvar_tmask")
status = NF90_GET_VAR(fidMSH,nav_lon_ID,glamt_CHLD)      ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidMSH,nav_lat_ID,gphit_CHLD)      ; call erreur(status,.TRUE.,"getvar_gphit")

tmask_CHLD(:,:)=tmask(:,:,1,1)
DEALLOCATE(tmask)

status = NF90_CLOSE(fidMSH) ; call erreur(status,.TRUE.,"fin_lecture")     

zglamt_CHLD(:,:)=glamt_CHLD(:,:)
if ( ln_datelinex ) then
  where( glamt_CHLD(:,:) .lt. 0.e0 )
    zglamt_CHLD(:,:)=glamt_CHLD(:,:)+360.e0
  endwhere
else
  where( glamt_CHLD(:,:) .gt. 180.e0 )
    zglamt_CHLD(:,:)=glamt_CHLD(:,:)-360.e0
  endwhere
endif

!- lat/lon CHLD boundaries to limit loops on RTOPO's lon/lat
minlon = MINVAL(zglamt_CHLD-1.0)
maxlon = MAXVAL(zglamt_CHLD+1.0)
minlat = MINVAL( gphit_CHLD-1.0)
maxlat = MAXVAL( gphit_CHLD+1.0)

ALLOCATE( tmp1(mlon), tmp2(mlon) )
ALLOCATE( tmp3(mlat), tmp4(mlat) )
do ii=1,mlon
  tmp1(ii) = abs( zlon(ii) - minlon )
  tmp2(ii) = abs( zlon(ii) - maxlon )
enddo
do jj=1,mlat
  tmp3(jj) = abs(  lat(jj) - minlat )
  tmp4(jj) = abs(  lat(jj) - maxlat )
enddo
imin = MINLOC( tmp1(:), 1 )
imax = MINLOC( tmp2(:), 1 )
jmin = MINLOC( tmp3(:), 1 )
jmax = MINLOC( tmp4(:), 1 )
DEALLOCATE( lon, tmp1, tmp2, tmp3, tmp4 )

write(*,527) imin, imax, jmin, jmax
527 FORMAT('Restricting search on (lon,lat) chlorophyl file to i=',i5,':',i5,' and j=',i5,':',i5)

!=================================================================================
! 3- Projection onto child grid :
!=================================================================================

! Default value :
CHLA_CHLD(:,:,:) = rn_chla


do iCHLD=1,mx_CHLD
do jCHLD=1,my_CHLD

  do ii=imin,imax
    iip1=ii+1
    if ( iip1 .gt. mlon .and. nn_periox .eq. 1 ) then
       iip1=1
    elseif ( iip1 .gt. mlon ) then
       iip1=ii
    endif
    if ( zglamt_CHLD(iCHLD,jCHLD) .ge. zlon(ii) .and. zglamt_CHLD(iCHLD,jCHLD) .lt. zlon(iip1) ) then
      iinf=ii
      isup=iip1
    endif
  enddo

  do jj=jmin,jmax
    if ( gphit_CHLD(iCHLD,jCHLD) .ge. lat(jj) .and. gphit_CHLD(iCHLD,jCHLD) .lt. lat(MIN(jj+1,mlat)) ) then
      jinf=jj
      jsup=MIN(jj+1,mlat)
    endif
  enddo

  CHLA_CHLD(iCHLD,jCHLD,:) = (   CHLA(iinf,jinf,:) * (zlon(isup)-zglamt_CHLD(iCHLD,jCHLD)) * (lat(jsup)-gphit_CHLD(iCHLD,jCHLD))   &
  &                            + CHLA(iinf,jsup,:) * (zlon(isup)-zglamt_CHLD(iCHLD,jCHLD)) * (gphit_CHLD(iCHLD,jCHLD)-lat(jinf))   &
  &                            + CHLA(isup,jinf,:) * (zglamt_CHLD(iCHLD,jCHLD)-zlon(iinf)) * (lat(jsup)-gphit_CHLD(iCHLD,jCHLD))   &
  &                            + CHLA(isup,jsup,:) * (zglamt_CHLD(iCHLD,jCHLD)-zlon(iinf)) * (gphit_CHLD(iCHLD,jCHLD)-lat(jinf)) ) &
  &                        / (                       (zlon(isup)-zglamt_CHLD(iCHLD,jCHLD)) * (lat(jsup)-gphit_CHLD(iCHLD,jCHLD))   &
  &                            +                     (zlon(isup)-zglamt_CHLD(iCHLD,jCHLD)) * (gphit_CHLD(iCHLD,jCHLD)-lat(jinf))   &
  &                            +                     (zglamt_CHLD(iCHLD,jCHLD)-zlon(iinf)) * (lat(jsup)-gphit_CHLD(iCHLD,jCHLD))   &
  &                            +                     (zglamt_CHLD(iCHLD,jCHLD)-zlon(iinf)) * (gphit_CHLD(iCHLD,jCHLD)-lat(jinf)) )

enddo
enddo

!=================================================================================
! 5- Writing new chlorophyll file on CHLDIONAL grid
!=================================================================================
                                      
status = NF90_CREATE(TRIM(file_chloro_out),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create new chloro file')                     
                                        
status = NF90_DEF_DIM(fidM,"x",mx_CHLD,dimID_x)                              ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"y",my_CHLD,dimID_y)                              ; call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
                                      
status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID) ; call erreur(status,.TRUE.,"def_var_time_counter_ID")
status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lon_ID)              ; call erreur(status,.TRUE.,"def_var_nav_lon_ID")
status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lat_ID)              ; call erreur(status,.TRUE.,"def_var_nav_lat_ID")
status = NF90_DEF_VAR(fidM,"CHLA",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),CHLA_ID)
call erreur(status,.TRUE.,"def_var_CHLA_ID")
           
status = NF90_PUT_ATT(fidM,time_counter_ID,"long_name","Time axis")                  ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"title","Time")                           ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
!status = NF90_PUT_ATT(fidM,time_counter_ID,"time_origin","01-JAN-0000 00:00:00")     ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
!status = NF90_PUT_ATT(fidM,time_counter_ID,"units","hour since 0000-01-01 00:00:00") ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
!status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar","noleap")                      ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
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
status = NF90_PUT_ATT(fidM,CHLA_ID,"long_name","Chlorophyll-A")   ; call erreur(status,.TRUE.,"put_att_CHLA_ID")
status = NF90_PUT_ATT(fidM,CHLA_ID,"units","-")                   ; call erreur(status,.TRUE.,"put_att_CHLA_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_chloro_from_lonlat.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
                              
status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"fin_definition") 
                              
status = NF90_PUT_VAR(fidM,time_counter_ID,time_counter) ; call erreur(status,.TRUE.,"var_time_counter_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,glamt_CHLD)          ; call erreur(status,.TRUE.,"var_nav_lon_ID")
status = NF90_PUT_VAR(fidM,nav_lat_ID,gphit_CHLD)          ; call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,CHLA_ID,CHLA_CHLD)            ; call erreur(status,.TRUE.,"var_CHLA_ID")
                              
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
