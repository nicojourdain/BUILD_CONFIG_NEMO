program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Script to extract the sea ice initial state from the parent grid
!
! To be used with nn_iceini_file = 1 in NEMO4.2 (sea ice initial state from
! single category or multi-category-averaged sea ice data).
!
! 1- Read PARENT ("_PAR") mask
! 2- Read CHILD ("_CHLD") mask
! 3- Read PARENT sea ice initial state
! 4- Extract CHILD from PARENT
! 5- Writing CHILD initial state for sea ice
!
! History: - Nov. 2021: initial version (N. Jourdain, CNRS-IGE)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /init/ file_in_mask_extract, file_in_T, file_in_S, nn_eosmatch, nn_iter, nn_rsmax, nn_rzmax, &
&               rn_temp, rn_sal, nn_smooth, file_in_SI
INTEGER                               :: nn_iter, nn_rsmax, nn_rzmax, nn_eosmatch, nn_smooth
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: file_in_mask_extract, config_dir, file_in_T, file_in_S, file_in_SI
REAL(KIND=4)                          :: rn_temp, rn_sal

INTEGER :: fidMSKIN, fidMSKCHLD, status, dimID_z, dimID_y, dimID_x, my_PAR, mx_PAR, tmask_PAR_ID, &
&          mx_CHLD, my_CHLD, tmask_CHLD_ID, fidSAL, fidTEMP, siconc_ID, sithic_ID, snthic_ID,       &
&          dimID_time_counter, ai, bi, aj, bj, iii, jjj, iPAR, jPAR, iCHLD, jCHLD, fidTin, fidSin, &
&          kiter, rs, time_counter_ID, fidCOORD, imin_EXT, jmin_EXT, lon_ID, lat_ID, dij,  &
&          dep_ID, kPAR, ntest, im1, ip1, jm1, jp1

CHARACTER(LEN=180) :: file_in_mask_CHLD, file_in_coord_CHLD, file_out_SI, file_out_sal 

INTEGER*1,ALLOCATABLE,DIMENSION(:,:) :: tmask_PAR, tmask_CHLD, missing, tmp_missing

REAL(KIND=4),ALLOCATABLE,DIMENSION(:) ::  dep_PAR

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: lon_PAR, lat_PAR

REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: siconc_PAR, sithic_PAR, siconc_CHLD, sithic_CHLD, tmp_siconc_CHLD, &
&                                          tmp_sithic_CHLD, snthic_PAR, snthic_CHLD, tmp_snthic_CHLD

LOGICAL :: iout

!=================================================================================
! 0- Initializations 
!=================================================================================

call gsw_saar_init (.true.)

! Default values (replaced with namelist values if specified):
config_dir        = '.'
nn_iter           = 100
nn_rsmax          =   5
nn_rzmax          =   1
nn_eosmatch       =   1
nn_smooth         =   1
file_in_T         = 'NOT USED'
file_in_S         = 'NOT USED'

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=init)
CLOSE(1)

!- name of regional mesh_mask (input) :
write(file_in_mask_CHLD,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

!- name of regional coordinates (input) :
write(file_in_coord_CHLD,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

!- output file names :
write(file_out_SI,201)  TRIM(config_dir), TRIM(config)
201 FORMAT(a,'/istate_sea_ice_',a,'.nc')

!=================================================================================
! 1- Read PARENT mask :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_in_mask_extract),0,fidMSKIN); call erreur(status,.TRUE.,"read mask input") 

status = NF90_INQ_DIMID(fidMSKIN,"z",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"nav_lev",dimID_z)
call erreur(status,.TRUE.,"inq_dimID_z_PAR")
status = NF90_INQ_DIMID(fidMSKIN,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_PAR")
status = NF90_INQ_DIMID(fidMSKIN,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_PAR")

status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_y,len=my_PAR); call erreur(status,.TRUE.,"inq_dim_y_PAR")
status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_x,len=mx_PAR); call erreur(status,.TRUE.,"inq_dim_x_PAR")

ALLOCATE(  tmask_PAR(mx_PAR,my_PAR)  ) 
ALLOCATE(  lon_PAR  (mx_PAR,my_PAR)  )
ALLOCATE(  lat_PAR  (mx_PAR,my_PAR)  )

status = NF90_INQ_VARID(fidMSKIN,"tmask",tmask_PAR_ID); call erreur(status,.TRUE.,"inq_tmask_PAR_ID")
status = NF90_INQ_VARID(fidMSKIN,"nav_lon",lon_ID)    ; call erreur(status,.TRUE.,"inq_lon_PAR_ID")
status = NF90_INQ_VARID(fidMSKIN,"nav_lat",lat_ID)    ; call erreur(status,.TRUE.,"inq_lat_PAR_ID")

status = NF90_GET_VAR(fidMSKIN,tmask_PAR_ID,tmask_PAR,start=(/1,1,1/),count=(/mx_PAR,my_PAR,1/))
call erreur(status,.TRUE.,"getvar_tmask_PAR")
status = NF90_GET_VAR(fidMSKIN,lon_ID,lon_PAR)        ; call erreur(status,.TRUE.,"getvar_lon_PAR")
status = NF90_GET_VAR(fidMSKIN,lat_ID,lat_PAR)        ; call erreur(status,.TRUE.,"getvar_lat_PAR")

status = NF90_CLOSE(fidMSKIN); call erreur(status,.TRUE.,"end read mask_PAR")

!=================================================================================
! 2- Read CHILD mask :
!=================================================================================

status = NF90_OPEN(TRIM(file_in_mask_CHLD),0,fidMSKCHLD); call erreur(status,.TRUE.,"read regional mask") 

status = NF90_INQ_DIMID(fidMSKCHLD,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_CHLD")
status = NF90_INQ_DIMID(fidMSKCHLD,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_CHLD")

status = NF90_INQUIRE_DIMENSION(fidMSKCHLD,dimID_y,len=my_CHLD); call erreur(status,.TRUE.,"inq_dim_y_CHLD")
status = NF90_INQUIRE_DIMENSION(fidMSKCHLD,dimID_x,len=mx_CHLD); call erreur(status,.TRUE.,"inq_dim_x_CHLD")

ALLOCATE(  tmask_CHLD(mx_CHLD,my_CHLD)  ) 

status = NF90_INQ_VARID(fidMSKCHLD,"tmask",tmask_CHLD_ID); call erreur(status,.TRUE.,"inq_tmask_CHLD_ID")

status = NF90_GET_VAR(fidMSKCHLD,tmask_CHLD_ID,tmask_CHLD,start=(/1,1,1/),count=(/mx_CHLD,my_CHLD,1/))
call erreur(status,.TRUE.,"getvar_tmask_CHLD")

status = NF90_CLOSE(fidMSKCHLD); call erreur(status,.TRUE.,"end read fidMSKCHLD")

!=================================================================================
! 3- Read PARENT sea ice initial state
!=================================================================================

write(*,*) 'Reading ', TRIM(file_in_SI)
status = NF90_OPEN(TRIM(file_in_SI),0,fidTin); call erreur(status,.TRUE.,"Read PARENT sea ice file")

ALLOCATE(  siconc_PAR(mx_PAR,my_PAR)  )   
ALLOCATE(  sithic_PAR(mx_PAR,my_PAR)  )   
ALLOCATE(  snthic_PAR(mx_PAR,my_PAR)  )   

status = NF90_INQ_VARID(fidTin,"siconc",siconc_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"ileadfra",siconc_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"ice_cover",siconc_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"sic",siconc_ID)
call erreur(status,.TRUE.,"No variable identified as sea ice concentration (add other possible names in extract_istate_sea_ice.f90)")

status = NF90_INQ_VARID(fidTin,"sithic",sithic_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"icethic_cea",sithic_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"ice_thickness",sithic_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"iicethic",sithic_ID)
call erreur(status,.TRUE.,"No variable identified as sea ice thickness (add other possible names in extract_istate_sea_ice.f90)")

status = NF90_INQ_VARID(fidTin,"snthic",snthic_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"icesnow_cea",snthic_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"snow_thickness",snthic_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"isnothic",snthic_ID)
if (status .ne. 0) status = NF90_INQ_VARID(fidTin,"isnowthi",snthic_ID)
call erreur(status,.TRUE.,"No variable identified as snow thickness over sea ice (add other possible names in extract_istate_sea_ice.f90)")

status = NF90_GET_VAR(fidTin,siconc_ID,siconc_PAR); call erreur(status,.TRUE.,"getvar_siconc_PAR")
status = NF90_GET_VAR(fidTin,sithic_ID,sithic_PAR); call erreur(status,.TRUE.,"getvar_sithic_PAR")
status = NF90_GET_VAR(fidTin,snthic_ID,snthic_PAR); call erreur(status,.TRUE.,"getvar_snthic_PAR")

status = NF90_CLOSE(fidTin); call erreur(status,.TRUE.,"Close PARENT sea ice file") 

!=================================================================================
! 4- Extract CHILD from PARENT
!=================================================================================

ALLOCATE( siconc_CHLD(mx_CHLD,my_CHLD) )
ALLOCATE( sithic_CHLD(mx_CHLD,my_CHLD) )
ALLOCATE( snthic_CHLD(mx_CHLD,my_CHLD) )

!- Read global attributes of coordinate file to get grid correspondance :
!       i_EXT = ai * i_PAR + bi
!       j_EXT = aj * j_PAR + bj
status = NF90_OPEN(TRIM(file_in_coord_CHLD),0,fidCOORD); call erreur(status,.TRUE.,"read coord input")
status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "ai", ai); call erreur(status,.TRUE.,"read att1")
status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "bi", bi); call erreur(status,.TRUE.,"read att2")
status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "aj", aj); call erreur(status,.TRUE.,"read att3")
status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "bj", bj); call erreur(status,.TRUE.,"read att4")
status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "imin_extraction", imin_EXT); call erreur(status,.TRUE.,"read att5")
status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "jmin_extraction", jmin_EXT); call erreur(status,.TRUE.,"read att6")
status = NF90_CLOSE(fidCOORD)                         ; call erreur(status,.TRUE.,"end read fidCOORD")

! Just extract where ocean points on both grids :
ALLOCATE( missing(mx_CHLD,my_CHLD) )
ALLOCATE( tmp_missing(mx_CHLD,my_CHLD) )
ALLOCATE( tmp_siconc_CHLD(mx_CHLD,my_CHLD) )
ALLOCATE( tmp_sithic_CHLD(mx_CHLD,my_CHLD) )
ALLOCATE( tmp_snthic_CHLD(mx_CHLD,my_CHLD) )
missing(:,:)=0
siconc_CHLD(:,:)=0.d0
sithic_CHLD(:,:)=0.d0
snthic_CHLD(:,:)=0.d0
do iCHLD=1,mx_CHLD
do jCHLD=1,my_CHLD
   iPAR=NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai)
   jPAR=NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
   if ( iPAR .ge. 1 .and. jPAR .ge. 1 ) then
       if ( tmask_PAR(iPAR,jPAR) .eq. 1 ) then
         siconc_CHLD(iCHLD,jCHLD) = siconc_PAR(iPAR,jPAR) * tmask_CHLD(iCHLD,jCHLD) 
         sithic_CHLD(iCHLD,jCHLD) = sithic_PAR(iPAR,jPAR) * tmask_CHLD(iCHLD,jCHLD) 
         snthic_CHLD(iCHLD,jCHLD) = snthic_PAR(iPAR,jPAR) * tmask_CHLD(iCHLD,jCHLD) 
       elseif ( tmask_CHLD(iCHLD,jCHLD) .eq. 1 ) then ! unmasked CHLD but masked PAR
         missing(iCHLD,jCHLD) = 1
       endif
   else ! part of the regional domain not covered by the global domain
       if ( tmask_CHLD(iCHLD,jCHLD) .eq. 1 ) missing(iCHLD,jCHLD) = 1
   endif
enddo
enddo

! Look for closest neighbours where we have missing values:
do kiter=1,nn_iter
  ntest = NINT(sum(sum(FLOAT(missing),2),1))
  write(*,*) '  kiter = ', kiter
  write(*,*) '     nb of pts with missing value: ', ntest
  if ( ntest .eq. 0 ) exit
  tmp_siconc_CHLD(:,:)=siconc_CHLD(:,:)
  tmp_sithic_CHLD(:,:)=sithic_CHLD(:,:)
  tmp_snthic_CHLD(:,:)=snthic_CHLD(:,:)
  tmp_missing(:,:)=missing(:,:)
  do iCHLD=1,mx_CHLD
  do jCHLD=1,my_CHLD
    if ( missing(iCHLD,jCHLD) .eq. 1 ) then
      iout=.FALSE.
      do rs=1,nn_rsmax,1
          iii=iCHLD               ; jjj=jCHLD
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          iii=MIN(iCHLD+rs,mx_CHLD); jjj=jCHLD
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          iii=MAX(iCHLD-rs,1)     ; jjj=jCHLD
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          iii=iCHLD               ; jjj=MIN(jCHLD+rs,my_CHLD)
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          iii=iCHLD               ; jjj=MAX(jCHLD-rs,1)
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          iii=MIN(iCHLD+rs,mx_CHLD); jjj=MIN(jCHLD+rs,my_CHLD)
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          iii=MIN(iCHLD+rs,mx_CHLD); jjj=MAX(jCHLD-rs,1)
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          iii=MAX(iCHLD-rs,1)     ; jjj=MIN(jCHLD+rs,my_CHLD)
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          iii=MAX(iCHLD-rs,1)     ; jjj=MAX(jCHLD-rs,1)
          if ( tmask_CHLD(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then
            iout=.TRUE.
            exit
          endif
      enddo !- rs
      if (iout) then
        tmp_missing(iCHLD,jCHLD) = 0
        tmp_siconc_CHLD(iCHLD,jCHLD) = siconc_CHLD(iii,jjj)
        tmp_sithic_CHLD(iCHLD,jCHLD) = sithic_CHLD(iii,jjj)
        tmp_snthic_CHLD(iCHLD,jCHLD) = snthic_CHLD(iii,jjj)
        !exit
      elseif ( kiter .eq. nn_iter ) then
        write(*,953) iCHLD, jCHLD
        953 FORMAT(' >>> WARNING for point (',2I5,') --> filled with zero sea ice concentration and thickness')
        tmp_missing(iCHLD,jCHLD) = 0
        tmp_siconc_CHLD(iCHLD,jCHLD) = 0.0
        tmp_sithic_CHLD(iCHLD,jCHLD) = 0.0
        tmp_snthic_CHLD(iCHLD,jCHLD) = 0.0
        !exit
      endif
    endif !-if ( missing(iCHLD,jCHLD) .eq. 1 )
  enddo !- jCHLD
  enddo !- iCHLD
  missing(:,:)=tmp_missing(:,:)
  siconc_CHLD(:,:)=tmp_siconc_CHLD(:,:)
  sithic_CHLD(:,:)=tmp_sithic_CHLD(:,:)
  snthic_CHLD(:,:)=tmp_snthic_CHLD(:,:)
enddo !- kiter

!- Smoothing :
if ( nn_smooth .gt. 1 ) then
  write(*,*) 'Smoothing window width = ', nn_smooth
  dij=INT(nn_smooth*0.5)
  tmp_siconc_CHLD(:,:)=siconc_CHLD(:,:)
  tmp_sithic_CHLD(:,:)=sithic_CHLD(:,:)
  tmp_snthic_CHLD(:,:)=snthic_CHLD(:,:)
  do iCHLD=1,mx_CHLD
  do jCHLD=1,my_CHLD
    im1=MAX(iCHLD-dij,1) ; ip1=MIN(iCHLD+dij,mx_CHLD) 
    jm1=MAX(jCHLD-dij,1) ; jp1=MIN(jCHLD+dij,my_CHLD)
    if ( tmask_CHLD(iCHLD,jCHLD) .eq. 1 ) then 
      tmp_siconc_CHLD(iCHLD,jCHLD) =   SUM( SUM( siconc_CHLD(im1:ip1,jm1:jp1) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
      &                           / SUM( SUM(                        1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
      tmp_sithic_CHLD(iCHLD,jCHLD) =   SUM( SUM( sithic_CHLD(im1:ip1,jm1:jp1) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
      &                           / SUM( SUM(                        1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
      tmp_snthic_CHLD(iCHLD,jCHLD) =   SUM( SUM( snthic_CHLD(im1:ip1,jm1:jp1) * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1) &
      &                           / SUM( SUM(                        1.0  * tmask_CHLD(im1:ip1,jm1:jp1), 2), 1)
    else
      tmp_siconc_CHLD(iCHLD,jCHLD) = 0.d0
      tmp_sithic_CHLD(iCHLD,jCHLD) = 0.d0
      tmp_snthic_CHLD(iCHLD,jCHLD) = 0.d0
    endif
  enddo
  enddo
  siconc_CHLD(:,:)=tmp_siconc_CHLD(:,:)
  sithic_CHLD(:,:)=tmp_sithic_CHLD(:,:)
  snthic_CHLD(:,:)=tmp_snthic_CHLD(:,:)
else
  write(*,*) 'No Smoothing'
endif

!--  
DEALLOCATE( tmp_siconc_CHLD, tmp_sithic_CHLD, missing )

!=================================================================================
! 5- Writing initial state for sea ice
!=================================================================================

write(*,*) 'Writing ', TRIM(file_out_SI)

status = NF90_CREATE(TRIM(file_out_SI),NF90_NOCLOBBER,fidTEMP) ; call erreur(status,.TRUE.,'create output temp')

status = NF90_DEF_DIM(fidTEMP,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
status = NF90_DEF_DIM(fidTEMP,"x",mx_CHLD,dimID_x)                               ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidTEMP,"y",my_CHLD,dimID_y)                               ; call erreur(status,.TRUE.,"def_dimID_y")

status = NF90_DEF_VAR(fidTEMP,"time_counter",NF90_DOUBLE,(/dimID_time_counter/),time_counter_ID)
call erreur(status,.TRUE.,"def_var_time_counter_ID")
status = NF90_DEF_VAR(fidTEMP,"siconc",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),siconc_ID)
call erreur(status,.TRUE.,"def_var_siconc_ID")
status = NF90_DEF_VAR(fidTEMP,"sithic",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),sithic_ID)
call erreur(status,.TRUE.,"def_var_sithic_ID")
status = NF90_DEF_VAR(fidTEMP,"snthic",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),snthic_ID)
call erreur(status,.TRUE.,"def_var_snthic_ID")

status = NF90_PUT_ATT(fidTEMP,siconc_ID,"associate","time_counter, y, x")    ; call erreur(status,.TRUE.,"put_att_siconc_ID")
status = NF90_PUT_ATT(fidTEMP,siconc_ID,"long_name","sea ice concentration") ; call erreur(status,.TRUE.,"put_att_siconc_ID")
status = NF90_PUT_ATT(fidTEMP,siconc_ID,"missing_value",0.)                  ; call erreur(status,.TRUE.,"put_att_siconc_ID")
status = NF90_PUT_ATT(fidTEMP,siconc_ID,"_FillValue",0.)                     ; call erreur(status,.TRUE.,"put_att_siconc_ID")
status = NF90_PUT_ATT(fidTEMP,siconc_ID,"units","-")                         ; call erreur(status,.TRUE.,"put_att_siconc_ID")

status = NF90_PUT_ATT(fidTEMP,sithic_ID,"associate","time_counter, y, x")    ; call erreur(status,.TRUE.,"put_att_sithic_ID")
status = NF90_PUT_ATT(fidTEMP,sithic_ID,"long_name","sea ice thickness")     ; call erreur(status,.TRUE.,"put_att_sithic_ID")
status = NF90_PUT_ATT(fidTEMP,sithic_ID,"missing_value",0.)                  ; call erreur(status,.TRUE.,"put_att_sithic_ID")
status = NF90_PUT_ATT(fidTEMP,sithic_ID,"_FillValue",0.)                     ; call erreur(status,.TRUE.,"put_att_sithic_ID")
status = NF90_PUT_ATT(fidTEMP,sithic_ID,"units","m")                         ; call erreur(status,.TRUE.,"put_att_sithic_ID")

status = NF90_PUT_ATT(fidTEMP,snthic_ID,"associate","time_counter, y, x")    ; call erreur(status,.TRUE.,"put_att_snthic_ID")
status = NF90_PUT_ATT(fidTEMP,snthic_ID,"long_name","snow thickness on sea ice") ; call erreur(status,.TRUE.,"put_att_snthic_ID")
status = NF90_PUT_ATT(fidTEMP,snthic_ID,"missing_value",0.)                  ; call erreur(status,.TRUE.,"put_att_snthic_ID")
status = NF90_PUT_ATT(fidTEMP,snthic_ID,"_FillValue",0.)                     ; call erreur(status,.TRUE.,"put_att_snthic_ID")
status = NF90_PUT_ATT(fidTEMP,snthic_ID,"units","m")                         ; call erreur(status,.TRUE.,"put_att_snthic_ID")

status = NF90_PUT_ATT(fidTEMP,time_counter_ID,"title","Time")                  ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidTEMP,time_counter_ID,"long_name","Time axis")         ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidTEMP,time_counter_ID,"standard_name","time")          ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidTEMP,time_counter_ID,"axis","T")                      ; call erreur(status,.TRUE.,"put_att_time_counter_ID")

status = NF90_PUT_ATT(fidTEMP,NF90_GLOBAL,"history","Created using extract_istate_sea_ice.f90")
status = NF90_PUT_ATT(fidTEMP,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL")

status = NF90_ENDDEF(fidTEMP) ; call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidTEMP,time_counter_ID,1.0)   ; call erreur(status,.TRUE.,"var_time_counter_ID")
status = NF90_PUT_VAR(fidTEMP,siconc_ID,siconc_CHLD)  ; call erreur(status,.TRUE.,"var_siconc_ID")
status = NF90_PUT_VAR(fidTEMP,sithic_ID,sithic_CHLD)  ; call erreur(status,.TRUE.,"var_sithic_ID")
status = NF90_PUT_VAR(fidTEMP,snthic_ID,snthic_CHLD)  ; call erreur(status,.TRUE.,"var_snthic_ID")

status = NF90_CLOSE(fidTEMP) ; call erreur(status,.TRUE.,"final")         


end program modif



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
WRITE(*,*) 'MESSAGE: ', TRIM(chaine)
WRITE(*,*) 'ERROR: ', iret
message=NF90_STRERROR(iret)
WRITE(*,*) 'WHICH MEANS:',TRIM(message)
IF ( lstop ) STOP
ENDIF
!
END SUBROUTINE erreur
