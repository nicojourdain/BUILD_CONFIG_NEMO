program modif                                         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, IGE-CNRS, Jan. 2017
!
! Script to extract the bathymetry from a lon/lat dataset (e.g. RTOPO2).
!
! The bathymetry along the boundaries (over a NINT(ai*1.5)-pts halo) is the same as in the
! dataset used as lateral boundary conditions (referred to as "CRS").
!
! 0- Initializations 
! 1- Read RTopo bathymetry and ice shelf draft
! 2- Read grid correspondance with GLO (i.e. extraction coordinates)
! 3- Read coarse bathymetry ("CRS") used for consistent bathymetry along boundaries
! 4- Calculate bathy/isf draft on the REG grid
! 5- Writing new REG bathymetry file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /griddata/ inputdir, file_in_coord_extract, file_in_bathy_extract, file_bathy_out, ln_coarse_bdy, file_in_bathy_bdy, ln_isfcav, &
& nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract, rn_latref, rn_lonref
namelist /bathy_special/ file_spe_bathy, file_spe_isf_draft, ln_dateline, nn_perio
INTEGER                               :: nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract, nn_perio
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: inputdir, file_in_bathy_extract, file_in_coord_extract, file_in_bathy_bdy, config_dir, &
&                                        file_spe_bathy, file_spe_isf_draft
LOGICAL                               :: ln_coarse_bdy, ln_isfcav, ln_dateline
REAL(KIND=4)                          :: rn_latref, rn_lonref

!-- RTOPO variables :
INTEGER :: fidRTOPO1, fidRTOPO2, dimID_latdim, dimID_londim, my_RTOPO, mx_RTOPO, bedrock_topography_ID, lat_ID, lon_ID, &
&          bathy_RTOPO_ID, isf_draft_RTOPO_ID

REAL*4,ALLOCATABLE,DIMENSION(:) :: lat_RTOPO, lon_RTOPO, zlon_RTOPO, tmp1, tmp2, tmp3, tmp4

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: bathy_RTOPO, isf_draft_RTOPO


!-- local variables :
INTEGER :: fidGLO, fidCRS, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, isf_draft_ID, Bathymetry_isf_ID, Bathymetry_ID, &
&          my_GLO, mx_GLO,  my_CRS, mx_CRS,  my_REG, mx_REG, imin_GLO, imax_GLO, jmin_GLO, jmax_GLO, npts, jtmp, Nbox,             &
&          fidCOORDreg, fidCOORDpar, minlon, maxlon, minlat, maxlat, imin_RTOPO, imax_RTOPO, jmin_RTOPO, jmax_RTOPO, iREG, jREG,   &
&          iRTOPO, jRTOPO, iREGm1, iREGp1, jREGm1, jREGp1, kk, mx_tmp, my_tmp, i0, j0, rs, ai, aj, bi, bj

CHARACTER(LEN=150) :: file_bathy_out, file_in_coord_REG

INTEGER, ALLOCATABLE,DIMENSION(:,:) :: nn

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: gphit_REG, glamt_REG, zglamt_REG, isf_draft_REG, Bathymetry_isf_REG, Bathymetry_REG, &
&                                          Bathymetry_CRS, Bathymetry_isf_CRS, isf_draft_CRS, Bathymetry_0, Bathymetry_isf_0, isf_draft_0

INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: imin, imax, jmin, jmax

!-- Regional initial state
REAL*8                                 :: eps
 
!=================================================================================
! 0- Initializations 
!=================================================================================

! Default values (replaced with namelist values if specified):
config_dir        = '.'
file_in_bathy_bdy = 'not_used'

!- read namelist values
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=griddata)
READ (UNIT=1, NML=bathy_special)
CLOSE(1)

! name of regional bathymetry file (output file) :
write(file_bathy_out,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/bathy_meter_',a,'.nc')

! name of regional coordinates file (output file) :
write(file_in_coord_REG,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

eps = 1.d-9

!=================================================================================
! 1- Read lon/lat bathymetry and ice shelf draft
!=================================================================================

write(*,*) 'Reading bathymetry lon/lat dataset in ', TRIM(file_spe_bathy)

status = NF90_OPEN(TRIM(file_spe_bathy),0,fidRTOPO1); call erreur(status,.TRUE.,"Sart read RTopo") 

status = NF90_INQ_DIMID(fidRTOPO1,"latdim",dimID_latdim); call erreur(status,.TRUE.,"inq_dimID_latdim")
status = NF90_INQ_DIMID(fidRTOPO1,"londim",dimID_londim); call erreur(status,.TRUE.,"inq_dimID_londim")

status = NF90_INQUIRE_DIMENSION(fidRTOPO1,dimID_latdim,len=my_RTOPO); call erreur(status,.TRUE.,"inq_dim_latdim")
status = NF90_INQUIRE_DIMENSION(fidRTOPO1,dimID_londim,len=mx_RTOPO); call erreur(status,.TRUE.,"inq_dim_londim")

ALLOCATE(  bathy_RTOPO     (mx_RTOPO,my_RTOPO)  )
ALLOCATE(  isf_draft_RTOPO (mx_RTOPO,my_RTOPO)  ) 
ALLOCATE(  lat_RTOPO (my_RTOPO)  ) 
ALLOCATE(  lon_RTOPO (mx_RTOPO)  ) 
ALLOCATE( zlon_RTOPO (mx_RTOPO)  ) 

status = NF90_INQ_VARID(fidRTOPO1,"bedrock_topography",bathy_RTOPO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRTOPO1,"bathymetry",bathy_RTOPO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRTOPO1,"Bathymetry",bathy_RTOPO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRTOPO1,"BATHYMETRY",bathy_RTOPO_ID)
call erreur(status,.TRUE.,"inq_bathy_RTOPO_ID")
status = NF90_INQ_VARID(fidRTOPO1,"lat",lat_ID);                        call erreur(status,.TRUE.,"inq_lat_ID")
status = NF90_INQ_VARID(fidRTOPO1,"lon",lon_ID);                        call erreur(status,.TRUE.,"inq_lon_ID")

status = NF90_GET_VAR(fidRTOPO1,bathy_RTOPO_ID,bathy_RTOPO); call erreur(status,.TRUE.,"getvar_bathy_RTOPO")
status = NF90_GET_VAR(fidRTOPO1,lat_ID,lat_RTOPO);           call erreur(status,.TRUE.,"getvar_lat")
status = NF90_GET_VAR(fidRTOPO1,lon_ID,lon_RTOPO);           call erreur(status,.TRUE.,"getvar_lon")

status = NF90_CLOSE(fidRTOPO1); call erreur(status,.TRUE.,"End read bathy RTopo")     

zlon_RTOPO(:) = lon_RTOPO(:)
if ( ln_dateline ) then
  where( lon_RTOPO(:) .lt. 0.e0 )
    zlon_RTOPO(:) = 360.e0 + lon_RTOPO(:)
  endwhere
else
  where( lon_RTOPO(:) .gt. 180.e0 )
    zlon_RTOPO(:) = lon_RTOPO(:) - 360.e0
  endwhere
endif

!-----

if ( ln_isfcav ) then

  write(*,*) 'Reading RTopo ice shelf draft in ', TRIM(file_spe_isf_draft)

  status = NF90_OPEN(TRIM(file_spe_isf_draft),0,fidRTOPO2); call erreur(status,.TRUE.,"read isf_draft RTopo") 
  status = NF90_INQ_VARID(fidRTOPO2,"ice_base_topography",isf_draft_RTOPO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRTOPO2,"ice_draft",isf_draft_RTOPO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRTOPO2,"ICE_DRAFT",isf_draft_RTOPO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidRTOPO2,"isf_draft",isf_draft_RTOPO_ID)
  call erreur(status,.TRUE.,"inq_isf_draft_RTOPO_ID")
  status = NF90_GET_VAR(fidRTOPO2,isf_draft_RTOPO_ID,isf_draft_RTOPO); call erreur(status,.TRUE.,"getvar_isf_draft_RTOPO")
  status = NF90_CLOSE(fidRTOPO2); call erreur(status,.TRUE.,"End read isf_draft RTopo")     

endif

!=================================================================================
! 2- Read grid correspondance with GLO (i.e. extraction coordinates)
!=================================================================================

write(*,*) 'Reading REGIONAL lon,lat in ', TRIM(file_in_coord_REG)

status = NF90_OPEN(TRIM(file_in_coord_REG),0,fidCOORDreg); call erreur(status,.TRUE.,"read coord input")

status = NF90_INQ_DIMID(fidCOORDreg,"x",dimID_x)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidCOORDreg,"X",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidCOORDreg,"y",dimID_y)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidCOORDreg,"Y",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")
                                                     
status = NF90_INQUIRE_DIMENSION(fidCOORDreg,dimID_y,len=my_REG); call erreur(status,.TRUE.,"inq_dim_y_REG")
status = NF90_INQUIRE_DIMENSION(fidCOORDreg,dimID_x,len=mx_REG); call erreur(status,.TRUE.,"inq_dim_x_REG")

ALLOCATE(  gphit_REG (mx_REG,my_REG)  )
ALLOCATE(  glamt_REG (mx_REG,my_REG)  )
ALLOCATE( zglamt_REG (mx_REG,my_REG)  )

status = NF90_INQ_VARID(fidCOORDreg,"gphit",nav_lat_ID); call erreur(status,.TRUE.,"inq_gphit_REG_ID")
status = NF90_INQ_VARID(fidCOORDreg,"glamt",nav_lon_ID); call erreur(status,.TRUE.,"inq_glamt_REG_ID")

status = NF90_GET_VAR(fidCOORDreg,nav_lat_ID,gphit_REG); call erreur(status,.TRUE.,"getvar_gphit_REG")
status = NF90_GET_VAR(fidCOORDreg,nav_lon_ID,glamt_REG); call erreur(status,.TRUE.,"getvar_glamt_REG")

status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "ai", ai); call erreur(status,.TRUE.,"read att1")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "bi", bi); call erreur(status,.TRUE.,"read att2")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "aj", aj); call erreur(status,.TRUE.,"read att3")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "bj", bj); call erreur(status,.TRUE.,"read att4")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "imin_extraction", imin_GLO); call erreur(status,.TRUE.,"read att5")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "jmin_extraction", jmin_GLO); call erreur(status,.TRUE.,"read att6")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "imax_extraction", imax_GLO); call erreur(status,.TRUE.,"read att7")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "jmax_extraction", jmax_GLO); call erreur(status,.TRUE.,"read att8")

status = NF90_CLOSE(fidCOORDreg)                       ; call erreur(status,.TRUE.,"end read fidCOORDreg")

zglamt_REG(:,:)=glamt_REG(:,:)
if ( ln_dateline ) then
  where( glamt_REG(:,:) .lt. 0.e0 )
    zglamt_REG(:,:)=glamt_REG(:,:)+360.e0
  endwhere
else
  where( glamt_REG(:,:) .gt. 180.e0 )
    zglamt_REG(:,:)=glamt_REG(:,:)-360.e0
  endwhere
endif

!- lat/lon REG boundaries to limit loops on RTOPO's lon/lat
minlon = MINVAL(zglamt_REG-0.5)
maxlon = MAXVAL(zglamt_REG+0.5)
minlat = MINVAL( gphit_REG-0.5)
maxlat = MAXVAL( gphit_REG+0.5)

ALLOCATE( tmp1(mx_RTOPO), tmp2(mx_RTOPO) )
ALLOCATE( tmp3(my_RTOPO), tmp4(my_RTOPO) )
do iRTOPO=1,mx_RTOPO
  tmp1(iRTOPO) = abs( zlon_RTOPO(iRTOPO) - minlon )
  tmp2(iRTOPO) = abs( zlon_RTOPO(iRTOPO) - maxlon )
enddo
do jRTOPO=1,my_RTOPO
  tmp3(jRTOPO) = abs(  lat_RTOPO(jRTOPO) - minlat )
  tmp4(jRTOPO) = abs(  lat_RTOPO(jRTOPO) - maxlat )
enddo
imin_RTOPO = MINLOC( tmp1(:), 1 ) - 50
imax_RTOPO = MINLOC( tmp2(:), 1 ) + 50
jmin_RTOPO = MINLOC( tmp3(:), 1 ) - 50
jmax_RTOPO = MINLOC( tmp4(:), 1 ) + 50
DEALLOCATE( lon_RTOPO, tmp1, tmp2, tmp3, tmp4 )

write(*,527) imin_RTOPO, imax_RTOPO, jmin_RTOPO, jmax_RTOPO
527 FORMAT('Restricting search on RTOPO to i=',i5,':',i5,' and j=',i5,':',i5)

!=================================================================================
! 3- Read coarse bathymetry used for consistent bathymetry along boundaries
!=================================================================================

write(*,*) 'Reading coarse bathymetry for consistent boundaries: ', TRIM(file_in_bathy_bdy)

status = NF90_OPEN(TRIM(file_in_bathy_bdy),0,fidCRS); call erreur(status,.TRUE.,"read_coarse_bathymetry") 

status = NF90_INQ_DIMID(fidCRS,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_CRS")
status = NF90_INQ_DIMID(fidCRS,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_CRS")

status = NF90_INQUIRE_DIMENSION(fidCRS,dimID_y,len=my_CRS); call erreur(status,.TRUE.,"inq_dim_y_CRS")
status = NF90_INQUIRE_DIMENSION(fidCRS,dimID_x,len=mx_CRS); call erreur(status,.TRUE.,"inq_dim_x_CRS")

ALLOCATE(  Bathymetry_CRS (mx_CRS,my_CRS)  ) 
if ( ln_isfcav ) then
  ALLOCATE(  Bathymetry_isf_CRS (mx_CRS,my_CRS)  )
  ALLOCATE(  isf_draft_CRS (mx_CRS,my_CRS)  )
endif

status = NF90_INQ_VARID(fidCRS,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"latitude",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID_CRS")
status = NF90_INQ_VARID(fidCRS,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"longitude",nav_lon_ID)
call erreur(status,.TRUE.,"inq_nav_lon_ID_CRS")

!NB: here we try to read two bathymetry variables: one without ice shelf cavities ("Bathymetry")
!    and, if ln_isfcav is true,  one with ice shelf cavities ("Bathymetry_isf").
!    If no "Bathymetry_isf" is found, we use a single bathymetry variable.
status = NF90_INQ_VARID(fidCRS,"Bathymetry",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"bathymetry",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"bathy_metry",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"bathy",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"Bathy",Bathymetry_ID)
call erreur(status,.TRUE.,"inq_Bathymetry_ID_CRS")
status = NF90_GET_VAR(fidCRS,Bathymetry_ID,Bathymetry_CRS); call erreur(status,.TRUE.,"getvar_Bathymetry_CRS")
if ( ln_isfcav ) then
  status = NF90_INQ_VARID(fidCRS,"Bathymetry_isf",Bathymetry_isf_ID)
  if ( status .eq. 0 ) then
    status = NF90_GET_VAR(fidCRS,Bathymetry_isf_ID,Bathymetry_isf_CRS); call erreur(status,.TRUE.,"getvar_Bathymetry_isf_CRS")
  else
    Bathymetry_isf_CRS(:,:) = Bathymetry_CRS(:,:)
  endif
  status = NF90_INQ_VARID(fidCRS,"isf_draft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"isfdraft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"draft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"ice_shelf_draft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"Ice_shelf_draft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidCRS,"ice_shelf_base",isf_draft_ID)
  if ( status .eq. 0 ) then
    status = NF90_GET_VAR(fidCRS,isf_draft_ID,isf_draft_CRS); call erreur(status,.TRUE.,"getvar_isf_draft_CRS")
  else
    isf_draft_CRS(:,:) = 0.e0
     write(*,*) 'WARNING: no ice shelf draft found. Allowed names are :'
     write(*,*) '         "isf_draft", "isfdraft", "draft", "ice_shelf_draft", "Ice_shelf_draft", "ice_shelf_base"'
     write(*,*) '         >>>>> ASSUMING ice_draft = 0 for the COARSE data'
  endif
endif

status = NF90_CLOSE(fidCRS); call erreur(status,.TRUE.,"close_coarse_bathy_file")

!===================================================================
! put COARSE grid in a npts-pts halo (+ transition in another halo)

npts=CEILING(ai*1.5)
write(*,*) 'put COARSE grid bathymetry in a npts-pts halo, with npts = ', npts

!=================================================================================
! 4- Calculate bathy/isf draft on the REG grid :
!=================================================================================

write(*,*) 'Calculating bathymetry on the regional grid...'
    
ALLOCATE(  nn                 (mx_REG,my_REG)  )
ALLOCATE(  isf_draft_REG      (mx_REG,my_REG)  )
ALLOCATE(  Bathymetry_isf_REG (mx_REG,my_REG)  )
ALLOCATE(  Bathymetry_REG     (mx_REG,my_REG)  )
ALLOCATE(  isf_draft_0        (mx_REG,my_REG)  )
ALLOCATE(  Bathymetry_isf_0   (mx_REG,my_REG)  )
ALLOCATE(  Bathymetry_0       (mx_REG,my_REG)  ) 
 
nn(:,:) = 0 
Bathymetry_isf_REG (:,:) = 0.e0
isf_draft_REG      (:,:) = 0.e0
Bathymetry_REG     (:,:) = 0.e0

do iREG=1,mx_REG
write(*,*) 'iREG = ', iREG
do jREG=1,my_REG

  if ( nn_perio .eq. 1 ) then !- periodic
    if    ( iREG+1 .gt. mx_REG ) then
      iREGp1 = 3
      iREGm1 = iREG - 1
    elseif ( iREG-1 .lt. 1     ) then
      iREGp1 = iREG + 1
      iREGm1 = mx_REG-2
    else
      iREGp1 = iREG + 1
      iREGm1 = iREG - 1
    endif
  elseif ( nn_perio .eq. 0 ) then
    if    ( iREG+1 .gt. mx_REG ) then
      iREGp1 = iREG
      iREGm1 = iREG - 1
    elseif ( iREG-1 .lt. 1     ) then
      iREGp1 = iREG + 1
      iREGm1 = iREG
    else
      iREGp1 = iREG + 1
      iREGm1 = iREG - 1
    endif
  else
    write(*,*) '~!@#$%^* nn_perio must be either 0 or 1 >>>>> stop !!'
    stop
  endif
  !-
  if    ( jREG+1 .gt. my_REG ) then
    jREGp1 = jREG
    jREGm1 = jREG - 1
  elseif ( iREG-1 .lt. 1     ) then
    jREGp1 = jREG + 1
    jREGm1 = jREG
  else
    jREGp1 = jREG + 1
    jREGm1 = jREG - 1
  endif
  !-
  iREGm1 = MAX( MIN( iREGm1, mx_REG ), 1 )
  iREGp1 = MAX( MIN( iREGp1, mx_REG ), 1 )
  jREGm1 = MAX( MIN( jREGm1, my_REG ), 1 )
  jREGp1 = MAX( MIN( jREGp1, my_REG ), 1 )
  !-
  
  do iRTOPO=imin_RTOPO,imax_RTOPO
  do jRTOPO=jmin_RTOPO,jmax_RTOPO
 
    !-NB: we don't care too much about i=1 and i=mx (same for j) because it is masked anyway...    
    if (       zlon_RTOPO(iRTOPO) .ge. zglamt_REG(iREG,jREG) - 0.5*(zglamt_REG(iREG,jREG)-zglamt_REG(iREGm1,jREG)) &
    &    .and. zlon_RTOPO(iRTOPO) .lt. zglamt_REG(iREG,jREG) + 0.5*(zglamt_REG(iREGp1,jREG)-zglamt_REG(iREG,jREG)) &
    &    .and.  lat_RTOPO(jRTOPO) .ge.  gphit_REG(iREG,jREG) - 0.5*( gphit_REG(iREG,jREG)- gphit_REG(iREG,jREGm1)) &
    &    .and.  lat_RTOPO(jRTOPO) .lt.  gphit_REG(iREG,jREG) + 0.5*( gphit_REG(iREG,jREGp1)- gphit_REG(iREG,jREG)) ) then

      !Bathymetry_isf_REG (iREG,jREG) = Bathymetry_isf_REG (iREG,jREG) - MIN(0.0, bathy_RTOPO     (iRTOPO,jRTOPO))
      !isf_draft_REG      (iREG,jREG) = isf_draft_REG      (iREG,jREG) - MIN(0.0, isf_draft_RTOPO (iRTOPO,jRTOPO))
      Bathymetry_isf_REG (iREG,jREG) = Bathymetry_isf_REG (iREG,jREG) - bathy_RTOPO     (iRTOPO,jRTOPO)
      isf_draft_REG      (iREG,jREG) = isf_draft_REG      (iREG,jREG) - isf_draft_RTOPO (iRTOPO,jRTOPO)

      nn(iREG,jREG) = nn(iREG,jREG) + 1

    endif 

  enddo
  enddo

enddo
enddo

!-----

where ( nn(:,:) .ge. 1 )
  Bathymetry_isf_REG (:,:) = Bathymetry_isf_REG (:,:) / nn(:,:)
  isf_draft_REG      (:,:) = isf_draft_REG      (:,:) / nn(:,:)
elsewhere
  Bathymetry_isf_REG (:,:) = 99999999.9
  isf_draft_REG      (:,:) = 99999999.9
endwhere

!-----
! Interpolation at points with nn=0 :
! Only look for first neighbour
! needs to be adapted if many missing points

do iREG=1,mx_REG
do jREG=1,my_REG

  if ( nn(iREG,jREG) .eq. 0 ) then

    do rs=1,10
      if ( nn_perio .eq. 1 ) then !- periodic
        if    ( iREG+rs .gt. mx_REG ) then
          iREGp1 = 2+rs
        else
          iREGp1 = iREG + rs
        endif
      elseif ( nn_perio .eq. 0 ) then
        if    ( iREG+rs .gt. mx_REG ) then
          iREGp1 = mx_REG
          exit ! with nn(iREGp1,jREG)=0
        else
          iREGp1 = iREG + rs
        endif
      endif
      if ( nn(iREGp1,jREG) .ne. 0 ) exit ! with nn(iREGp1,jREG)>0
    enddo

    do rs=1,10
      if ( nn_perio .eq. 1 ) then !- periodic
        if ( iREG-rs .lt. 1     ) then
          iREGm1 = mx_REG-rs-1
        else
          iREGm1 = iREG - rs
        endif
      elseif ( nn_perio .eq. 0 ) then
        if ( iREG-rs .lt. 1     ) then
          iREGm1 = 1
          exit ! with nn(iREGm1,jREG)=0
        else
          iREGm1 = iREG - rs
        endif
      endif
      if ( nn(iREGm1,jREG) .gt. 0 ) exit ! with nn(iREGm1,jREG)>0
    enddo

    do rs=1,10
      if    ( jREG+rs .gt. my_REG ) then
        jREGp1 = my_REG
        exit ! with nn(iREG,jREGp1)=0
      else
        jREGp1 = jREG + rs
      endif
      if ( nn(iREG,jREGp1) .gt. 0 ) exit ! with nn(iREG,jREGp1)>0
    enddo

    do rs=1,10
      if    ( jREG-rs .lt. 1 ) then
        jREGm1 = 1
        exit ! with nn(iREG,jREGm1)=0
      else
        jREGm1 = jREG - rs
      endif
      if ( nn(iREG,jREGm1) .gt. 0 ) exit ! with nn(iREG,jREGm1)>0
    enddo

    Bathymetry_isf_REG(iREG,jREG) = (   nn(iREGm1,jREG  ) * (iREGp1-iREG) * Bathymetry_isf_REG(iREGm1,jREG  )   &
    &                                 + nn(iREGp1,jREG  ) * (iREG-iREGm1) * Bathymetry_isf_REG(iREGp1,jREG  )   &
    &                                 + nn(iREG  ,jREGm1) * (jREGp1-jREG) * Bathymetry_isf_REG(iREG  ,jREGm1)   &
    &                                 + nn(iREG  ,jREGp1) * (jREG-jREGm1) * Bathymetry_isf_REG(iREG  ,jREGp1) ) &
    &                             / (   nn(iREGm1,jREG  ) * (iREGp1-iREG)                                       &
    &                                 + nn(iREGp1,jREG  ) * (iREG-iREGm1)                                       &
    &                                 + nn(iREG  ,jREGm1) * (jREGp1-jREG)                                       &
    &                                 + nn(iREG  ,jREGp1) * (jREG-jREGm1)                                     )

    isf_draft_REG(iREG,jREG) = (   nn(iREGm1,jREG  ) * (iREGp1-iREG) * isf_draft_REG(iREGm1,jREG  )   &
    &                            + nn(iREGp1,jREG  ) * (iREG-iREGm1) * isf_draft_REG(iREGp1,jREG  )   &
    &                            + nn(iREG  ,jREGm1) * (jREGp1-jREG) * isf_draft_REG(iREG  ,jREGm1)   &
    &                            + nn(iREG  ,jREGp1) * (jREG-jREGm1) * isf_draft_REG(iREG  ,jREGp1) ) &
    &                        / (   nn(iREGm1,jREG  ) * (iREGp1-iREG)                                       &
    &                            + nn(iREGp1,jREG  ) * (iREG-iREGm1)                                       &
    &                            + nn(iREG  ,jREGm1) * (jREGp1-jREG)                                       &
    &                            + nn(iREG  ,jREGp1) * (jREG-jREGm1)                                     )

   endif

enddo
enddo

!---------

write(*,*) '>i=125 ', nn(125,20), Bathymetry_isf_REG(125,20), isf_draft_REG(125,20)
write(*,*) '>i=624 ', nn(624,20), Bathymetry_isf_REG(624,20), isf_draft_REG(624,20)

where( isf_draft_REG(:,:) .lt. 0.e0 )
  isf_draft_REG(:,:) = 0.e0
endwhere

where( Bathymetry_isf_REG(:,:) .lt. 0.e0 )
  Bathymetry_isf_REG(:,:) = 0.e0
endwhere

where( isf_draft_REG(:,:) .gt. 1.e0 .and. isf_draft_REG(:,:) .lt. 1.e4 )
  Bathymetry_REG(:,:) = 0.e0
elsewhere
  Bathymetry_REG(:,:) = Bathymetry_isf_REG(:,:)
endwhere

write(*,*) '#i=125 ', Bathymetry_REG(125,20), Bathymetry_isf_REG(125,20), isf_draft_REG(125,20)
write(*,*) '#i=624 ', Bathymetry_REG(624,20), Bathymetry_isf_REG(624,20), isf_draft_REG(624,20)

!---------------------------------------
! Adjust bathy along the edges of the REG grid :

write(*,*) 'Halo with bathy from the dataset used as lateral boundaries...'

if ( ln_coarse_bdy ) then
    
    !---------------------------------------
    ! Adjust bathy along the edges of the REG grid :
    
    write(*,*) 'Halo with bathy from coarse resolution...'
    
    !=== put exactly the coarse data (CRS) over a npts-point halo :
    do jREG=1,my_REG
      jtmp = NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        !-- Western BDY :
        do iREG=1,npts
          if ( ln_isfcav) then
            isf_draft_REG     (iREG,jREG) = isf_draft_CRS      ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
            Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
          endif
          Bathymetry_REG    (iREG,jREG) = Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
        enddo
        !--- Eastern BDY :
        do iREG=mx_REG-npts+1,mx_REG
          if ( ln_isfcav) then
            isf_draft_REG     (iREG,jREG) = isf_draft_CRS      ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
            Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
          endif
          Bathymetry_REG    (iREG,jREG) = Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    !-- Southern BDY :
    do jREG=1,npts
      jtmp = NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iREG=npts+1,mx_REG-npts
          if ( ln_isfcav) then
            isf_draft_REG     (iREG,jREG) = isf_draft_CRS      ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
            Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
          endif
          Bathymetry_REG    (iREG,jREG) = Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    !--- Northern BDY :
    do jREG=my_REG-npts+1,my_REG
      jtmp = NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iREG=npts+1,mx_REG-npts
          if ( ln_isfcav) then
            isf_draft_REG     (iREG,jREG) = isf_draft_CRS      ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
            Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
          endif
          Bathymetry_REG    (iREG,jREG) = Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    
    write(*,*) 'Smooth transition...'

    Bathymetry_0     = Bathymetry_REG
    Bathymetry_isf_0 = Bathymetry_isf_REG
    isf_draft_0      = isf_draft_REG

    !=== smooth transition from BDY to the interpolated fields (over npts points again) :
    do jREG=npts+1,my_REG-npts
      jtmp = NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj)
      !-- Western BDY :
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iREG=npts+1,2*npts
          if ( ln_isfcav) then
            if ( isf_draft_0(iREG,jREG) .gt. 1.0 .and. isf_draft_CRS(NINT(FLOAT(iREG+imin_GLO-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_REG(iREG,jREG) = (  (2*npts+1-iREG) * isf_draft_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
              &                           + (iREG-npts) * isf_draft_0(iREG,jREG) ) / (npts+1)
            !else
            !  isf_draft_REG   (iREG,jREG) = 0.0
            endif
            if ( Bathymetry_isf_0(iREG,jREG) .gt. 1.0 .and. Bathymetry_isf_CRS(NINT(FLOAT(iREG+imin_GLO-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_REG(iREG,jREG) = (  (2*npts+1-iREG) * Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
              &                                + (iREG-npts) * Bathymetry_isf_0(iREG,jREG) ) / (npts+1)
            !else
            !  Bathymetry_isf_REG   (iREG,jREG) = 0.0
            endif
          endif
          Bathymetry_REG(iREG,jREG) = (   (2*npts+1-iREG) * Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
          &                             + (iREG-npts) * Bathymetry_0(iREG,jREG) ) / (npts+1)
        enddo
        !-- Eastern BDY
        do iREG=mx_REG-2*npts+1,mx_REG-npts
          if ( ln_isfcav) then
            if ( isf_draft_0(iREG,jREG) .gt. 1.0 .and. isf_draft_CRS(NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp) .gt. 1.0 ) then
              isf_draft_REG(iREG,jREG) = (   (iREG-mx_REG+2*npts) * isf_draft_CRS( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
              &                            + (mx_REG-npts+1-iREG) * isf_draft_0(iREG,jREG) ) / (npts+1)
            !else
            !  isf_draft_REG   (iREG,jREG) = 0.0
            endif
            if ( Bathymetry_isf_0(iREG,jREG) .gt. 1.0 .and. Bathymetry_isf_CRS(NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp) .gt. 1.0 ) then
              Bathymetry_isf_REG(iREG,jREG) = (   (iREG-mx_REG+2*npts) * Bathymetry_isf_CRS( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
              &                                 + (mx_REG-npts+1-iREG) * Bathymetry_isf_0(iREG,jREG) ) / (npts+1)
            !else
            !  Bathymetry_isf_REG   (iREG,jREG) = 0.0
            endif
          endif
          Bathymetry_REG(iREG,jREG) = (   (iREG-mx_REG+2*npts) * Bathymetry_CRS( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
          &                             + (mx_REG-npts+1-iREG) * Bathymetry_0(iREG,jREG) ) / (npts+1)
        enddo
      endif
    enddo
    !-- Southern BDY :
    do jREG=npts+1,2*npts
      jtmp = NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iREG=2*npts+1,mx_REG-2*npts
          if ( ln_isfcav) then
            if ( isf_draft_0(iREG,jREG) .gt. 1.0 .and. isf_draft_CRS(NINT(FLOAT(iREG+imin_GLO-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_REG(iREG,jREG) = (   (2*npts+1-jREG) * isf_draft_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
              &                            + (jREG-npts) * isf_draft_0(iREG,jREG) ) / (npts+1)
            !else
            !  isf_draft_REG   (iREG,jREG) = 0.0
            endif
            if ( Bathymetry_isf_0(iREG,jREG) .gt. 1.0 .and. Bathymetry_isf_CRS(NINT(FLOAT(iREG+imin_GLO-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_REG(iREG,jREG) = (   (2*npts+1-jREG) * Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
              &                                 + (jREG-npts) * Bathymetry_isf_0(iREG,jREG) ) / (npts+1)
            !else
            !  Bathymetry_isf_REG   (iREG,jREG) = 0.0
            endif
          endif
          Bathymetry_REG(iREG,jREG) = (   (2*npts+1-jREG) * Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
          &                             + (jREG-npts) * Bathymetry_0(iREG,jREG) ) / (npts+1)
        enddo
      endif
    enddo
    !-- Northern BDY :
    do jREG=my_REG-2*npts+1,my_REG-npts
      jtmp = NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iREG=2*npts+1,mx_REG-2*npts
          if ( ln_isfcav) then
            if ( isf_draft_0(iREG,jREG) .gt. 1.0 .and. isf_draft_CRS(NINT(FLOAT(iREG+imin_GLO-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_REG(iREG,jREG) = (   (jREG-my_REG+2*npts) * isf_draft_CRS( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
              &                            + (my_REG-npts+1-jREG) * isf_draft_0(iREG,jREG) ) / (npts+1)
            !else
            !  isf_draft_REG   (iREG,jREG) = 0.0
            endif
            if ( Bathymetry_isf_0(iREG,jREG) .gt. 1.0 .and. Bathymetry_isf_CRS(NINT(FLOAT(iREG+imin_GLO-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_REG(iREG,jREG) = (   (jREG-my_REG+2*npts) * Bathymetry_isf_CRS( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
              &                                 + (my_REG-npts+1-jREG) * Bathymetry_isf_0(iREG,jREG) ) / (npts+1)
            !else
            !  Bathymetry_isf_REG   (iREG,jREG) = 0.0
            endif
            Bathymetry_isf_REG(iREG,jREG) = (   (jREG-my_REG+2*npts) * Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
            &                                 + (my_REG-npts+1-jREG) * Bathymetry_isf_0(iREG,jREG) ) / (npts+1)
          endif
          Bathymetry_REG(iREG,jREG) = (   (jREG-my_REG+2*npts) * Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), jtmp ) &
          &                             + (my_REG-npts+1-jREG) * Bathymetry_0(iREG,jREG) ) / (npts+1)
        enddo
      endif
    enddo
   
endif
   
!=================================================================================
! 4- Manual corrections for WED12 :
!=================================================================================

if ( TRIM(config) == 'WED12' ) then
  
    write(*,*) 'Special correction for config ', TRIM(config)
  
    ! Note that you DO NOT have to change the following in you change the domain size
    ! through modifications of nn_imin_extract, nn_imax_extract, ... in the namelist 

    ! To keep the boxes at the same position:
    i0 = 2464 - imin_GLO
    j0 =  151 - jmin_GLO

    !! correction to avoid a closed cavity of 2x2x2 pts (identified after first
    !! mesh_mask creation)
    !isf_draft_REG     (i0+241:i0+242,j0+667:j0+668) = 0.0
    !Bathymetry_isf_REG(i0+241:i0+242,j0+667:j0+668) = 0.0
    !
    !! no isf along eastern boundary (adapt manually to adjust more accurately) :
    !isf_draft_REG     (i0+1095:mx_REG,j0+668:j0+703) = 0.0
    !Bathymetry_isf_REG(i0+1095:mx_REG,j0+668:j0+703) = Bathymetry_REG(i0+1095:mx_REG,j0+668:j0+703)

    ! boxes to fill the Bellingshausen Sea : filled over [imin:imax,jmin:my]
    Nbox = 8
    ALLOCATE( imin(Nbox), imax(Nbox), jmin(Nbox), jmax(Nbox) )
    imin = (/   1+i0 , 192+i0 , 213+i0 , 237+i0 , 254+i0 , 275+i0 , 287+i0 , 299+i0 /)
    imax = (/ 191+i0 , 212+i0 , 236+i0 , 253+i0 , 274+i0 , 286+i0 , 298+i0 , &
    &           NINT(FLOAT(325+i0+imin_GLO-1-bi)/ai)*ai+bi-imin_GLO+1-1 /) ! WARNING: last number must match with CRS (i.e. next unmasked point neads to be on CRS) !!
    jmin = (/ 494+j0 , 807+j0 , 835+j0 , 862+j0 , 876+j0 , 894+j0 ,            &
    &           NINT(FLOAT(899+j0+jmin_GLO-1-bj)/aj)*aj+bj-jmin_GLO+1+1 ,&
    &           NINT(FLOAT(903+j0+jmin_GLO-1-bj)/aj)*aj+bj-jmin_GLO+1+1 /) ! WARNING: last two numbers must match with CRS (i.e. next unmasked point neads to be on CRS) !!
    jmax(:) = my_REG

    do kk=1,Nbox
      imin(kk) = MIN( MAX( 1, imin(kk) ), mx_REG )
      imax(kk) = MIN( MAX( 1, imax(kk) ), mx_REG )
      jmin(kk) = MIN( MAX( 1, jmin(kk) ), my_REG )
      jmax(kk) = MIN( MAX( 1, jmax(kk) ), my_REG )
    enddo     
 
    write(*,*) 'Note for future bdy building:'
    write(*,*) '                    '
    write(*,*) '  ii_bdy_west(1)  = ', imax(8)+1
    write(*,*) '  j1_bdy_west(1)  = ', jmin(8)
    write(*,*) '  j2_bdy_west(1)  = ', my_REG-1 ! given that jmax(8)=my_REG-1
    write(*,*) '                    '
    write(*,*) '  i1_bdy_north(1) = ', imin(8)
    write(*,*) '  i2_bdy_north(1) = ', imax(8)
    write(*,*) '  jj_bdy_north(1) = ', jmin(8)-1
    write(*,*) '                    '
    write(*,*) '  i1_bdy_north(2) = ', imax(8)+2
    write(*,*) '  i2_bdy_north(2) = ', mx_REG-2 ! -2 because east bdy already contains mx_REG-1
    write(*,*) '  jj_bdy_north(2) = ', my_REG-1
    write(*,*) '                    '

    !----- put coarse data (CRS) along modified North-Western corner :
    !- bdy_west(1) :
    do jREG=jmin(8),my_REG
      do iREG=imax(8)+1,imax(8)+npts
        isf_draft_REG     (iREG,jREG) = isf_draft_CRS      ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) )
        Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) )
        Bathymetry_REG    (iREG,jREG) = Bathymetry_CRS     ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) )
      enddo
    enddo
    !- bdy_north(1) :
    do jREG=jmin(8)-npts, jmin(8)-1
      do iREG=imin(8),imax(8)
        isf_draft_REG     (iREG,jREG) = isf_draft_CRS      ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) )
        Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) )
        Bathymetry_REG    (iREG,jREG) = Bathymetry_CRS     ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) )
      enddo
    enddo

    !---- smooth transitions :
    !- bdy_west(1) :
    do jREG=jmin(8),my_REG
      do iREG=imax(8)+npts+1,imax(8)+2*npts
        isf_draft_REG     (iREG,jREG) = (   (imax(8)+2*npts-iREG+1) * isf_draft_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                 +     (iREG-imax(8)-npts) * isf_draft_0 ( iREG, jREG ) ) / (npts+1)
        Bathymetry_isf_REG(iREG,jREG) = (   (imax(8)+2*npts-iREG+1) * Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                 +     (iREG-imax(8)-npts) * Bathymetry_isf_0 ( iREG, jREG ) ) / (npts+1)
        Bathymetry_REG    (iREG,jREG) = (   (imax(8)+2*npts-iREG+1) * Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                 +     (iREG-imax(8)-npts) * Bathymetry_0 ( iREG, jREG ) ) / (npts+1)
      enddo
    enddo
    !- bdy_north(1) :
    do jREG=jmin(8)-2*npts,jmin(8)-npts-1
      do iREG=imin(8),imax(8)
        isf_draft_REG     (iREG,jREG) = ( (jREG-jmin(8)+2*npts+1)    * isf_draft_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                      + (jmin(8)-npts-jREG) * isf_draft_0 ( iREG, jREG ) ) / (npts+1)                     
        Bathymetry_isf_REG(iREG,jREG) = ( (jREG-jmin(8)+2*npts+1)    * Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                      + (jmin(8)-npts-jREG) * Bathymetry_isf_0 ( iREG, jREG ) ) / (npts+1)                     
        Bathymetry_REG    (iREG,jREG) = ( (jREG-jmin(8)+2*npts+1)    * Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                      + (jmin(8)-npts-jREG) * Bathymetry_0 ( iREG, jREG ) ) / (npts+1)                     
      enddo
    enddo

    !- fill bellingshausen:
    do kk=1,Nbox
      isf_draft_REG      (imin(kk):imax(kk),jmin(kk):my_REG) = 0.0
      Bathymetry_isf_REG (imin(kk):imax(kk),jmin(kk):my_REG) = 0.0
      Bathymetry_REG     (imin(kk):imax(kk),jmin(kk):my_REG) = 0.0
    enddo

elseif ( TRIM(config) == 'AMUXL12' ) then

    write(*,*) 'Special correction for config ', TRIM(config)
 
    ! Note that you DO NOT have to change the following in you change the domain size
    ! through modifications of nn_imin_extract, nn_imax_extract, ... in the namelist 

    ! To keep the boxes at the same position:
    i0 = 1771 - imin_GLO
    j0 =   30 - jmin_GLO
    write(*,*) '   i0 = ', i0
    write(*,*) '   j0 = ', j0

    ! manual correction to avoid isolated open ocean:
    Nbox = 1 
    ALLOCATE( imin(Nbox), imax(Nbox), jmin(Nbox), jmax(Nbox) )
    imin(1) = i0+311
    imax(1) = i0+312
    jmin(1) = j0+73
    jmax(1) = j0+73    
    !-
    do kk=1,Nbox
      imin(kk) = MIN( MAX( 1, imin(kk) ), mx_REG )
      imax(kk) = MIN( MAX( 1, imax(kk) ), mx_REG )
      jmin(kk) = MIN( MAX( 1, jmin(kk) ), my_REG )
      jmax(kk) = MIN( MAX( 1, jmax(kk) ), my_REG )
    enddo
    !-
    isf_draft_REG(imin(1):imax(1),jmin(1):jmax(1)) = Bathymetry_isf_REG(imin(1):imax(1),jmin(1):jmax(1))

    ! no isf over a safety zone (2*npts wide halo) from the eastern and western BDY :
    isf_draft_REG     (1:2*npts,:) = 0.0
    Bathymetry_isf_REG(1:2*npts,:) = Bathymetry_REG(1:2*npts,:)
    isf_draft_REG     (mx_REG-2*npts+1:mx_REG,:) = 0.0
    Bathymetry_isf_REG(mx_REG-2*npts+1:mx_REG,:) = Bathymetry_REG(mx_REG-2*npts+1:mx_REG,:)

endif

!=================================================================================
! 5- Writing new REG bathymetry file :
!=================================================================================

write(*,*) 'Creating ', TRIM(file_bathy_out)

status = NF90_CREATE(TRIM(file_bathy_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')                     

status = NF90_DEF_DIM(fidM,"y",my_REG,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx_REG,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")

status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lat_ID); call erreur(status,.TRUE.,"def_var_nav_lat_ID")
status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lon_ID); call erreur(status,.TRUE.,"def_var_nav_lon_ID")
if ( ln_isfcav ) then
  status = NF90_DEF_VAR(fidM,"isf_draft",NF90_FLOAT,(/dimID_x,dimID_y/),isf_draft_ID); call erreur(status,.TRUE.,"def_var_isf_draft_ID")
  status = NF90_DEF_VAR(fidM,"Bathymetry_isf",NF90_FLOAT,(/dimID_x,dimID_y/),Bathymetry_isf_ID); call erreur(status,.TRUE.,"def_var_Bathymetry_isf_ID")
endif
status = NF90_DEF_VAR(fidM,"Bathymetry",NF90_FLOAT,(/dimID_x,dimID_y/),Bathymetry_ID); call erreur(status,.TRUE.,"def_var_Bathymetry_ID")

status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")
call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"longname","Latitude")
call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"standard_name","latitude")
call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")
call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"longname","Longitude")
call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"standard_name","longitude")
call erreur(status,.TRUE.,"put_att_nav_lon_ID")
if ( ln_isfcav ) then
  status = NF90_PUT_ATT(fidM,isf_draft_ID,"long_name","ice-shelf draft")
  call erreur(status,.TRUE.,"put_att_isf_draft_ID")
  status = NF90_PUT_ATT(fidM,isf_draft_ID,"units","m")
  call erreur(status,.TRUE.,"put_att_isf_draft_ID")
  status = NF90_PUT_ATT(fidM,isf_draft_ID,"coordinates","nav_lat nav_lon")
  call erreur(status,.TRUE.,"put_att_isf_draft_ID")
  status = NF90_PUT_ATT(fidM,Bathymetry_isf_ID,"long_name","bathymetry")
  call erreur(status,.TRUE.,"put_att_Bathymetry_isf_ID")
  status = NF90_PUT_ATT(fidM,Bathymetry_isf_ID,"units","m")
  call erreur(status,.TRUE.,"put_att_Bathymetry_isf_ID")
  status = NF90_PUT_ATT(fidM,Bathymetry_isf_ID,"coordinates","nav_lat nav_lon")
  call erreur(status,.TRUE.,"put_att_Bathymetry_isf_ID")
endif
status = NF90_PUT_ATT(fidM,Bathymetry_ID,"long_name","bathymetry with masked ice-shelf cavities")
call erreur(status,.TRUE.,"put_att_Bathymetry_ID")
status = NF90_PUT_ATT(fidM,Bathymetry_ID,"units","m")
call erreur(status,.TRUE.,"put_att_Bathymetry_ID")
status = NF90_PUT_ATT(fidM,Bathymetry_ID,"coordinates","nav_lat nav_lon")
call erreur(status,.TRUE.,"put_att_Bathymetry_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bathy_special_lonlat.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL1")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"data",TRIM(file_spe_bathy)) ; call erreur(status,.TRUE.,"put_att_GLOBAL2")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imin_extraction",imin_GLO)  ; call erreur(status,.TRUE.,"put_att_GLOBAL3")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imax_extraction",imax_GLO)  ; call erreur(status,.TRUE.,"put_att_GLOBAL4")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmin_extraction",jmin_GLO)  ; call erreur(status,.TRUE.,"put_att_GLOBAL5")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmax_extraction",jmax_GLO)  ; call erreur(status,.TRUE.,"put_att_GLOBAL6")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidM,nav_lat_ID,gphit_REG); call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,glamt_REG); call erreur(status,.TRUE.,"var_nav_lon_ID")
if ( ln_isfcav ) then
  status = NF90_PUT_VAR(fidM,isf_draft_ID,isf_draft_REG);           call erreur(status,.TRUE.,"var_isf_draft_ID")
  status = NF90_PUT_VAR(fidM,Bathymetry_isf_ID,Bathymetry_isf_REG); call erreur(status,.TRUE.,"var_Bathymetry_isf_ID")
endif
status = NF90_PUT_VAR(fidM,Bathymetry_ID,Bathymetry_REG); call erreur(status,.TRUE.,"var_Bathymetry_ID")

status = NF90_CLOSE(fidM)                    
call erreur(status,.TRUE.,"final")         

write(*,*) '[done]'

end program modif


!=======================================================
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
    WRITE(*,*) 'ERROR:   ', iret
    message=NF90_STRERROR(iret)
    WRITE(*,*) 'WHICH MEANS: ',TRIM(message)
    IF ( lstop ) STOP
  ENDIF
  !
END SUBROUTINE erreur
