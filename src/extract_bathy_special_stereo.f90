program modif                                         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, IGE-CNRS, Jan. 2017
!
! Script to extract the bathymetry from a stereographic dataset (e.g. BEDMAP).
!
! The bathymetry along the boundaries (over a NINT(ai*1.5)-pts halo) is the same as in the
! dataset used as lateral boundary conditions (referred to as "CRS").
!
! 0- Initializations 
! 1- Read stereographic bathymetry and ice shelf draft
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

!-- STEREO variables :
INTEGER :: fidSTEREO1, fidSTEREO2, x_ID, y_ID, my_STEREO, mx_STEREO, bedrock_topography_ID, lat_ID, lon_ID, &
&          bathy_STEREO_ID, isf_draft_STEREO_ID, rq, thickness_STEREO_ID, surface_STEREO_ID

REAL*8 :: a, e, lat_c, pm, lon_0, lon_STEREO, chi, m_c, t_c, t, x, y, res_STEREO, res_REG

REAL*4,ALLOCATABLE,DIMENSION(:) :: x_STEREO, y_STEREO

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: bathy_STEREO, isf_draft_STEREO, surface_STEREO, lat_STEREO, zlon_STEREO


!-- local variables :
INTEGER :: fidGLO, fidCRS, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, isf_draft_ID, Bathymetry_isf_ID, Bathymetry_ID,    &
&          my_GLO, mx_GLO,  my_CRS, mx_CRS,  my_REG, mx_REG, imin_GLO, imax_GLO, jmin_GLO, jmax_GLO, npts, jtmp, ai, aj, bi, bj,      &
&          fidCOORDreg, fidCOORDpar, minlon, maxlon, minlat, maxlat, imin_STEREO, imax_STEREO, jmin_STEREO, jmax_STEREO, iREG, jREG,  &
&          iSTEREO, jSTEREO, iREGm1, iREGp1, jREGm1, jREGp1, kk, mx_tmp, my_tmp, i0, j0, rs, imin_STEREO2, imax_STEREO2, Nbox,        &
&          jmin_STEREO2, jmax_STEREO2

CHARACTER(LEN=150) :: file_bathy_out, file_in_coord_REG

INTEGER, ALLOCATABLE,DIMENSION(:,:) :: nn

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: gphit_REG, glamt_REG, zglamt_REG, isf_draft_REG, Bathymetry_isf_REG, Bathymetry_REG, &
&                                          Bathymetry_CRS, Bathymetry_isf_CRS, isf_draft_CRS, Bathymetry_0, Bathymetry_isf_0, isf_draft_0

INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: imin, imax, jmin, jmax

!-- Regional initial state
REAL*8                                 :: eps, pi, rad2deg, deg2rad
 
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

pi = ABS(ACOS(-1.d0))
deg2rad=pi/180.d0
rad2deg=180d0/pi

!=================================================================================
! 1- Read bathymetry and ice shelf draft on the stereographic grid
!=================================================================================

write(*,*) 'Reading bathy in ', TRIM(file_spe_bathy)

status = NF90_OPEN(TRIM(file_spe_bathy),0,fidSTEREO1); call erreur(status,.TRUE.,"Sart read STEREO") 

status = NF90_INQ_DIMID(fidSTEREO1,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidSTEREO1,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")

status = NF90_INQUIRE_DIMENSION(fidSTEREO1,dimID_y,len=my_STEREO); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidSTEREO1,dimID_x,len=mx_STEREO); call erreur(status,.TRUE.,"inq_dim_x")

ALLOCATE(  bathy_STEREO     (mx_STEREO,my_STEREO)  )
ALLOCATE(  isf_draft_STEREO (mx_STEREO,my_STEREO)  ) 
ALLOCATE(  x_STEREO (mx_STEREO)  ) 
ALLOCATE(  y_STEREO (my_STEREO)  ) 

status = NF90_INQ_VARID(fidSTEREO1,"bedrock",bathy_STEREO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO1,"bed",bathy_STEREO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO1,"BED",bathy_STEREO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO1,"BEDROCK",bathy_STEREO_ID)
call erreur(status,.TRUE.,"inq_bathy_STEREO_ID")
status = NF90_INQ_VARID(fidSTEREO1,"y",y_ID);                   call erreur(status,.TRUE.,"inq_lat_ID")
status = NF90_INQ_VARID(fidSTEREO1,"x",x_ID);                   call erreur(status,.TRUE.,"inq_lon_ID")

status = NF90_GET_VAR(fidSTEREO1,bathy_STEREO_ID,bathy_STEREO); call erreur(status,.TRUE.,"getvar_bathy_STEREO")
status = NF90_GET_VAR(fidSTEREO1,y_ID,y_STEREO);                call erreur(status,.TRUE.,"getvar_lat")
status = NF90_GET_VAR(fidSTEREO1,x_ID,x_STEREO);                call erreur(status,.TRUE.,"getvar_lon")

status = NF90_CLOSE(fidSTEREO1); call erreur(status,.TRUE.,"End read bathy STEREO")

!-----

write(*,*) 'Converting x,y to lon,lat'

a     = 6378137.00000000
e     =       0.08181919
lat_c =     -71.00000000
lon_0 =       0.d0

write(*,*) '        Earth Radius                          = ', a
write(*,*) '        Earth misshapenness (excentricity)    = ', e
write(*,*) '        latitude of true scale in degrees     = ', lat_c
write(*,*) '        meridian in along the positive Y axis = ', lon_0 

!- if the standard parallel is in S.Hemi., switch signs.
if ( lat_c .lt. 0.0 ) then
  pm     = -1.d0    ! plus or minus, north lat. or south
  lat_c  = -lat_c
  lon_0  = -lon_0
else
  pm     = 1.d0
endif

!- convert to radians :
lat_c    = deg2rad * lat_c
lon_0    = deg2rad * lon_0

ALLOCATE( lat_STEREO(mx_STEREO,my_STEREO), zlon_STEREO(mx_STEREO,my_STEREO) )

do iSTEREO=1,mx_STEREO
write(*,*) iSTEREO
do jSTEREO=1,my_STEREO

  if ( pm .lt. 0.0 ) then
    x      = -x_STEREO(iSTEREO)
    y      = -y_STEREO(jSTEREO)
  else
    x      =  x_STEREO(iSTEREO)
    y      =  y_STEREO(jSTEREO)
  endif

  !- See Snyder for details.
  t_c  = tan(pi/4-lat_c/2) / ( (1-e*sin(lat_c)) / (1+e*sin(lat_c)))**(e/2)
  m_c  = cos(lat_c) / sqrt( 1 - e**2 * (sin(lat_c))**2 )
  t    = sqrt(x**2+y**2) * t_c / ( a * m_c )

  chi = 0.5*pi - 2 * atan(t) !- find lat with a series instead of iterating.

  lat_STEREO(iSTEREO,jSTEREO) = chi + ( (1./2.) * e**2 + (5./24.) * e**4 + ( 1./ 12.) * e**6 + (  13./   360.) * e**8 ) * sin(2*chi) &
  &                                 + (                  (7./48.) * e**4 + (29./240.) * e**6 + ( 811./ 11520.) * e**8 ) * sin(4*chi) &
  &                                 + (                                  + ( 7./120.) * e**6 + (  81./  1120.) * e**8 ) * sin(6*chi) &
  &                                 + (                                                      + (4279./161280.) * e**8 ) * sin(8*chi)

  lon_STEREO = lon_0 + atan2(x,-y)

  !- correct the signs and phasing :
  lat_STEREO(iSTEREO,jSTEREO) = pm * lat_STEREO(iSTEREO,jSTEREO)
  lon_STEREO                  = pm * lon_STEREO
  lon_STEREO                  = mod(lon_STEREO+pi,2*pi)-pi !- want longitude in the range -pi to pi

  !- convert back to degrees :
  lat_STEREO (iSTEREO,jSTEREO) = rad2deg * lat_STEREO(iSTEREO,jSTEREO)
  lon_STEREO                   = rad2deg * lon_STEREO

  zlon_STEREO(iSTEREO,jSTEREO) = lon_STEREO
  if ( ln_dateline ) then
    if( lon_STEREO .lt. 0.e0 ) then
      zlon_STEREO(iSTEREO,jSTEREO) = 360.e0 + lon_STEREO
    endif
  else
    if( lon_STEREO .gt. 180.e0 ) then
      zlon_STEREO(iSTEREO,jSTEREO) = lon_STEREO - 360.e0
    endif
  endif

enddo
enddo

!check :
write(*,*) 'min lon = ', MINVAL(zlon_STEREO)
write(*,*) 'max lon = ', MAXVAL(zlon_STEREO)
write(*,*) 'min lat = ', MINVAL( lat_STEREO)
write(*,*) 'max lat = ', MAXVAL( lat_STEREO)

!-----

write(*,*) 'Reading ice shelf draft in ', TRIM(file_spe_isf_draft)

status = NF90_OPEN(TRIM(file_spe_isf_draft),0,fidSTEREO2); call erreur(status,.TRUE.,"read isf_draft STEREO") 
status = NF90_INQ_VARID(fidSTEREO2,"ice_draft",isf_draft_STEREO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"draft",isf_draft_STEREO_ID)
if ( status .ne. 0 ) then
  ! if ice draft not found, looking for surface and thickness :
  write(*,*) '  ... no ice-draft => reading thickness and surface height'
  ALLOCATE( surface_STEREO(mx_STEREO,my_STEREO) )
  status = NF90_INQ_VARID(fidSTEREO2,"thickness",thickness_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"thick",thickness_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"THICK",thickness_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"THICKNESS",thickness_STEREO_ID)
  call erreur(status,.TRUE.,"read thickness STEREO")
  status = NF90_INQ_VARID(fidSTEREO2,"surface",surface_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"surf",surface_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"SURF",surface_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"SURFACE",surface_STEREO_ID)
  call erreur(status,.TRUE.,"read surface STEREO")
  status = NF90_GET_VAR(fidSTEREO2,thickness_STEREO_ID,isf_draft_STEREO); call erreur(status,.TRUE.,"getvar_thickness_STEREO")
  status = NF90_GET_VAR(fidSTEREO2,surface_STEREO_ID,surface_STEREO); call erreur(status,.TRUE.,"getvar_surface_STEREO")
  isf_draft_STEREO(:,:) = surface_STEREO(:,:) - isf_draft_STEREO(:,:)
  DEALLOCATE( surface_STEREO )
else
  write(*,*) '  ... reading ice draft'
  status = NF90_GET_VAR(fidSTEREO2,isf_draft_STEREO_ID,isf_draft_STEREO); call erreur(status,.TRUE.,"getvar_isf_draft_STEREO")
endif
status = NF90_CLOSE(fidSTEREO2); call erreur(status,.TRUE.,"End read isf_draft STEREO")     


!=================================================================================
! 2- Read grid correspondance with GLO (i.e. extraction coordinates)
!=================================================================================

write(*,*) 'Reading REGIONAL lon,lat in ', TRIM(file_in_coord_REG)

status = NF90_OPEN(TRIM(file_in_coord_REG),0,fidCOORDreg); call erreur(status,.TRUE.,"read coord input")

status = NF90_INQ_DIMID(fidCOORDreg,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidCOORDreg,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
                                                     
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

status = NF90_CLOSE(fidCOORDreg); call erreur(status,.TRUE.,"end read fidCOORDreg")

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

status = NF90_INQ_VARID(fidCRS,"Bathymetry",Bathymetry_ID); call erreur(status,.TRUE.,"inq_Bathymetry_ID_CRS")
status = NF90_GET_VAR(fidCRS,Bathymetry_ID,Bathymetry_CRS); call erreur(status,.TRUE.,"getvar_Bathymetry_CRS")

if ( ln_isfcav ) then
  status = NF90_INQ_VARID(fidCRS,"Bathymetry_isf",Bathymetry_isf_ID)
  if ( status .eq. 0 ) then
    status = NF90_GET_VAR(fidCRS,Bathymetry_isf_ID,Bathymetry_isf_CRS); call erreur(status,.TRUE.,"getvar_Bathymetry_isf_CRS")
  else
    Bathymetry_isf_CRS(:,:) = Bathymetry_CRS(:,:)
  endif
  status = NF90_INQ_VARID(fidCRS,"isf_draft",isf_draft_ID)
  if ( status .eq. 0 ) then
    status = NF90_GET_VAR(fidCRS,isf_draft_ID,isf_draft_CRS); call erreur(status,.TRUE.,"getvar_isf_draft_CRS")
  else
    isf_draft_CRS(:,:) = 0.e0
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

res_REG = a * deg2rad * ABS( gphit_REG(1,2) - gphit_REG(1,1) ) 
res_STEREO = ABS( y_STEREO(2) - y_STEREO(1) )
rq = MIN( mx_STEREO , MIN( my_STEREO, 20 * CEILING( res_REG / res_STEREO ) ) )
write(*,*) ' '
write(*,*) 'Once a correspondance (iSTEREO,jSTEREO) on the stereographic grid is found for (iREG,jREG),'
write(*,*) ' the next stereographic point of the regional grid is searched in a square of'
write(*,*) ' radius ', rq, ' around (iSTEREO,jSTEREO).'
write(*,*) ' '

imin_STEREO2 = 1 
imax_STEREO2 = mx_STEREO
jmin_STEREO2 = 1
jmax_STEREO2 = my_STEREO

do iREG=1,mx_REG
!write(*,*) 'iREG = ', iREG
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

  ! restrict search for corresponding points on the stereographic grid (RAZ after a column)
  if ( ( jREG .eq. 1 ) .or. ( nn_perio .eq. 1 .and. ( iREG .lt. 11 .or. iREG .gt. mx_REG - 10 ) ) ) then
    imin_STEREO = 1
    imax_STEREO = mx_STEREO
    jmin_STEREO = 1
    jmax_STEREO = my_STEREO
  else
    imin_STEREO = MAX(         1, imin_STEREO2 )
    imax_STEREO = MIN( mx_STEREO, imax_STEREO2 )
    jmin_STEREO = MAX(         1, jmin_STEREO2 )
    jmax_STEREO = MIN( my_STEREO, jmax_STEREO2 )
  endif

  do iSTEREO=imin_STEREO,imax_STEREO
  do jSTEREO=jmin_STEREO,jmax_STEREO
 
    !-NB: we don't care too much about i=1 and i=mx (same for j) because it is masked anyway...    
    if (       zlon_STEREO(iSTEREO,jSTEREO) .ge. zglamt_REG(iREG,jREG) - 0.5*(zglamt_REG(iREG,jREG)-zglamt_REG(iREGm1,jREG)) &
    &    .and. zlon_STEREO(iSTEREO,jSTEREO) .lt. zglamt_REG(iREG,jREG) + 0.5*(zglamt_REG(iREGp1,jREG)-zglamt_REG(iREG,jREG)) &
    &    .and.  lat_STEREO(iSTEREO,jSTEREO) .ge.  gphit_REG(iREG,jREG) - 0.5*( gphit_REG(iREG,jREG)- gphit_REG(iREG,jREGm1)) &
    &    .and.  lat_STEREO(iSTEREO,jSTEREO) .lt.  gphit_REG(iREG,jREG) + 0.5*( gphit_REG(iREG,jREGp1)- gphit_REG(iREG,jREG)) ) then

      if ( nn(iREG,jREG) .eq. 0 ) then
        imin_STEREO2 = mx_STEREO ! starting value for the MIN
        imax_STEREO2 = 1         ! starting value for the MAX
        jmin_STEREO2 = my_STEREO
        jmax_STEREO2 = 1
      endif

      Bathymetry_isf_REG (iREG,jREG) = Bathymetry_isf_REG (iREG,jREG) - bathy_STEREO     (iSTEREO,jSTEREO)
      isf_draft_REG      (iREG,jREG) = isf_draft_REG      (iREG,jREG) - isf_draft_STEREO (iSTEREO,jSTEREO)

      nn(iREG,jREG) = nn(iREG,jREG) + 1

      imin_STEREO2 = MIN( imin_STEREO2, iSTEREO - rq )
      imax_STEREO2 = MAX( imax_STEREO2, iSTEREO + rq )
      jmin_STEREO2 = MIN( jmin_STEREO2, jSTEREO - rq )
      jmax_STEREO2 = MAX( jmax_STEREO2, jSTEREO + rq )

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
    
    !=== put exactly CRS data over a npts-point halo :
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

    !=== smooth transition from CRS to the interpolated fields (over npts points again) :
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

    !----- put CRS along modified North-Western corner :
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
        isf_draft_REG     (iREG,jREG) = ( (jREG-jmin(8)+2*npts+1) * isf_draft_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jREG) * isf_draft_0 ( iREG, jREG ) ) / (npts+1)                     
        Bathymetry_isf_REG(iREG,jREG) = ( (jREG-jmin(8)+2*npts+1) * Bathymetry_isf_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jREG) * Bathymetry_isf_0 ( iREG, jREG ) ) / (npts+1)                     
        Bathymetry_REG    (iREG,jREG) = ( (jREG-jmin(8)+2*npts+1) * Bathymetry_CRS ( NINT(FLOAT(iREG+imin_GLO-1-bi)/ai), NINT(FLOAT(jREG+jmin_GLO-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jREG) * Bathymetry_0 ( iREG, jREG ) ) / (npts+1)                     
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

    ! correction to avoid isolated open ocean:
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
    write(*,*) '   imin, imax = ', imin(:), imax(:)
    write(*,*) '   jmin, jmax = ', jmin(:), jmax(:)

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

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bathy_specila_stereo.f90")
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
