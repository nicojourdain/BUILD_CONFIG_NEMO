program modif                                         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Script to extract the bathymetry of the child domain (CHLD) from a stereographic dataset
! (e.g. BEDMAP, BedMachine).
!
! The bathymetry along the boundaries (over a NINT(ai*1.5)-pts halo) is the same as in the
! dataset used as lateral boundary conditions (referred to as "PAR" for parent).
!
! 0- Initializations 
! 1- Read stereographic bathymetry and ice shelf draft
! 2- Read grid correspondance with EXT (i.e. extraction coordinates)
! 3- Read coarse bathymetry ("PAR") used for consistent bathymetry along boundaries
! 4- Calculate bathy/isf draft on the CHLD grid
! 5- Writing new CHLD bathymetry file
!
! History: - Jan. 2017: Initial version (N. Jourdain, CNRS-IGE)
!          - Jan. 2022: cleaning and new naming conventions (EXT/PAR/CHLD)
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

REAL*8 :: a, e, lat_c, pm, lon_0, lon_STEREO, chi, m_c, t_c, t, x, y, res_STEREO, res_CHLD

REAL*4,ALLOCATABLE,DIMENSION(:) :: x_STEREO, y_STEREO

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: bathy_STEREO, isf_draft_STEREO, surface_STEREO, lat_STEREO, zlon_STEREO


!-- local variables :
INTEGER :: fidEXT, fidPAR, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, isf_draft_ID, Bathymetry_isf_ID, Bathymetry_ID,    &
&          my_EXT, mx_EXT,  my_PAR, mx_PAR,  my_CHLD, mx_CHLD, imin_EXT, imax_EXT, jmin_EXT, jmax_EXT, npts, jtmp, ai, aj, bi, bj,      &
&          fidCOORDchld, fidCOORDpar, minlon, maxlon, minlat, maxlat, imin_STEREO, imax_STEREO, jmin_STEREO, jmax_STEREO, iCHLD, jCHLD,  &
&          iSTEREO, jSTEREO, iCHLDm1, iCHLDp1, jCHLDm1, jCHLDp1, kk, mx_tmp, my_tmp, i0, j0, rs, imin_STEREO2, imax_STEREO2, Nbox,        &
&          jmin_STEREO2, jmax_STEREO2

CHARACTER(LEN=150) :: file_bathy_out, file_in_coord_CHLD

INTEGER, ALLOCATABLE,DIMENSION(:,:) :: nn

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: gphit_CHLD, glamt_CHLD, zglamt_CHLD, isf_draft_CHLD, Bathymetry_isf_CHLD, Bathymetry_CHLD, &
&                                          Bathymetry_PAR, Bathymetry_isf_PAR, isf_draft_PAR, Bathymetry_0, Bathymetry_isf_0, isf_draft_0

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

! name of child bathymetry file (output file) :
write(file_bathy_out,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/bathy_meter_',a,'.nc')

! name of child coordinates file (output file) :
write(file_in_coord_CHLD,102) TRIM(config_dir), TRIM(config)
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
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"ice_base_topography",isf_draft_STEREO_ID)
if ( status .ne. 0 ) then
  ! if ice draft not found, looking for surface and thickness :
  write(*,*) '  ... no ice-draft => reading thickness and surface height'
  ALLOCATE( surface_STEREO(mx_STEREO,my_STEREO) )
  status = NF90_INQ_VARID(fidSTEREO2,"thickness",thickness_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"thi",thickness_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"thick",thickness_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"THICK",thickness_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"THICKNESS",thickness_STEREO_ID)
  call erreur(status,.TRUE.,"read thickness STEREO")
  status = NF90_INQ_VARID(fidSTEREO2,"surface",surface_STEREO_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSTEREO2,"srf",surface_STEREO_ID)
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
! 2- Read grid correspondance with EXT (i.e. extraction coordinates)
!=================================================================================

write(*,*) 'Reading CHLDIONAL lon,lat in ', TRIM(file_in_coord_CHLD)

status = NF90_OPEN(TRIM(file_in_coord_CHLD),0,fidCOORDchld); call erreur(status,.TRUE.,"read coord input")

status = NF90_INQ_DIMID(fidCOORDchld,"x",dimID_x)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidCOORDchld,"X",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidCOORDchld,"y",dimID_y)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidCOORDchld,"Y",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")
                                                     
status = NF90_INQUIRE_DIMENSION(fidCOORDchld,dimID_y,len=my_CHLD); call erreur(status,.TRUE.,"inq_dim_y_CHLD")
status = NF90_INQUIRE_DIMENSION(fidCOORDchld,dimID_x,len=mx_CHLD); call erreur(status,.TRUE.,"inq_dim_x_CHLD")

ALLOCATE(  gphit_CHLD (mx_CHLD,my_CHLD)  )
ALLOCATE(  glamt_CHLD (mx_CHLD,my_CHLD)  )
ALLOCATE( zglamt_CHLD (mx_CHLD,my_CHLD)  )

status = NF90_INQ_VARID(fidCOORDchld,"gphit",nav_lat_ID); call erreur(status,.TRUE.,"inq_gphit_CHLD_ID")
status = NF90_INQ_VARID(fidCOORDchld,"glamt",nav_lon_ID); call erreur(status,.TRUE.,"inq_glamt_CHLD_ID")

status = NF90_GET_VAR(fidCOORDchld,nav_lat_ID,gphit_CHLD); call erreur(status,.TRUE.,"getvar_gphit_CHLD")
status = NF90_GET_VAR(fidCOORDchld,nav_lon_ID,glamt_CHLD); call erreur(status,.TRUE.,"getvar_glamt_CHLD")

status = NF90_GET_ATT(fidCOORDchld, NF90_GLOBAL, "ai", ai); call erreur(status,.TRUE.,"read att1")
status = NF90_GET_ATT(fidCOORDchld, NF90_GLOBAL, "bi", bi); call erreur(status,.TRUE.,"read att2")
status = NF90_GET_ATT(fidCOORDchld, NF90_GLOBAL, "aj", aj); call erreur(status,.TRUE.,"read att3")
status = NF90_GET_ATT(fidCOORDchld, NF90_GLOBAL, "bj", bj); call erreur(status,.TRUE.,"read att4")
status = NF90_GET_ATT(fidCOORDchld, NF90_GLOBAL, "imin_extraction", imin_EXT); call erreur(status,.TRUE.,"read att5")
status = NF90_GET_ATT(fidCOORDchld, NF90_GLOBAL, "jmin_extraction", jmin_EXT); call erreur(status,.TRUE.,"read att6")
status = NF90_GET_ATT(fidCOORDchld, NF90_GLOBAL, "imax_extraction", imax_EXT); call erreur(status,.TRUE.,"read att7")
status = NF90_GET_ATT(fidCOORDchld, NF90_GLOBAL, "jmax_extraction", jmax_EXT); call erreur(status,.TRUE.,"read att8")

status = NF90_CLOSE(fidCOORDchld); call erreur(status,.TRUE.,"end read fidCOORDchld")

zglamt_CHLD(:,:)=glamt_CHLD(:,:)
if ( ln_dateline ) then
  where( glamt_CHLD(:,:) .lt. 0.e0 )
    zglamt_CHLD(:,:)=glamt_CHLD(:,:)+360.e0
  endwhere
else
  where( glamt_CHLD(:,:) .gt. 180.e0 )
    zglamt_CHLD(:,:)=glamt_CHLD(:,:)-360.e0
  endwhere
endif

!=================================================================================
! 3- Read coarse bathymetry used for consistent bathymetry along boundaries
!=================================================================================

write(*,*) 'Reading coarse bathymetry for consistent boundaries: ', TRIM(file_in_bathy_bdy)

status = NF90_OPEN(TRIM(file_in_bathy_bdy),0,fidPAR); call erreur(status,.TRUE.,"read_coarse_bathymetry") 

status = NF90_INQ_DIMID(fidPAR,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_PAR")
status = NF90_INQ_DIMID(fidPAR,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_PAR")

status = NF90_INQUIRE_DIMENSION(fidPAR,dimID_y,len=my_PAR); call erreur(status,.TRUE.,"inq_dim_y_PAR")
status = NF90_INQUIRE_DIMENSION(fidPAR,dimID_x,len=mx_PAR); call erreur(status,.TRUE.,"inq_dim_x_PAR")

ALLOCATE(  Bathymetry_PAR (mx_PAR,my_PAR)  ) 
if ( ln_isfcav ) then
  ALLOCATE(  Bathymetry_isf_PAR (mx_PAR,my_PAR)  )
  ALLOCATE(  isf_draft_PAR (mx_PAR,my_PAR)  )
endif

status = NF90_INQ_VARID(fidPAR,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"latitude",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID_PAR")
status = NF90_INQ_VARID(fidPAR,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"longitude",nav_lon_ID)
call erreur(status,.TRUE.,"inq_nav_lon_ID_PAR")

!NB: here we try to read two bathymetry variables: one without ice shelf cavities ("Bathymetry")
!    and, if ln_isfcav is true,  one with ice shelf cavities ("Bathymetry_isf").
!    If no "Bathymetry_isf" is found, we use a single bathymetry variable.
status = NF90_INQ_VARID(fidPAR,"Bathymetry",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"bathymetry",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"bathy_metry",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"bathy",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"Bathy",Bathymetry_ID)
call erreur(status,.TRUE.,"inq_Bathymetry_ID_PAR")
status = NF90_GET_VAR(fidPAR,Bathymetry_ID,Bathymetry_PAR); call erreur(status,.TRUE.,"getvar_Bathymetry_PAR")
if ( ln_isfcav ) then
  status = NF90_INQ_VARID(fidPAR,"Bathymetry_isf",Bathymetry_isf_ID)
  if ( status .eq. 0 ) then
    status = NF90_GET_VAR(fidPAR,Bathymetry_isf_ID,Bathymetry_isf_PAR); call erreur(status,.TRUE.,"getvar_Bathymetry_isf_PAR")
  else
    Bathymetry_isf_PAR(:,:) = Bathymetry_PAR(:,:)
  endif
  status = NF90_INQ_VARID(fidPAR,"isf_draft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"isfdraft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"draft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"ice_shelf_draft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"Ice_shelf_draft",isf_draft_ID)
  if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"ice_shelf_base",isf_draft_ID)
  if ( status .eq. 0 ) then
    status = NF90_GET_VAR(fidPAR,isf_draft_ID,isf_draft_PAR); call erreur(status,.TRUE.,"getvar_isf_draft_PAR")
  else
    isf_draft_PAR(:,:) = 0.e0
     write(*,*) 'WARNING: no ice shelf draft found. Allowed names are :'
     write(*,*) '         "isf_draft", "isfdraft", "draft", "ice_shelf_draft", "Ice_shelf_draft", "ice_shelf_base"'
     write(*,*) '         >>>>> ASSUMING ice_draft = 0 for the COARSE data'
  endif
endif

status = NF90_CLOSE(fidPAR); call erreur(status,.TRUE.,"close_coarse_bathy_file")

!===================================================================
! put COARSE grid in a npts-pts halo (+ transition in another halo)

npts=CEILING(ai*1.5) 
write(*,*) 'put COARSE grid bathymetry in a npts-pts halo, with npts = ', npts

!=================================================================================
! 4- Calculate bathy/isf draft on the CHLD grid :
!=================================================================================

write(*,*) 'Calculating bathymetry on the child grid...'
    
ALLOCATE(  nn                 (mx_CHLD,my_CHLD)  )
ALLOCATE(  isf_draft_CHLD      (mx_CHLD,my_CHLD)  )
ALLOCATE(  Bathymetry_isf_CHLD (mx_CHLD,my_CHLD)  )
ALLOCATE(  Bathymetry_CHLD     (mx_CHLD,my_CHLD)  )
ALLOCATE(  isf_draft_0        (mx_CHLD,my_CHLD)  )
ALLOCATE(  Bathymetry_isf_0   (mx_CHLD,my_CHLD)  )
ALLOCATE(  Bathymetry_0       (mx_CHLD,my_CHLD)  )
  
nn(:,:) = 0 
Bathymetry_isf_CHLD (:,:) = 0.e0
isf_draft_CHLD      (:,:) = 0.e0
Bathymetry_CHLD     (:,:) = 0.e0

res_CHLD = a * deg2rad * ABS( gphit_CHLD(1,2) - gphit_CHLD(1,1) ) 
res_STEREO = ABS( y_STEREO(2) - y_STEREO(1) )
rq = MIN( mx_STEREO , MIN( my_STEREO, 20 * CEILING( res_CHLD / res_STEREO ) ) )
write(*,*) ' '
write(*,*) 'Once a correspondance (iSTEREO,jSTEREO) on the stereographic grid is found for (iCHLD,jCHLD),'
write(*,*) ' the next stereographic point of the child grid is searched in a square of'
write(*,*) ' radius ', rq, ' around (iSTEREO,jSTEREO).'
write(*,*) ' '

imin_STEREO2 = 1 
imax_STEREO2 = mx_STEREO
jmin_STEREO2 = 1
jmax_STEREO2 = my_STEREO

do iCHLD=1,mx_CHLD
!write(*,*) 'iCHLD = ', iCHLD
do jCHLD=1,my_CHLD

  if ( nn_perio .eq. 1 ) then !- periodic
    if    ( iCHLD+1 .gt. mx_CHLD ) then
      iCHLDp1 = 3
      iCHLDm1 = iCHLD - 1
    elseif ( iCHLD-1 .lt. 1     ) then
      iCHLDp1 = iCHLD + 1
      iCHLDm1 = mx_CHLD-2
    else
      iCHLDp1 = iCHLD + 1
      iCHLDm1 = iCHLD - 1
    endif
  elseif ( nn_perio .eq. 0 ) then
    if    ( iCHLD+1 .gt. mx_CHLD ) then
      iCHLDp1 = iCHLD
      iCHLDm1 = iCHLD - 1
    elseif ( iCHLD-1 .lt. 1     ) then
      iCHLDp1 = iCHLD + 1
      iCHLDm1 = iCHLD
    else
      iCHLDp1 = iCHLD + 1
      iCHLDm1 = iCHLD - 1
    endif
  else
    write(*,*) '~!@#$%^* nn_perio must be either 0 or 1 >>>>> stop !!'
    stop
  endif
  !-
  if    ( jCHLD+1 .gt. my_CHLD ) then
    jCHLDp1 = jCHLD
    jCHLDm1 = jCHLD - 1
  elseif ( iCHLD-1 .lt. 1     ) then
    jCHLDp1 = jCHLD + 1
    jCHLDm1 = jCHLD
  else
    jCHLDp1 = jCHLD + 1
    jCHLDm1 = jCHLD - 1
  endif
  !-
  iCHLDm1 = MAX( MIN( iCHLDm1, mx_CHLD ), 1 )
  iCHLDp1 = MAX( MIN( iCHLDp1, mx_CHLD ), 1 )
  jCHLDm1 = MAX( MIN( jCHLDm1, my_CHLD ), 1 )
  jCHLDp1 = MAX( MIN( jCHLDp1, my_CHLD ), 1 )
  !-

  ! restrict search for corresponding points on the stereographic grid (RAZ after a column)
  if ( ( jCHLD .eq. 1 ) .or. ( nn_perio .eq. 1 .and. ( iCHLD .lt. 11 .or. iCHLD .gt. mx_CHLD - 10 ) ) ) then
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
    if (       zlon_STEREO(iSTEREO,jSTEREO) .ge. zglamt_CHLD(iCHLD,jCHLD) - 0.5*(zglamt_CHLD(iCHLD,jCHLD)-zglamt_CHLD(iCHLDm1,jCHLD)) &
    &    .and. zlon_STEREO(iSTEREO,jSTEREO) .lt. zglamt_CHLD(iCHLD,jCHLD) + 0.5*(zglamt_CHLD(iCHLDp1,jCHLD)-zglamt_CHLD(iCHLD,jCHLD)) &
    &    .and.  lat_STEREO(iSTEREO,jSTEREO) .ge.  gphit_CHLD(iCHLD,jCHLD) - 0.5*( gphit_CHLD(iCHLD,jCHLD)- gphit_CHLD(iCHLD,jCHLDm1)) &
    &    .and.  lat_STEREO(iSTEREO,jSTEREO) .lt.  gphit_CHLD(iCHLD,jCHLD) + 0.5*( gphit_CHLD(iCHLD,jCHLDp1)- gphit_CHLD(iCHLD,jCHLD)) ) then

      if ( nn(iCHLD,jCHLD) .eq. 0 ) then
        imin_STEREO2 = mx_STEREO ! starting value for the MIN
        imax_STEREO2 = 1         ! starting value for the MAX
        jmin_STEREO2 = my_STEREO
        jmax_STEREO2 = 1
      endif

      Bathymetry_isf_CHLD (iCHLD,jCHLD) = Bathymetry_isf_CHLD (iCHLD,jCHLD) - bathy_STEREO     (iSTEREO,jSTEREO)
      isf_draft_CHLD      (iCHLD,jCHLD) = isf_draft_CHLD      (iCHLD,jCHLD) - isf_draft_STEREO (iSTEREO,jSTEREO)

      nn(iCHLD,jCHLD) = nn(iCHLD,jCHLD) + 1

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
  Bathymetry_isf_CHLD (:,:) = Bathymetry_isf_CHLD (:,:) / nn(:,:)
  isf_draft_CHLD      (:,:) = isf_draft_CHLD      (:,:) / nn(:,:)
elsewhere
  Bathymetry_isf_CHLD (:,:) = 99999999.9
  isf_draft_CHLD      (:,:) = 99999999.9
endwhere

!-----
! Interpolation at points with nn=0 :
! Only look for first neighbour
! needs to be adapted if many missing points

do iCHLD=1,mx_CHLD
do jCHLD=1,my_CHLD

  if ( nn(iCHLD,jCHLD) .eq. 0 ) then

    do rs=1,10
      if ( nn_perio .eq. 1 ) then !- periodic
        if    ( iCHLD+rs .gt. mx_CHLD ) then
          iCHLDp1 = 2+rs
        else
          iCHLDp1 = iCHLD + rs
        endif
      elseif ( nn_perio .eq. 0 ) then
        if    ( iCHLD+rs .gt. mx_CHLD ) then
          iCHLDp1 = mx_CHLD
          exit ! with nn(iCHLDp1,jCHLD)=0
        else
          iCHLDp1 = iCHLD + rs
        endif
      endif
      if ( nn(iCHLDp1,jCHLD) .ne. 0 ) exit ! with nn(iCHLDp1,jCHLD)>0
    enddo

    do rs=1,10
      if ( nn_perio .eq. 1 ) then !- periodic
        if ( iCHLD-rs .lt. 1     ) then
          iCHLDm1 = mx_CHLD-rs-1
        else
          iCHLDm1 = iCHLD - rs
        endif
      elseif ( nn_perio .eq. 0 ) then
        if ( iCHLD-rs .lt. 1     ) then
          iCHLDm1 = 1
          exit ! with nn(iCHLDm1,jCHLD)=0
        else
          iCHLDm1 = iCHLD - rs
        endif
      endif
      if ( nn(iCHLDm1,jCHLD) .gt. 0 ) exit ! with nn(iCHLDm1,jCHLD)>0
    enddo

    do rs=1,10
      if    ( jCHLD+rs .gt. my_CHLD ) then
        jCHLDp1 = my_CHLD
        exit ! with nn(iCHLD,jCHLDp1)=0
      else
        jCHLDp1 = jCHLD + rs
      endif
      if ( nn(iCHLD,jCHLDp1) .gt. 0 ) exit ! with nn(iCHLD,jCHLDp1)>0
    enddo

    do rs=1,10
      if    ( jCHLD-rs .lt. 1 ) then
        jCHLDm1 = 1
        exit ! with nn(iCHLD,jCHLDm1)=0
      else
        jCHLDm1 = jCHLD - rs
      endif
      if ( nn(iCHLD,jCHLDm1) .gt. 0 ) exit ! with nn(iCHLD,jCHLDm1)>0
    enddo

    Bathymetry_isf_CHLD(iCHLD,jCHLD) = (   nn(iCHLDm1,jCHLD  ) * (iCHLDp1-iCHLD) * Bathymetry_isf_CHLD(iCHLDm1,jCHLD  )   &
    &                                 + nn(iCHLDp1,jCHLD  ) * (iCHLD-iCHLDm1) * Bathymetry_isf_CHLD(iCHLDp1,jCHLD  )   &
    &                                 + nn(iCHLD  ,jCHLDm1) * (jCHLDp1-jCHLD) * Bathymetry_isf_CHLD(iCHLD  ,jCHLDm1)   &
    &                                 + nn(iCHLD  ,jCHLDp1) * (jCHLD-jCHLDm1) * Bathymetry_isf_CHLD(iCHLD  ,jCHLDp1) ) &
    &                             / (   nn(iCHLDm1,jCHLD  ) * (iCHLDp1-iCHLD)                                       &
    &                                 + nn(iCHLDp1,jCHLD  ) * (iCHLD-iCHLDm1)                                       &
    &                                 + nn(iCHLD  ,jCHLDm1) * (jCHLDp1-jCHLD)                                       &
    &                                 + nn(iCHLD  ,jCHLDp1) * (jCHLD-jCHLDm1)                                     )

    isf_draft_CHLD(iCHLD,jCHLD) = (   nn(iCHLDm1,jCHLD  ) * (iCHLDp1-iCHLD) * isf_draft_CHLD(iCHLDm1,jCHLD  )   &
    &                            + nn(iCHLDp1,jCHLD  ) * (iCHLD-iCHLDm1) * isf_draft_CHLD(iCHLDp1,jCHLD  )   &
    &                            + nn(iCHLD  ,jCHLDm1) * (jCHLDp1-jCHLD) * isf_draft_CHLD(iCHLD  ,jCHLDm1)   &
    &                            + nn(iCHLD  ,jCHLDp1) * (jCHLD-jCHLDm1) * isf_draft_CHLD(iCHLD  ,jCHLDp1) ) &
    &                        / (   nn(iCHLDm1,jCHLD  ) * (iCHLDp1-iCHLD)                                       &
    &                            + nn(iCHLDp1,jCHLD  ) * (iCHLD-iCHLDm1)                                       &
    &                            + nn(iCHLD  ,jCHLDm1) * (jCHLDp1-jCHLD)                                       &
    &                            + nn(iCHLD  ,jCHLDp1) * (jCHLD-jCHLDm1)                                     )

   endif

enddo
enddo

!---------

write(*,*) '>i=125 ', nn(125,20), Bathymetry_isf_CHLD(125,20), isf_draft_CHLD(125,20)
write(*,*) '>i=624 ', nn(624,20), Bathymetry_isf_CHLD(624,20), isf_draft_CHLD(624,20)

where( isf_draft_CHLD(:,:) .lt. 0.e0 )
  isf_draft_CHLD(:,:) = 0.e0
endwhere

where( Bathymetry_isf_CHLD(:,:) .lt. 0.e0 )
  Bathymetry_isf_CHLD(:,:) = 0.e0
endwhere

where( isf_draft_CHLD(:,:) .gt. 1.e0 .and. isf_draft_CHLD(:,:) .lt. 1.e4 )
  Bathymetry_CHLD(:,:) = 0.e0
elsewhere
  Bathymetry_CHLD(:,:) = Bathymetry_isf_CHLD(:,:)
endwhere

!write(*,*) '#i=125 ', Bathymetry_CHLD(125,20), Bathymetry_isf_CHLD(125,20), isf_draft_CHLD(125,20)
!write(*,*) '#i=624 ', Bathymetry_CHLD(624,20), Bathymetry_isf_CHLD(624,20), isf_draft_CHLD(624,20)

!---------------------------------------
! Adjust bathy along the edges of the CHLD grid :

write(*,*) 'Halo with bathy from the dataset used as lateral boundaries...'

if ( ln_coarse_bdy ) then
    
    !---------------------------------------
    ! Adjust bathy along the edges of the CHLD grid :
    
    write(*,*) 'Halo with bathy from parent grid...'
    
    !=== put exactly PAR data over a npts-point halo :
    do jCHLD=1,my_CHLD
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        !-- Western BDY :
        do iCHLD=1,npts
          if ( ln_isfcav) then
            isf_draft_CHLD     (iCHLD,jCHLD) = isf_draft_PAR      ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
            Bathymetry_isf_CHLD(iCHLD,jCHLD) = Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
          endif
          Bathymetry_CHLD    (iCHLD,jCHLD) = Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
        enddo
        !--- Eastern BDY :
        do iCHLD=mx_CHLD-npts+1,mx_CHLD
          if ( ln_isfcav) then
            isf_draft_CHLD     (iCHLD,jCHLD) = isf_draft_PAR      ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
            Bathymetry_isf_CHLD(iCHLD,jCHLD) = Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
          endif
          Bathymetry_CHLD    (iCHLD,jCHLD) = Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    !-- Southern BDY :
    do jCHLD=1,npts
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iCHLD=npts+1,mx_CHLD-npts
          if ( ln_isfcav) then
            isf_draft_CHLD     (iCHLD,jCHLD) = isf_draft_PAR      ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
            Bathymetry_isf_CHLD(iCHLD,jCHLD) = Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
          endif
          Bathymetry_CHLD    (iCHLD,jCHLD) = Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    !--- Northern BDY :
    do jCHLD=my_CHLD-npts+1,my_CHLD
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iCHLD=npts+1,mx_CHLD-npts
          if ( ln_isfcav) then
            isf_draft_CHLD     (iCHLD,jCHLD) = isf_draft_PAR      ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
            Bathymetry_isf_CHLD(iCHLD,jCHLD) = Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
          endif
          Bathymetry_CHLD    (iCHLD,jCHLD) = Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    
    write(*,*) 'Smooth transition...'

    Bathymetry_0     = Bathymetry_CHLD
    Bathymetry_isf_0 = Bathymetry_isf_CHLD
    isf_draft_0      = isf_draft_CHLD

    !=== smooth transition from PAR to the interpolated fields (over npts points again) :
    do jCHLD=npts+1,my_CHLD-npts
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      !-- Western BDY :
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iCHLD=npts+1,2*npts
          if ( ln_isfcav) then
            if ( isf_draft_0(iCHLD,jCHLD) .gt. 1.0 .and. isf_draft_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_CHLD(iCHLD,jCHLD) = (  (2*npts+1-iCHLD) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                           + (iCHLD-npts) * isf_draft_0(iCHLD,jCHLD) ) / (npts+1)
            !else
            !  isf_draft_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            if ( Bathymetry_isf_0(iCHLD,jCHLD) .gt. 1.0 .and. Bathymetry_isf_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_CHLD(iCHLD,jCHLD) = (  (2*npts+1-iCHLD) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                                + (iCHLD-npts) * Bathymetry_isf_0(iCHLD,jCHLD) ) / (npts+1)
            !else
            !  Bathymetry_isf_CHLD   (iCHLD,jCHLD) = 0.0
            endif
          endif
          Bathymetry_CHLD(iCHLD,jCHLD) = (   (2*npts+1-iCHLD) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
          &                             + (iCHLD-npts) * Bathymetry_0(iCHLD,jCHLD) ) / (npts+1)
        enddo
        !-- Eastern BDY
        do iCHLD=mx_CHLD-2*npts+1,mx_CHLD-npts
          if ( ln_isfcav) then
            if ( isf_draft_0(iCHLD,jCHLD) .gt. 1.0 .and. isf_draft_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp) .gt. 1.0 ) then
              isf_draft_CHLD(iCHLD,jCHLD) = (   (iCHLD-mx_CHLD+2*npts) * isf_draft_PAR( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                            + (mx_CHLD-npts+1-iCHLD) * isf_draft_0(iCHLD,jCHLD) ) / (npts+1)
            !else
            !  isf_draft_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            if ( Bathymetry_isf_0(iCHLD,jCHLD) .gt. 1.0 .and. Bathymetry_isf_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp) .gt. 1.0 ) then
              Bathymetry_isf_CHLD(iCHLD,jCHLD) = (   (iCHLD-mx_CHLD+2*npts) * Bathymetry_isf_PAR( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                                 + (mx_CHLD-npts+1-iCHLD) * Bathymetry_isf_0(iCHLD,jCHLD) ) / (npts+1)
            !else
            !  Bathymetry_isf_CHLD   (iCHLD,jCHLD) = 0.0
            endif
          endif
          Bathymetry_CHLD(iCHLD,jCHLD) = (   (iCHLD-mx_CHLD+2*npts) * Bathymetry_PAR( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
          &                             + (mx_CHLD-npts+1-iCHLD) * Bathymetry_0(iCHLD,jCHLD) ) / (npts+1)
        enddo
      endif
    enddo
    !-- Southern BDY :
    do jCHLD=npts+1,2*npts
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iCHLD=2*npts+1,mx_CHLD-2*npts
          if ( ln_isfcav) then
            if ( isf_draft_0(iCHLD,jCHLD) .gt. 1.0 .and. isf_draft_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_CHLD(iCHLD,jCHLD) = (   (2*npts+1-jCHLD) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                            + (jCHLD-npts) * isf_draft_0(iCHLD,jCHLD) ) / (npts+1)
            !else
            !  isf_draft_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            if ( Bathymetry_isf_0(iCHLD,jCHLD) .gt. 1.0 .and. Bathymetry_isf_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_CHLD(iCHLD,jCHLD) = (   (2*npts+1-jCHLD) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                                 + (jCHLD-npts) * Bathymetry_isf_0(iCHLD,jCHLD) ) / (npts+1)
            !else
            !  Bathymetry_isf_CHLD   (iCHLD,jCHLD) = 0.0
            endif
          endif
          Bathymetry_CHLD(iCHLD,jCHLD) = (   (2*npts+1-jCHLD) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
          &                             + (jCHLD-npts) * Bathymetry_0(iCHLD,jCHLD) ) / (npts+1)
        enddo
      endif
    enddo
    !-- Northern BDY :
    do jCHLD=my_CHLD-2*npts+1,my_CHLD-npts
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of the grid
        do iCHLD=2*npts+1,mx_CHLD-2*npts
          if ( ln_isfcav) then
            if ( isf_draft_0(iCHLD,jCHLD) .gt. 1.0 .and. isf_draft_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_CHLD(iCHLD,jCHLD) = (   (jCHLD-my_CHLD+2*npts) * isf_draft_PAR( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                            + (my_CHLD-npts+1-jCHLD) * isf_draft_0(iCHLD,jCHLD) ) / (npts+1)
            !else
            !  isf_draft_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            if ( Bathymetry_isf_0(iCHLD,jCHLD) .gt. 1.0 .and. Bathymetry_isf_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_CHLD(iCHLD,jCHLD) = (   (jCHLD-my_CHLD+2*npts) * Bathymetry_isf_PAR( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                                 + (my_CHLD-npts+1-jCHLD) * Bathymetry_isf_0(iCHLD,jCHLD) ) / (npts+1)
            !else
            !  Bathymetry_isf_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            Bathymetry_isf_CHLD(iCHLD,jCHLD) = (   (jCHLD-my_CHLD+2*npts) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
            &                                 + (my_CHLD-npts+1-jCHLD) * Bathymetry_isf_0(iCHLD,jCHLD) ) / (npts+1)
          endif
          Bathymetry_CHLD(iCHLD,jCHLD) = (   (jCHLD-my_CHLD+2*npts) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
          &                             + (my_CHLD-npts+1-jCHLD) * Bathymetry_0(iCHLD,jCHLD) ) / (npts+1)
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
    i0 = 2464 - imin_EXT
    j0 =  151 - jmin_EXT

    !! correction to avoid a closed cavity of 2x2x2 pts (identified after first
    !! mesh_mask creation)
    !isf_draft_CHLD     (i0+241:i0+242,j0+667:j0+668) = 0.0
    !Bathymetry_isf_CHLD(i0+241:i0+242,j0+667:j0+668) = 0.0
    !
    !! no isf along eastern boundary (adapt manually to adjust more accurately) :
    !isf_draft_CHLD     (i0+1095:mx_CHLD,j0+668:j0+703) = 0.0
    !Bathymetry_isf_CHLD(i0+1095:mx_CHLD,j0+668:j0+703) = Bathymetry_CHLD(i0+1095:mx_CHLD,j0+668:j0+703)

    ! boxes to fill the Bellingshausen Sea : filled over [imin:imax,jmin:my]
    Nbox = 8
    ALLOCATE( imin(Nbox), imax(Nbox), jmin(Nbox), jmax(Nbox) )
    imin = (/   1+i0 , 192+i0 , 213+i0 , 237+i0 , 254+i0 , 275+i0 , 287+i0 , 299+i0 /)
    imax = (/ 191+i0 , 212+i0 , 236+i0 , 253+i0 , 274+i0 , 286+i0 , 298+i0 , &
    &           NINT(FLOAT(325+i0+imin_EXT-1-bi)/ai)*ai+bi-imin_EXT+1-1 /) ! WARNING: last number must match with PAR (i.e. next unmasked point neads to be on PAR) !!
    jmin = (/ 494+j0 , 807+j0 , 835+j0 , 862+j0 , 876+j0 , 894+j0 ,            &
    &           NINT(FLOAT(899+j0+jmin_EXT-1-bj)/aj)*aj+bj-jmin_EXT+1+1 ,&
    &           NINT(FLOAT(903+j0+jmin_EXT-1-bj)/aj)*aj+bj-jmin_EXT+1+1 /) ! WARNING: last two numbers must match with PAR (i.e. next unmasked point neads to be on PAR) !!
    jmax(:) = my_CHLD

    do kk=1,Nbox
      imin(kk) = MIN( MAX( 1, imin(kk) ), mx_CHLD )
      imax(kk) = MIN( MAX( 1, imax(kk) ), mx_CHLD )
      jmin(kk) = MIN( MAX( 1, jmin(kk) ), my_CHLD )
      jmax(kk) = MIN( MAX( 1, jmax(kk) ), my_CHLD )
    enddo     
 
    write(*,*) 'Note for future bdy building:'
    write(*,*) '                    '
    write(*,*) '  ii_bdy_west(1)  = ', imax(8)+1
    write(*,*) '  j1_bdy_west(1)  = ', jmin(8)
    write(*,*) '  j2_bdy_west(1)  = ', my_CHLD-1 ! given that jmax(8)=my_CHLD-1
    write(*,*) '                    '
    write(*,*) '  i1_bdy_north(1) = ', imin(8)
    write(*,*) '  i2_bdy_north(1) = ', imax(8)
    write(*,*) '  jj_bdy_north(1) = ', jmin(8)-1
    write(*,*) '                    '
    write(*,*) '  i1_bdy_north(2) = ', imax(8)+2
    write(*,*) '  i2_bdy_north(2) = ', mx_CHLD-2 ! -2 because east bdy already contains mx_CHLD-1
    write(*,*) '  jj_bdy_north(2) = ', my_CHLD-1
    write(*,*) '                    '

    !----- put PAR along modified North-Western corner :
    !- bdy_west(1) :
    do jCHLD=jmin(8),my_CHLD
      do iCHLD=imax(8)+1,imax(8)+npts
        isf_draft_CHLD     (iCHLD,jCHLD) = isf_draft_PAR      ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) )
        Bathymetry_isf_CHLD(iCHLD,jCHLD) = Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) )
        Bathymetry_CHLD    (iCHLD,jCHLD) = Bathymetry_PAR     ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) )
      enddo
    enddo
    !- bdy_north(1) :
    do jCHLD=jmin(8)-npts, jmin(8)-1
      do iCHLD=imin(8),imax(8)
        isf_draft_CHLD     (iCHLD,jCHLD) = isf_draft_PAR      ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) )
        Bathymetry_isf_CHLD(iCHLD,jCHLD) = Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) )
        Bathymetry_CHLD    (iCHLD,jCHLD) = Bathymetry_PAR     ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) )
      enddo
    enddo

    !---- smooth transitions :
    !- bdy_west(1) :
    do jCHLD=jmin(8),my_CHLD
      do iCHLD=imax(8)+npts+1,imax(8)+2*npts
        isf_draft_CHLD     (iCHLD,jCHLD) = (   (imax(8)+2*npts-iCHLD+1) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                 +     (iCHLD-imax(8)-npts) * isf_draft_0 ( iCHLD, jCHLD ) ) / (npts+1)
        Bathymetry_isf_CHLD(iCHLD,jCHLD) = (   (imax(8)+2*npts-iCHLD+1) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                 +     (iCHLD-imax(8)-npts) * Bathymetry_isf_0 ( iCHLD, jCHLD ) ) / (npts+1)
        Bathymetry_CHLD    (iCHLD,jCHLD) = (   (imax(8)+2*npts-iCHLD+1) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                 +     (iCHLD-imax(8)-npts) * Bathymetry_0 ( iCHLD, jCHLD ) ) / (npts+1)
      enddo
    enddo
    !- bdy_north(1) :
    do jCHLD=jmin(8)-2*npts,jmin(8)-npts-1
      do iCHLD=imin(8),imax(8)
        isf_draft_CHLD     (iCHLD,jCHLD) = ( (jCHLD-jmin(8)+2*npts+1) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jCHLD) * isf_draft_0 ( iCHLD, jCHLD ) ) / (npts+1)                     
        Bathymetry_isf_CHLD(iCHLD,jCHLD) = ( (jCHLD-jmin(8)+2*npts+1) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jCHLD) * Bathymetry_isf_0 ( iCHLD, jCHLD ) ) / (npts+1)                     
        Bathymetry_CHLD    (iCHLD,jCHLD) = ( (jCHLD-jmin(8)+2*npts+1) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jCHLD) * Bathymetry_0 ( iCHLD, jCHLD ) ) / (npts+1)                     
      enddo
    enddo

    !- fill bellingshausen:
    do kk=1,Nbox
      isf_draft_CHLD      (imin(kk):imax(kk),jmin(kk):my_CHLD) = 0.0
      Bathymetry_isf_CHLD (imin(kk):imax(kk),jmin(kk):my_CHLD) = 0.0
      Bathymetry_CHLD     (imin(kk):imax(kk),jmin(kk):my_CHLD) = 0.0
    enddo

elseif ( TRIM(config) == 'AMUXL12' ) then
 
    write(*,*) 'Special correction for config ', TRIM(config)
 
    ! Note that you DO NOT have to change the following in you change the domain size
    ! through modifications of nn_imin_extract, nn_imax_extract, ... in the namelist 

    ! To keep the boxes at the same position if we slightly change the domain:
    i0 = 1771 - imin_EXT
    j0 =   30 - jmin_EXT
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
      imin(kk) = MIN( MAX( 1, imin(kk) ), mx_CHLD )
      imax(kk) = MIN( MAX( 1, imax(kk) ), mx_CHLD )
      jmin(kk) = MIN( MAX( 1, jmin(kk) ), my_CHLD )
      jmax(kk) = MIN( MAX( 1, jmax(kk) ), my_CHLD )
    enddo
    !-
    isf_draft_CHLD(imin(1):imax(1),jmin(1):jmax(1)) = Bathymetry_isf_CHLD(imin(1):imax(1),jmin(1):jmax(1))
    write(*,*) '   imin, imax = ', imin(:), imax(:)
    write(*,*) '   jmin, jmax = ', jmin(:), jmax(:)

    ! no isf over a safety zone (2*npts wide halo) from the eastern and western BDY :
    isf_draft_CHLD     (1:2*npts,:) = 0.0
    Bathymetry_isf_CHLD(1:2*npts,:) = Bathymetry_CHLD(1:2*npts,:)
    isf_draft_CHLD     (mx_CHLD-2*npts+1:mx_CHLD,:) = 0.0
    Bathymetry_isf_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:) = Bathymetry_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:)

    ! remove a narrow channel in the halo :
    Bathymetry_isf_CHLD(678:682,144:149)=0.e0
    isf_draft_CHLD(678:682,144:149)=0.e0
 
    ! remove ice shelves on islands and de-facto icebergs:
    isf_draft_CHLD(435:455,139:165)=0.e0
    isf_draft_CHLD(435:468,127:139)=0.e0
    isf_draft_CHLD(435:470,116:127)=0.e0
    isf_draft_CHLD(435:478,103:116)=0.e0

elseif ( TRIM(config(1:8)) == 'AMUXL025' ) then

    ! no isf over a safety zone (2*npts wide halo) from the eastern and western BDY :
    isf_draft_CHLD     (1:2*npts,:) = 0.0
    Bathymetry_isf_CHLD(1:2*npts,:) = Bathymetry_CHLD(1:2*npts,:)
    isf_draft_CHLD     (mx_CHLD-2*npts+1:mx_CHLD,:) = 0.0
    Bathymetry_isf_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:) = Bathymetry_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:)

    ! remove a narrow channel in the halo :
    Bathymetry_isf_CHLD(229,48:50)=0.e0
    isf_draft_CHLD(229,48:50)=0.e0

elseif ( TRIM(config) == 'CRODOT60' ) then

    ! no isf over a safety zone (2*npts wide halo) from the eastern and western BDY :
    isf_draft_CHLD     (1:2*npts,:) = 0.0
    Bathymetry_isf_CHLD(1:2*npts,:) = Bathymetry_CHLD(1:2*npts,:)
    isf_draft_CHLD     (mx_CHLD-2*npts+1:mx_CHLD,:) = 0.0
    Bathymetry_isf_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:) = Bathymetry_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:)

    ! put ice in the shear margin of Crosson:
    do iCHLD=260,271
    do jCHLD=158,202
      isf_draft_CHLD(iCHLD,jCHLD) = MAX(10.5,isf_draft_CHLD(iCHLD,jCHLD))
    enddo
    enddo

    ! prevent weird narrow cavities:
    isf_draft_CHLD(263:275,221:230)=Bathymetry_isf_CHLD(263:275,221:230)-Bathymetry_CHLD(263:275,221:230)

elseif ( TRIM(config) == 'eAMUXL12.L121' ) then

    write(*,*) 'Special correction for config ', TRIM(config)

    ! To keep the boxes at the same position if we slightly change the domain:
    i0 = 1737 - imin_EXT
    j0 =  562 - jmin_EXT
    write(*,*) '   i0 = ', i0
    write(*,*) '   j0 = ', j0

    ! to avoid hole in an ice shelf:
    isf_draft_CHLD(557+i0,199+j0:200+j0) = Bathymetry_isf_CHLD(557+i0,199+j0:200+j0)

    ! no isf over a safety zone (2*npts wide halo) from the eastern and western BDY :
    isf_draft_CHLD     (1:2*npts,:) = 0.0
    Bathymetry_isf_CHLD(1:2*npts,:) = Bathymetry_CHLD(1:2*npts,:)
    isf_draft_CHLD     (mx_CHLD-2*npts+1:mx_CHLD,:) = 0.0
    Bathymetry_isf_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:) = Bathymetry_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:)

endif

!=================================================================================
! 5- Writing new CHLD bathymetry file :
!=================================================================================

write(*,*) 'Creating ', TRIM(file_bathy_out)

status = NF90_CREATE(TRIM(file_bathy_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')                     

status = NF90_DEF_DIM(fidM,"y",my_CHLD,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx_CHLD,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")

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
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imin_extraction",imin_EXT)  ; call erreur(status,.TRUE.,"put_att_GLOBAL3")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imax_extraction",imax_EXT)  ; call erreur(status,.TRUE.,"put_att_GLOBAL4")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmin_extraction",jmin_EXT)  ; call erreur(status,.TRUE.,"put_att_GLOBAL5")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmax_extraction",jmax_EXT)  ; call erreur(status,.TRUE.,"put_att_GLOBAL6")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidM,nav_lat_ID,gphit_CHLD); call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,glamt_CHLD); call erreur(status,.TRUE.,"var_nav_lon_ID")
if ( ln_isfcav ) then
  status = NF90_PUT_VAR(fidM,isf_draft_ID,isf_draft_CHLD);           call erreur(status,.TRUE.,"var_isf_draft_ID")
  status = NF90_PUT_VAR(fidM,Bathymetry_isf_ID,Bathymetry_isf_CHLD); call erreur(status,.TRUE.,"var_Bathymetry_isf_ID")
endif
status = NF90_PUT_VAR(fidM,Bathymetry_ID,Bathymetry_CHLD); call erreur(status,.TRUE.,"var_Bathymetry_ID")

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
