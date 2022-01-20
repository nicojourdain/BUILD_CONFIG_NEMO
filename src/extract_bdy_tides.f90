program modif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Used to build netcdf tide lateral boundary conditions for the child domain (CHLD)
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read input file dimensions in first existing file for specified time window
! 3- Process all gridT files over specified period
!
! history : - Mar. 2015: initial version (N. Jourdain, CNRS-LGGE)
!           - Mar. 2017: version with namelist (N. Jourdain, CNRS-IGE)
!           - Jan. 2022: a bit of cleaning (N. Jourdain)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf

IMPLICIT NONE

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /bdy/ nn_bdy_east, nn_bdy_west, nn_bdy_north, nn_bdy_south, nn_harm
namelist /bdy_tide/ tide_dir, harm

CHARACTER(LEN=50)                         :: config
CHARACTER(LEN=150)                        :: config_dir, tide_dir
INTEGER                                   :: nn_bdy_east, nn_bdy_west, nn_bdy_north, nn_bdy_south, nn_harm
CHARACTER(LEN=4),ALLOCATABLE,DIMENSION(:) :: harm

INTEGER :: fidELEV, fidUV, status, dimID_y, dimID_x, my, mx, lat_ID, lon_ID, Hg_ID, Ha_ID, Vg_ID, Va_ID, Ug_ID, Ua_ID, knb,    & 
&          fidCOORD, dimID_yb, dimID_xbt, dimID_xbu, dimID_xbv, myb, mxbt, mxbu, mxbv, glamt_ID, gphit_ID, glamu_ID, gphiu_ID, &
&          glamv_ID, gphiv_ID, nbit_ID, nbjt_ID, nbrt_ID, nbiu_ID, nbju_ID, nbru_ID, nbiv_ID, nbjv_ID, nbrv_ID, z1_ID, z2_ID,  &
&          u1_ID, u2_ID, v1_ID, v2_ID, fidG, dimIDG_x, dimIDG_y, dimIDG_z, dimIDG_t, mxG, myG, mzG, mtG, tmask_ID, umask_ID,   &
&          vmask_ID, kbdy, iinf, isup, jinf, jsup, itmp, jtmp, kharm, fidoutT, fidoutU, fidoutV, i, j

CHARACTER(LEN=150) :: file_in_elev, file_in_uv, file_coord, file_mesh_mask, file_out_T, file_out_U, file_out_V, command_str

CHARACTER(LEN=4) :: harmstr

REAL*8,ALLOCATABLE,DIMENSION(:) :: lat, lon, zlon

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: Hg, Ha,  Vg, Va, Ug, Ua, z1, z2, u1, u2, v1, v2, glamt, gphit, glamu, gphiu, glamv, gphiv,&
&                                    zglamt, zglamu, zglamv, gphiu_CHLD, glamu_CHLD, gphiv_CHLD, glamv_CHLD, zglamu_CHLD, zglamv_CHLD

REAL*8 :: Hg_bdy, Ha_bdy, Vg_bdy, Va_bdy, Ug_bdy, Ua_bdy, div, zrad, angle, tmp1, tmp2, Hcs, Hsn, Ucs, Usn, Vcs, Vsn

INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbit, nbjt, nbrt, nbiu, nbju, nbru, nbiv, nbjv, nbrv

INTEGER*1,ALLOCATABLE,DIMENSION(:,:) :: mskfes

INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:,:) :: tmask, umask, vmask, fmask

LOGICAL :: idateline !! to put in the namelist at some point

!-----------------------------------------------------------------

! Default values (replaced with namelist values if specified):
config_dir        = '.'
tide_dir          = '.'
nn_harm           =  0

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=bdy)
CLOSE(1)

if ( nn_harm .lt. 1 ) then
  write(*,*) 'nn_harm < 1   NO TIDAL HARMONIC TO PRCESS'
  stop
else
  write(*,*) 'number of tidal harmonics: ', nn_harm
endif
ALLOCATE( harm(nn_harm) )

!- read remaining namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=bdy_tide)
CLOSE(1)

!- bdy coordinates :
write(file_coord,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/coordinates_bdy_',a,'.nc')

!- mesh/mask of child domain :
write(file_mesh_mask,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/mesh_mask_',a,'.nc')

!- format for file_in_uv and file_in_elev :
301 FORMAT(a,'/',a,'.FES2012.elev.nc')
302 FORMAT(a,'/',a,'.FES2012.uv.nc')

!- Structure of output files :
401 FORMAT(a,'/BDY/bdytide_',a,'_',a,'_grid_T.nc')
402 FORMAT(a,'/BDY/bdytide_',a,'_',a,'_grid_U.nc')
403 FORMAT(a,'/BDY/bdytide_',a,'_',a,'_grid_V.nc')

zrad = 3.14159265358979323846264338327 / 180.000000000000000000000000000

idateline = .false. !! = .true. if dateline withing the child domain

!--
write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir ',a,'/BDY')
CALL system(TRIM(command_str))

!---------------------------------------
! Read  bdy coordinates :

write(*,*) 'Reading ', TRIM(file_coord)

status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read bdy coord")

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbt",dimID_xbt) ; call erreur(status,.TRUE.,"inq_dimID_xbt")
status = NF90_INQ_DIMID(fidCOORD,"xbu",dimID_xbu) ; call erreur(status,.TRUE.,"inq_dimID_xbu")
status = NF90_INQ_DIMID(fidCOORD,"xbv",dimID_xbv) ; call erreur(status,.TRUE.,"inq_dimID_xbv")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbt,len=mxbt) ; call erreur(status,.TRUE.,"inq_dim_xbt")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbu,len=mxbu) ; call erreur(status,.TRUE.,"inq_dim_xbu")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbv,len=mxbv) ; call erreur(status,.TRUE.,"inq_dim_xbv")

ALLOCATE(  glamt(mxbt,myb), zglamt(mxbt,myb)  )
ALLOCATE(  gphit(mxbt,myb), zglamu(mxbt,myb)  )
ALLOCATE(  glamu(mxbu,myb), zglamv(mxbt,myb)  )
ALLOCATE(  gphiu(mxbu,myb) )
ALLOCATE(  glamv(mxbv,myb) )
ALLOCATE(  gphiv(mxbv,myb) )
ALLOCATE(  nbit(mxbt,myb)  )
ALLOCATE(  nbjt(mxbt,myb)  )
ALLOCATE(  nbrt(mxbt,myb)  )
ALLOCATE(  nbiu(mxbu,myb)  )
ALLOCATE(  nbju(mxbu,myb)  )
ALLOCATE(  nbru(mxbu,myb)  )
ALLOCATE(  nbiv(mxbv,myb)  )
ALLOCATE(  nbjv(mxbv,myb)  )
ALLOCATE(  nbrv(mxbv,myb)  )

status = NF90_INQ_VARID(fidCOORD,"glamt",glamt_ID) ;call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidCOORD,"gphit",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidCOORD,"glamu",glamu_ID) ; call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidCOORD,"gphiu",gphiu_ID) ; call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidCOORD,"glamv",glamv_ID) ; call erreur(status,.TRUE.,"inq_glamv_ID")
status = NF90_INQ_VARID(fidCOORD,"gphiv",gphiv_ID) ; call erreur(status,.TRUE.,"inq_gphiv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbit",nbit_ID)   ; call erreur(status,.TRUE.,"inq_nbit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjt",nbjt_ID)   ; call erreur(status,.TRUE.,"inq_nbjt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrt",nbrt_ID)   ; call erreur(status,.TRUE.,"inq_nbrt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiu",nbiu_ID)   ; call erreur(status,.TRUE.,"inq_nbiu_ID")
status = NF90_INQ_VARID(fidCOORD,"nbju",nbju_ID)   ; call erreur(status,.TRUE.,"inq_nbju_ID")
status = NF90_INQ_VARID(fidCOORD,"nbru",nbru_ID)   ; call erreur(status,.TRUE.,"inq_nbru_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiv",nbiv_ID)   ; call erreur(status,.TRUE.,"inq_nbiv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjv",nbjv_ID)   ; call erreur(status,.TRUE.,"inq_nbjv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrv",nbrv_ID)   ; call erreur(status,.TRUE.,"inq_nbrv_ID")

status = NF90_GET_VAR(fidCOORD,glamt_ID,glamt) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidCOORD,gphit_ID,gphit) ; call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidCOORD,glamu_ID,glamu) ; call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidCOORD,gphiu_ID,gphiu) ; call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidCOORD,glamv_ID,glamv) ; call erreur(status,.TRUE.,"getvar_glamv")
status = NF90_GET_VAR(fidCOORD,gphiv_ID,gphiv) ; call erreur(status,.TRUE.,"getvar_gphiv")
status = NF90_GET_VAR(fidCOORD,nbit_ID,nbit)   ; call erreur(status,.TRUE.,"getvar_nbit")
status = NF90_GET_VAR(fidCOORD,nbjt_ID,nbjt)   ; call erreur(status,.TRUE.,"getvar_nbjt")
status = NF90_GET_VAR(fidCOORD,nbrt_ID,nbrt)   ; call erreur(status,.TRUE.,"getvar_nbrt")
status = NF90_GET_VAR(fidCOORD,nbiu_ID,nbiu)   ; call erreur(status,.TRUE.,"getvar_nbiu")
status = NF90_GET_VAR(fidCOORD,nbju_ID,nbju)   ; call erreur(status,.TRUE.,"getvar_nbju")
status = NF90_GET_VAR(fidCOORD,nbru_ID,nbru)   ; call erreur(status,.TRUE.,"getvar_nbru")
status = NF90_GET_VAR(fidCOORD,nbiv_ID,nbiv)   ; call erreur(status,.TRUE.,"getvar_nbiv")
status = NF90_GET_VAR(fidCOORD,nbjv_ID,nbjv)   ; call erreur(status,.TRUE.,"getvar_nbjv")
status = NF90_GET_VAR(fidCOORD,nbrv_ID,nbrv)   ; call erreur(status,.TRUE.,"getvar_nbrv")

status = NF90_CLOSE(fidCOORD) ; call erreur(status,.TRUE.,"fin_lecture")

!---------------------------------------                   
! Read mesh_mask of child domain
                             
write(*,*) 'Reading ', TRIM(file_mesh_mask)
                              
status = NF90_OPEN(TRIM(file_mesh_mask),0,fidG) ; call erreur(status,.TRUE.,"read child mesh/mask") 
                                                           
status = NF90_INQ_DIMID(fidG,"x",dimIDG_x) ; call erreur(status,.TRUE.,"inq_dimIDG_x")
status = NF90_INQ_DIMID(fidG,"y",dimIDG_y) ; call erreur(status,.TRUE.,"inq_dimIDG_y")
status = NF90_INQ_DIMID(fidG,"z",dimIDG_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidG,"nav_lev",dimIDG_z)
call erreur(status,.TRUE.,"inq_dimIDG_z")
status = NF90_INQ_DIMID(fidG,"t",dimIDG_t)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidG,"time_counter",dimIDG_t)
call erreur(status,.TRUE.,"inq_dimIDG_t")
                                                               
status = NF90_INQUIRE_DIMENSION(fidG,dimIDG_x,len=mxG) ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidG,dimIDG_y,len=myG) ; call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidG,dimIDG_z,len=mzG) ; call erreur(status,.TRUE.,"inq_dim_z")
status = NF90_INQUIRE_DIMENSION(fidG,dimIDG_t,len=mtG) ; call erreur(status,.TRUE.,"inq_dim_t")
                               
ALLOCATE(  tmask(mxG,myG,mzG,mtG)  ) 
ALLOCATE(  umask(mxG,myG,mzG,mtG)  ) 
ALLOCATE(  vmask(mxG,myG,mzG,mtG)  )
ALLOCATE(  gphiu_CHLD(mxG,myG), glamu_CHLD(mxG,myG), zglamu_CHLD(mxG,myG) ) 
ALLOCATE(  gphiv_CHLD(mxG,myG), glamv_CHLD(mxG,myG), zglamv_CHLD(mxG,myG) ) 
                                 
status = NF90_INQ_VARID(fidG,"tmask",tmask_ID) ; call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_INQ_VARID(fidG,"umask",umask_ID) ; call erreur(status,.TRUE.,"inq_umask_ID")
status = NF90_INQ_VARID(fidG,"vmask",vmask_ID) ; call erreur(status,.TRUE.,"inq_vmask_ID")
status = NF90_INQ_VARID(fidG,"gphiu",gphiu_ID) ; call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidG,"glamu",glamu_ID) ; call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidG,"gphiv",gphiv_ID) ; call erreur(status,.TRUE.,"inq_gphiv_ID")
status = NF90_INQ_VARID(fidG,"glamv",glamv_ID) ; call erreur(status,.TRUE.,"inq_glamv_ID")
                                                              
status = NF90_GET_VAR(fidG,tmask_ID,tmask)     ; call erreur(status,.TRUE.,"getvar_tmask")
status = NF90_GET_VAR(fidG,umask_ID,umask)     ; call erreur(status,.TRUE.,"getvar_umask")
status = NF90_GET_VAR(fidG,vmask_ID,vmask)     ; call erreur(status,.TRUE.,"getvar_vmask")
status = NF90_GET_VAR(fidG,gphiu_ID,gphiu_CHLD) ; call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidG,glamu_ID,glamu_CHLD) ; call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidG,gphiv_ID,gphiv_CHLD) ; call erreur(status,.TRUE.,"getvar_gphiv")
status = NF90_GET_VAR(fidG,glamv_ID,glamv_CHLD) ; call erreur(status,.TRUE.,"getvar_glamv")
                                                      
status = NF90_CLOSE(fidG) ; call erreur(status,.TRUE.,"fin_lecture")     

if ( idateline) then

  write(*,*) 'swithc CHLD longitudes to [ 0; 360 ]'

  where( glamu_CHLD(:,:) .lt. 0.0 )
    zglamu_CHLD(:,:) = 360.0 + glamu_CHLD(:,:)
  elsewhere
    zglamu_CHLD(:,:) = glamu_CHLD(:,:)
  endwhere

  where( glamv_CHLD(:,:) .lt. 0.0 )
    zglamv_CHLD(:,:) = 360.0 + glamv_CHLD(:,:)
  elsewhere
    zglamv_CHLD(:,:) = glamv_CHLD(:,:)
  endwhere

else

  write(*,*) 'keep CHLD longitudes in [ -180; 180]'
  zglamu_CHLD(:,:) = glamu_CHLD(:,:)
  zglamv_CHLD(:,:) = glamv_CHLD(:,:)

endif

!----------------------------------------------------------------------------------------

if ( idateline ) then

  write(*,*) 'switch bdy longitudes to [ 0; 360 ]'
  
  where( glamt(:,:) .lt. 0.0 )
    zglamt(:,:) = 360.0 + glamt(:,:)
  elsewhere
    zglamt(:,:) = glamt(:,:)
  endwhere
  
  where( glamu(:,:) .lt. 0.0 )
    zglamu(:,:) = 360.0 + glamu(:,:)
  elsewhere
    zglamu(:,:) = glamu(:,:)
  endwhere
  
  where( glamv(:,:) .lt. 0.0 )
    zglamv(:,:) = 360.0 + glamv(:,:)
  elsewhere
    zglamv(:,:) = glamv(:,:)
  endwhere

else

  write(*,*) 'keep bdy longitudes in [ -180; 180 ]'

  zglamt(:,:) = glamt(:,:)
  zglamu(:,:) = glamu(:,:)
  zglamv(:,:) = glamv(:,:)

endif

!-----------------------------------------------------------------

DO kharm=1,nn_harm

   if ( TRIM(harm(kharm)) == 'Mu2' ) then
     harmstr = 'MU2 '
   elseif ( TRIM(harm(kharm)) == 'Nu2' ) then
     harmstr = 'NU2 '
   else
     harmstr = harm(kharm)
   endif

   !-----------------------------------------------------------------
   ! Read tides elevation :
  
   write(file_in_elev,301) TRIM(tide_dir), TRIM(harm(kharm))
   
   write(*,*) 'Reading ', TRIM(file_in_elev)
   
   status = NF90_OPEN(TRIM(file_in_elev),0,fidELEV) ; call erreur(status,.TRUE.,"read tide elevation") 
   
   status = NF90_INQ_DIMID(fidELEV,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
   status = NF90_INQ_DIMID(fidELEV,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
                                                                  
   status = NF90_INQUIRE_DIMENSION(fidELEV,dimID_y,len=my) ; call erreur(status,.TRUE.,"inq_dim_y")
   status = NF90_INQUIRE_DIMENSION(fidELEV,dimID_x,len=mx) ; call erreur(status,.TRUE.,"inq_dim_x")
                          
   ALLOCATE(  lat(my)     ) 
   ALLOCATE(  lon(mx), zlon(mx)  ) 
   ALLOCATE(  Hg (mx,my)  ) 
   ALLOCATE(  Ha (mx,my)  ) 
                                    
   status = NF90_INQ_VARID(fidELEV,"lat",lat_ID) ; call erreur(status,.TRUE.,"inq_lat_ID")
   status = NF90_INQ_VARID(fidELEV,"lon",lon_ID) ; call erreur(status,.TRUE.,"inq_lon_ID")
   status = NF90_INQ_VARID(fidELEV,"Hg",Hg_ID)   ; call erreur(status,.TRUE.,"inq_Hg_ID")
   status = NF90_INQ_VARID(fidELEV,"Ha",Ha_ID)   ; call erreur(status,.TRUE.,"inq_Ha_ID")
                                                                 
   status = NF90_GET_VAR(fidELEV,lat_ID,lat) ; call erreur(status,.TRUE.,"getvar_lat")
   status = NF90_GET_VAR(fidELEV,lon_ID,lon) ; call erreur(status,.TRUE.,"getvar_lon")
   status = NF90_GET_VAR(fidELEV,Hg_ID,Hg)   ; call erreur(status,.TRUE.,"getvar_Hg")
   status = NF90_GET_VAR(fidELEV,Ha_ID,Ha)   ; call erreur(status,.TRUE.,"getvar_Ha")
                                                         
   status = NF90_CLOSE(fidELEV) ; call erreur(status,.TRUE.,"fin_lecture")     
  
   if ( idateline ) then
     write(*,*) 'keep FES2012 longitude in [ 0; 360 ]'
     zlon(:) = lon(:)
   else
     write(*,*) 'switch FES2012 longitude to [ -180; 180]'
     where( lon(:) .gt. 180.00000000 )
       zlon(:) = lon(:) - 360.00000000
     elsewhere
       zlon(:) = lon(:)
     endwhere
   endif

   write(*,*) 'max(Ha) = ', maxval(maxval(abs(Ha),2),1)
   write(*,*) 'max(Hg) = ', maxval(maxval(abs(Hg),2),1)

   !---------------------------------------                   
   ! Read tidal u,v
   
   write(file_in_uv,302) TRIM(tide_dir), TRIM(harm(kharm))
   
   write(*,*) 'Reading ', TRIM(file_in_uv)
   
   status = NF90_OPEN(TRIM(file_in_uv),0,fidUV) ; call erreur(status,.TRUE.,"read tidal U,V") 
   
   ALLOCATE(  Vg(mx,my)  ) 
   ALLOCATE(  Va(mx,my)  ) 
   ALLOCATE(  Ug(mx,my)  ) 
   ALLOCATE(  Ua(mx,my)  ) 
   
   status = NF90_INQ_VARID(fidUV,"Vg",Vg_ID) ; call erreur(status,.TRUE.,"inq_Vg_ID")
   status = NF90_INQ_VARID(fidUV,"Va",Va_ID) ; call erreur(status,.TRUE.,"inq_Va_ID")
   status = NF90_INQ_VARID(fidUV,"Ug",Ug_ID) ; call erreur(status,.TRUE.,"inq_Ug_ID")
   status = NF90_INQ_VARID(fidUV,"Ua",Ua_ID) ; call erreur(status,.TRUE.,"inq_Ua_ID")
   
   status = NF90_GET_VAR(fidUV,Vg_ID,Vg) ; call erreur(status,.TRUE.,"getvar_Vg")
   status = NF90_GET_VAR(fidUV,Va_ID,Va) ; call erreur(status,.TRUE.,"getvar_Va")
   status = NF90_GET_VAR(fidUV,Ug_ID,Ug) ; call erreur(status,.TRUE.,"getvar_Ug")
   status = NF90_GET_VAR(fidUV,Ua_ID,Ua) ; call erreur(status,.TRUE.,"getvar_Ua")
   
   status = NF90_CLOSE(fidUV) ; call erreur(status,.TRUE.,"fin_lecture")     

   write(*,*) 'max(Va) = ', maxval(maxval(abs(Va),2),1)
   write(*,*) 'max(Vg) = ', maxval(maxval(abs(Vg),2),1)
   write(*,*) 'max(Ua) = ', maxval(maxval(abs(Ua),2),1)
   write(*,*) 'max(Ug) = ', maxval(maxval(abs(Ug),2),1)   

   !----------------------------------------------------------------------------------------
   ! Land/Sea mask for FES2012 :
  
   write(*,*) 'Define land/sea mask for FES2012'
 
   ALLOCATE( mskfes(mx,my) )
  
   do i=1,mx
   do j=1,my
     if (       abs(Ha(i,j)) .lt. 1.e5     &
     &    .and. abs(Hg(i,j)) .lt. 1.e5     &
     &    .and. abs(Ua(i,j)) .lt. 1.e5     &
     &    .and. abs(Ug(i,j)) .lt. 1.e5     &
     &    .and. abs(Va(i,j)) .lt. 1.e5     &
     &    .and. abs(Vg(i,j)) .lt. 1.e5     ) then
       mskfes(i,j)=1
     else
       mskfes(i,j)=0 
     endif
   enddo
   enddo

   !----------------------------------------------------------------------------------------
   
   write(*,*) 'Interpolation of z1 and z2'
   
   ALLOCATE( z1(mxbt,myb), z2(mxbt,myb) )
   
   do kbdy=1,mxbt
  
     write(*,*) 'nj1: ', kbdy, tmask(nbit(kbdy,1),nbjt(kbdy,1),1,1) 
 
     if ( tmask(nbit(kbdy,1),nbjt(kbdy,1),1,1) .ne. 0 ) then
 
       itmp=MINLOC(abs(zlon-zglamt(kbdy,1)),1)
       jtmp=MINLOC(abs(lat-gphit(kbdy,1)),1)
 
       do knb=1,10 ! looking for closest neighbour (knb=1), then expand if not found
  
         if ( zlon(itmp) .lt. zglamt(kbdy,1) ) then
           iinf=itmp-knb+1
           isup=itmp+knb
         else
           iinf=itmp-knb
           isup=itmp+knb-1
         endif
       
         if ( lat(jtmp) .lt. gphit(kbdy,1) ) then
           jinf=jtmp-knb+1
           jsup=jtmp+knb
         else
           jinf=jtmp-knb
           jsup=jtmp+knb-1
         endif
       
         Ha_bdy =   ( zlon(isup) - zglamt(kbdy,1) ) * ( lat(jsup) - gphit(kbdy,1) ) * Ha(iinf,jinf) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamt(kbdy,1) ) * ( gphit(kbdy,1) - lat(jinf) ) * Ha(iinf,jsup) * mskfes(iinf,jsup)  &
         &        + ( zglamt(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphit(kbdy,1) ) * Ha(isup,jinf) * mskfes(isup,jinf)  &
         &        + ( zglamt(kbdy,1) - zlon(iinf) ) * ( gphit(kbdy,1) - lat(jinf) ) * Ha(isup,jsup) * mskfes(isup,jsup)
       
         Hcs =      ( zlon(isup) - zglamt(kbdy,1) ) * ( lat(jsup) - gphit(kbdy,1) ) * COS(zrad*Hg(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamt(kbdy,1) ) * ( gphit(kbdy,1) - lat(jinf) ) * COS(zrad*Hg(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamt(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphit(kbdy,1) ) * COS(zrad*Hg(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamt(kbdy,1) - zlon(iinf) ) * ( gphit(kbdy,1) - lat(jinf) ) * COS(zrad*Hg(isup,jsup)) * mskfes(isup,jsup)
      
         Hsn =      ( zlon(isup) - zglamt(kbdy,1) ) * ( lat(jsup) - gphit(kbdy,1) ) * SIN(zrad*Hg(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamt(kbdy,1) ) * ( gphit(kbdy,1) - lat(jinf) ) * SIN(zrad*Hg(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamt(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphit(kbdy,1) ) * SIN(zrad*Hg(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamt(kbdy,1) - zlon(iinf) ) * ( gphit(kbdy,1) - lat(jinf) ) * SIN(zrad*Hg(isup,jsup)) * mskfes(isup,jsup)
  
         div =   ( zlon(isup) - zglamt(kbdy,1) ) * ( lat(jsup) - gphit(kbdy,1) ) * mskfes(iinf,jinf)  &
         &     + ( zlon(isup) - zglamt(kbdy,1) ) * ( gphit(kbdy,1) - lat(jinf) ) * mskfes(iinf,jsup)  &
         &     + ( zglamt(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphit(kbdy,1) ) * mskfes(isup,jinf)  &
         &     + ( zglamt(kbdy,1) - zlon(iinf) ) * ( gphit(kbdy,1) - lat(jinf) ) * mskfes(isup,jsup)
      
         write(*,*) '   nj2: ', knb, Ha_bdy, Hcs, Hsn, div
 
         if ( abs(div) .gt. 1.e-12 ) then
           Ha_bdy = Ha_bdy / div
           Hg_bdy = ATAN2( Hsn/div, Hcs/div ) / zrad !- to have a proper "modulo" interpolation of the phase
           exit
         endif
 
       enddo
 
       if ( .not. abs(div) .gt. 1.e-12 ) then
         write(*,*) '~!@#$%^ INCREASE NB OF LOOPS ON knb (bdyT)  >>>>>>>>>>  stop !!'
         write(*,*) kbdy, iinf, jinf, Ha_bdy, Hcs, Hsn
         stop
       endif

     else
    
        Ha_bdy = 0.0
        Hg_bdy = 0.0
   
     endif !--  if ( tmask(nbit(kbdy,1),nbjt(kbdy,1),1) .ne. 0 )
   
     z1(kbdy,1) = Ha_bdy * COS( zrad * Hg_bdy )
     z2(kbdy,1) = Ha_bdy * SIN( zrad * Hg_bdy )
   
   enddo
   
   !----------------------------------------------------------------------------------------
   
   write(*,*) 'Interpolation of u1 and u2 (from Utide,Vtide brought at U pt then rotated)'
   
   ALLOCATE( u1(mxbu,myb), u2(mxbu,myb) )
   
   do kbdy=1,mxbu
   
     if ( umask(nbiu(kbdy,1),nbju(kbdy,1),1,1) .ne. 0 ) then
   
       itmp=MINLOC(abs(zlon-zglamu(kbdy,1)),1)
       jtmp=MINLOC(abs(lat-gphiu(kbdy,1)),1)
 
       do knb=1,10 ! looking for closest neighbour (knb=1), then expand if not found
  
         if ( zlon(itmp) .lt. zglamu(kbdy,1) ) then
           iinf=itmp-knb+1
           isup=itmp+knb
         else
           iinf=itmp-knb
           isup=itmp+knb-1
         endif
       
         if ( lat(jtmp) .lt. gphiu(kbdy,1) ) then
           jinf=jtmp-knb+1
           jsup=jtmp+knb
         else
           jinf=jtmp-knb
           jsup=jtmp+knb-1
         endif

         Ua_bdy =   ( zlon(isup) - zglamu(kbdy,1) ) * ( lat(jsup) - gphiu(kbdy,1) ) * Ua(iinf,jinf) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamu(kbdy,1) ) * ( gphiu(kbdy,1) - lat(jinf) ) * Ua(iinf,jsup) * mskfes(iinf,jsup)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiu(kbdy,1) ) * Ua(isup,jinf) * mskfes(isup,jinf)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( gphiu(kbdy,1) - lat(jinf) ) * Ua(isup,jsup) * mskfes(isup,jsup)
       
         Ucs =      ( zlon(isup) - zglamu(kbdy,1) ) * ( lat(jsup) - gphiu(kbdy,1) ) * COS(zrad*Ug(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamu(kbdy,1) ) * ( gphiu(kbdy,1) - lat(jinf) ) * COS(zrad*Ug(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiu(kbdy,1) ) * COS(zrad*Ug(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( gphiu(kbdy,1) - lat(jinf) ) * COS(zrad*Ug(isup,jsup)) * mskfes(isup,jsup)
     
         Usn =      ( zlon(isup) - zglamu(kbdy,1) ) * ( lat(jsup) - gphiu(kbdy,1) ) * SIN(zrad*Ug(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamu(kbdy,1) ) * ( gphiu(kbdy,1) - lat(jinf) ) * SIN(zrad*Ug(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiu(kbdy,1) ) * SIN(zrad*Ug(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( gphiu(kbdy,1) - lat(jinf) ) * SIN(zrad*Ug(isup,jsup)) * mskfes(isup,jsup)
  
         Va_bdy =   ( zlon(isup) - zglamu(kbdy,1) ) * ( lat(jsup) - gphiu(kbdy,1) ) * Va(iinf,jinf) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamu(kbdy,1) ) * ( gphiu(kbdy,1) - lat(jinf) ) * Va(iinf,jsup) * mskfes(iinf,jsup)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiu(kbdy,1) ) * Va(isup,jinf) * mskfes(isup,jinf)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( gphiu(kbdy,1) - lat(jinf) ) * Va(isup,jsup) * mskfes(isup,jsup)
       
         Vcs =      ( zlon(isup) - zglamu(kbdy,1) ) * ( lat(jsup) - gphiu(kbdy,1) ) * COS(zrad*Vg(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamu(kbdy,1) ) * ( gphiu(kbdy,1) - lat(jinf) ) * COS(zrad*Vg(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiu(kbdy,1) ) * COS(zrad*Vg(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( gphiu(kbdy,1) - lat(jinf) ) * COS(zrad*Vg(isup,jsup)) * mskfes(isup,jsup)
 
         Vsn =      ( zlon(isup) - zglamu(kbdy,1) ) * ( lat(jsup) - gphiu(kbdy,1) ) * SIN(zrad*Vg(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamu(kbdy,1) ) * ( gphiu(kbdy,1) - lat(jinf) ) * SIN(zrad*Vg(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiu(kbdy,1) ) * SIN(zrad*Vg(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamu(kbdy,1) - zlon(iinf) ) * ( gphiu(kbdy,1) - lat(jinf) ) * SIN(zrad*Vg(isup,jsup)) * mskfes(isup,jsup)
 
         tmp1  =    gphiu_CHLD( MIN(mxG,nbiu(kbdy,1)+1), nbju(kbdy,1) ) -  gphiu_CHLD( MAX(1,nbiu(kbdy,1)-1), nbju(kbdy,1) )
         tmp2  = ( zglamu_CHLD( MIN(mxG,nbiu(kbdy,1)+1), nbju(kbdy,1) ) - zglamu_CHLD( MAX(1,nbiu(kbdy,1)-1), nbju(kbdy,1) ) ) &
         &       * cos( zrad * gphiu_CHLD(nbiu(kbdy,1),nbju(kbdy,1)) )
         angle = ATAN2( tmp1, tmp2 )
 
         div =   ( zlon(isup) - zglamu(kbdy,1) ) * ( lat(jsup) - gphiu(kbdy,1) ) * mskfes(iinf,jinf)  &
           &   + ( zlon(isup) - zglamu(kbdy,1) ) * ( gphiu(kbdy,1) - lat(jinf) ) * mskfes(iinf,jsup)  &
           &   + ( zglamu(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiu(kbdy,1) ) * mskfes(isup,jinf)  &
           &   + ( zglamu(kbdy,1) - zlon(iinf) ) * ( gphiu(kbdy,1) - lat(jinf) ) * mskfes(isup,jsup)
       
         if ( abs(div) .gt. 1.e-12 ) then
           Ua_bdy = Ua_bdy / div
           Ug_bdy = ATAN2( Usn/div, Ucs/div ) / zrad !- to have a proper "modulo" interpolation of the phase
           Va_bdy = Va_bdy / div
           Vg_bdy = ATAN2( Vsn/div, Vcs/div ) / zrad !- to have a proper "modulo" interpolation of the phase
           exit
         endif

       enddo  !- knb

       if ( .not. abs(div) .gt. 1.e-12 ) then
         write(*,*) '@@@@@@@ INCREASE NB OF LOOPS ON knb (bdyU) >>>>>>>>>>  stop !!'
         write(*,*) kbdy, iinf, jinf, Ua_bdy, Ug_bdy
         stop
       endif
       
     else
    
        Ua_bdy = 0.0
        Ug_bdy = 0.0
        Va_bdy = 0.0
        Vg_bdy = 0.0
        angle  = 0.0   
 
     endif !--  if ( tmask(nbit(kbdy,1),nbjt(kbdy,1),1) .ne. 0 )
  
     u1(kbdy,1) = Ua_bdy * COS( zrad * Ug_bdy ) * COS(angle) + Va_bdy * COS( zrad * Vg_bdy ) * SIN(angle)  ! real part
     u2(kbdy,1) = Ua_bdy * SIN( zrad * Ug_bdy ) * COS(angle) + Va_bdy * SIN( zrad * Vg_bdy ) * SIN(angle)  ! imaginary part
  
   enddo
   
   !----------------------------------------------------------------------------------------
   
   write(*,*) 'Interpolation of v1 and v2 (from Utide,Vtide brought at V pt then rotated)'
   
   ALLOCATE( v1(mxbv,myb), v2(mxbv,myb) )
   
   do kbdy=1,mxbv
   
     if ( vmask(nbiv(kbdy,1),nbjv(kbdy,1),1,1) .ne. 0 ) then
  
       itmp=MINLOC(abs(zlon-zglamv(kbdy,1)),1)
       jtmp=MINLOC(abs(lat-gphiv(kbdy,1)),1)
 
       do knb=1,10 ! looking for closest neighbour (knb=1), then expand if not found
  
         if ( zlon(itmp) .lt. zglamv(kbdy,1) ) then
           iinf=itmp-knb+1
           isup=itmp+knb
         else
           iinf=itmp-knb
           isup=itmp+knb-1
         endif
       
         if ( lat(jtmp) .lt. gphiv(kbdy,1) ) then
           jinf=jtmp-knb+1
           jsup=jtmp+knb
         else
           jinf=jtmp-knb
           jsup=jtmp+knb-1
         endif

         Ua_bdy =   ( zlon(isup) - zglamv(kbdy,1) ) * ( lat(jsup) - gphiv(kbdy,1) ) * Ua(iinf,jinf) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamv(kbdy,1) ) * ( gphiv(kbdy,1) - lat(jinf) ) * Ua(iinf,jsup) * mskfes(iinf,jsup)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiv(kbdy,1) ) * Ua(isup,jinf) * mskfes(isup,jinf)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( gphiv(kbdy,1) - lat(jinf) ) * Ua(isup,jsup) * mskfes(isup,jsup)
       
         Ucs =      ( zlon(isup) - zglamv(kbdy,1) ) * ( lat(jsup) - gphiv(kbdy,1) ) * COS(zrad*Ug(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamv(kbdy,1) ) * ( gphiv(kbdy,1) - lat(jinf) ) * COS(zrad*Ug(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiv(kbdy,1) ) * COS(zrad*Ug(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( gphiv(kbdy,1) - lat(jinf) ) * COS(zrad*Ug(isup,jsup)) * mskfes(isup,jsup)
 
         Usn =      ( zlon(isup) - zglamv(kbdy,1) ) * ( lat(jsup) - gphiv(kbdy,1) ) * SIN(zrad*Ug(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamv(kbdy,1) ) * ( gphiv(kbdy,1) - lat(jinf) ) * SIN(zrad*Ug(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiv(kbdy,1) ) * SIN(zrad*Ug(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( gphiv(kbdy,1) - lat(jinf) ) * SIN(zrad*Ug(isup,jsup)) * mskfes(isup,jsup)
  
         Va_bdy =   ( zlon(isup) - zglamv(kbdy,1) ) * ( lat(jsup) - gphiv(kbdy,1) ) * Va(iinf,jinf) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamv(kbdy,1) ) * ( gphiv(kbdy,1) - lat(jinf) ) * Va(iinf,jsup) * mskfes(iinf,jsup)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiv(kbdy,1) ) * Va(isup,jinf) * mskfes(isup,jinf)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( gphiv(kbdy,1) - lat(jinf) ) * Va(isup,jsup) * mskfes(isup,jsup)
       
         Vcs =      ( zlon(isup) - zglamv(kbdy,1) ) * ( lat(jsup) - gphiv(kbdy,1) ) * COS(zrad*Vg(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamv(kbdy,1) ) * ( gphiv(kbdy,1) - lat(jinf) ) * COS(zrad*Vg(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiv(kbdy,1) ) * COS(zrad*Vg(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( gphiv(kbdy,1) - lat(jinf) ) * COS(zrad*Vg(isup,jsup)) * mskfes(isup,jsup)
 
         Vsn =      ( zlon(isup) - zglamv(kbdy,1) ) * ( lat(jsup) - gphiv(kbdy,1) ) * SIN(zrad*Vg(iinf,jinf)) * mskfes(iinf,jinf)  &
         &        + ( zlon(isup) - zglamv(kbdy,1) ) * ( gphiv(kbdy,1) - lat(jinf) ) * SIN(zrad*Vg(iinf,jsup)) * mskfes(iinf,jsup)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiv(kbdy,1) ) * SIN(zrad*Vg(isup,jinf)) * mskfes(isup,jinf)  &
         &        + ( zglamv(kbdy,1) - zlon(iinf) ) * ( gphiv(kbdy,1) - lat(jinf) ) * SIN(zrad*Vg(isup,jsup)) * mskfes(isup,jsup)
 
         angle  = ATAN2(    gphiv_CHLD( MIN(mxG,nbiv(kbdy,1)+1), nbjv(kbdy,1) ) -  gphiv_CHLD( MAX(1,nbiv(kbdy,1)-1), nbjv(kbdy,1) )  , &
         &               ( zglamv_CHLD( MIN(mxG,nbiv(kbdy,1)+1), nbjv(kbdy,1) ) - zglamv_CHLD( MAX(1,nbiv(kbdy,1)-1), nbjv(kbdy,1) ) )  &
         &               * cos( zrad * gphiv_CHLD(nbiv(kbdy,1),nbjv(kbdy,1)) )  ) 
 
         div =   ( zlon(isup) - zglamv(kbdy,1) ) * ( lat(jsup) - gphiv(kbdy,1) ) * mskfes(iinf,jinf)  &
           &   + ( zlon(isup) - zglamv(kbdy,1) ) * ( gphiv(kbdy,1) - lat(jinf) ) * mskfes(iinf,jsup)  &
           &   + ( zglamv(kbdy,1) - zlon(iinf) ) * ( lat(jsup) - gphiv(kbdy,1) ) * mskfes(isup,jinf)  &
           &   + ( zglamv(kbdy,1) - zlon(iinf) ) * ( gphiv(kbdy,1) - lat(jinf) ) * mskfes(isup,jsup)
       
         if ( abs(div) .gt. 1.e-12 ) then
           Ua_bdy = Ua_bdy / div
           Ug_bdy = ATAN2( Usn/div, Ucs/div ) / zrad !- to have a proper "modulo" interpolation of the phase 
           Va_bdy = Va_bdy / div
           Vg_bdy = ATAN2( Vsn/div, Vcs/div ) / zrad !- to have a proper "modulo" interpolation of the phase
           exit
         endif

       enddo !- knb

       if ( .not. abs(div) .gt. 1.e-12 ) then
         write(*,*) '~!@#$%^ INCREASE NB OF LOOPS ON knb (bdyV)  >>>>>>>>>>  stop !!'
         write(*,*) kbdy, iinf, jinf, Va_bdy, Vg_bdy
         stop
       endif

     else
    
        Ua_bdy = 0.0
        Ug_bdy = 0.0
        Va_bdy = 0.0
        Vg_bdy = 0.0
        angle  = 0.0
   
     endif !--  if ( tmask(nbit(kbdy,1),nbjt(kbdy,1),1) .ne. 0 )
  
     v1(kbdy,1) = - Ua_bdy * COS( zrad * Ug_bdy ) * SIN(angle) + Va_bdy * COS( zrad * Vg_bdy ) * COS(angle) !- real part
     v2(kbdy,1) = - Ua_bdy * SIN( zrad * Ug_bdy ) * SIN(angle) + Va_bdy * SIN( zrad * Vg_bdy ) * COS(angle) !- imaginary part
   
   enddo

   !========================================
  
   !------------------------
   ! Create gridT bdy file :
   
   write(file_out_T,401) TRIM(config_dir), TRIM(config), TRIM(harmstr)
   
   write(*,*) 'Creating ', TRIM(file_out_T)
   
   status = NF90_CREATE(TRIM(file_out_T),NF90_NOCLOBBER,fidoutT) ; call erreur(status,.TRUE.,'create file_out_T')                     
   
   status = NF90_DEF_DIM(fidoutT,"yb",myb,dimID_yb)   ; call erreur(status,.TRUE.,"def_dimID_yb")
   status = NF90_DEF_DIM(fidoutT,"xb",mxbt,dimID_xbt) ; call erreur(status,.TRUE.,"def_dimID_xb")
   
   status = NF90_DEF_VAR(fidoutT,"z2",NF90_FLOAT,(/dimID_xbt,dimID_yb/),z2_ID)   ; call erreur(status,.TRUE.,"def_var_z2_ID")
   status = NF90_DEF_VAR(fidoutT,"z1",NF90_FLOAT,(/dimID_xbt,dimID_yb/),z1_ID)   ; call erreur(status,.TRUE.,"def_var_z1_ID")
   status = NF90_DEF_VAR(fidoutT,"nbrt",NF90_INT,(/dimID_xbt,dimID_yb/),nbrt_ID) ; call erreur(status,.TRUE.,"def_var_nbrt_ID")
   status = NF90_DEF_VAR(fidoutT,"nbjt",NF90_INT,(/dimID_xbt,dimID_yb/),nbjt_ID) ; call erreur(status,.TRUE.,"def_var_nbjt_ID")
   status = NF90_DEF_VAR(fidoutT,"nbit",NF90_INT,(/dimID_xbt,dimID_yb/),nbit_ID) ; call erreur(status,.TRUE.,"def_var_nbit_ID")
   
   status = NF90_PUT_ATT(fidoutT,z2_ID,"_FillValue",0.)                       ; call erreur(status,.TRUE.,"put_att_z2_ID")
   status = NF90_PUT_ATT(fidoutT,z2_ID,"long_name","(sin) tidal elevation")   ; call erreur(status,.TRUE.,"put_att_z2_ID")
   status = NF90_PUT_ATT(fidoutT,z2_ID,"units","m")                           ; call erreur(status,.TRUE.,"put_att_z2_ID")
   status = NF90_PUT_ATT(fidoutT,z1_ID,"_FillValue",0.)                       ; call erreur(status,.TRUE.,"put_att_z1_ID")
   status = NF90_PUT_ATT(fidoutT,z1_ID,"long_name","(cos) tidal elevation")   ; call erreur(status,.TRUE.,"put_att_z1_ID")
   status = NF90_PUT_ATT(fidoutT,z1_ID,"units","m")                           ; call erreur(status,.TRUE.,"put_att_z1_ID")
   status = NF90_PUT_ATT(fidoutT,nbrt_ID,"long_name","bdy discrete distance") ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
   status = NF90_PUT_ATT(fidoutT,nbrt_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
   status = NF90_PUT_ATT(fidoutT,nbjt_ID,"long_name","bdy j index")           ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
   status = NF90_PUT_ATT(fidoutT,nbjt_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
   status = NF90_PUT_ATT(fidoutT,nbit_ID,"long_name","bdy i index")           ; call erreur(status,.TRUE.,"put_att_nbit_ID")
   status = NF90_PUT_ATT(fidoutT,nbit_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbit_ID")
   
   status = NF90_PUT_ATT(fidoutT,NF90_GLOBAL,"history","Created using extract_bdy_tides.f90")
   status = NF90_PUT_ATT(fidoutT,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
   call erreur(status,.TRUE.,"put_att_GLOBAL")
   
   status = NF90_ENDDEF(fidoutT) ; call erreur(status,.TRUE.,"fin_definition") 
   
   status = NF90_PUT_VAR(fidoutT,z2_ID,z2)     ; call erreur(status,.TRUE.,"var_z2_ID")
   status = NF90_PUT_VAR(fidoutT,z1_ID,z1)     ; call erreur(status,.TRUE.,"var_z1_ID")
   status = NF90_PUT_VAR(fidoutT,nbrt_ID,nbrt) ; call erreur(status,.TRUE.,"var_nbrdta_ID")
   status = NF90_PUT_VAR(fidoutT,nbjt_ID,nbjt) ; call erreur(status,.TRUE.,"var_nbjdta_ID")
   status = NF90_PUT_VAR(fidoutT,nbit_ID,nbit) ; call erreur(status,.TRUE.,"var_nbidta_ID")
   
   status = NF90_CLOSE(fidoutT) ; call erreur(status,.TRUE.,"final")         
   
   !------------------------
   ! Create gridU bdy file :
   
   write(file_out_U,402) TRIM(config_dir), TRIM(config), TRIM(harmstr)
   
   write(*,*) 'Creating ', TRIM(file_out_U)
   
   status = NF90_CREATE(TRIM(file_out_U),NF90_NOCLOBBER,fidoutU) ; call erreur(status,.TRUE.,'create file_out_U')                     
   
   status = NF90_DEF_DIM(fidoutU,"yb",myb,dimID_yb)   ; call erreur(status,.TRUE.,"def_dimID_yb")
   status = NF90_DEF_DIM(fidoutU,"xb",mxbu,dimID_xbu) ; call erreur(status,.TRUE.,"def_dimID_xb")
   
   status = NF90_DEF_VAR(fidoutU,"u2",NF90_FLOAT,(/dimID_xbu,dimID_yb/),u2_ID)   ; call erreur(status,.TRUE.,"def_var_u2_ID")
   status = NF90_DEF_VAR(fidoutU,"u1",NF90_FLOAT,(/dimID_xbu,dimID_yb/),u1_ID)   ; call erreur(status,.TRUE.,"def_var_u1_ID")
   status = NF90_DEF_VAR(fidoutU,"nbru",NF90_INT,(/dimID_xbu,dimID_yb/),nbru_ID) ; call erreur(status,.TRUE.,"def_var_nbru_ID")
   status = NF90_DEF_VAR(fidoutU,"nbju",NF90_INT,(/dimID_xbu,dimID_yb/),nbju_ID) ; call erreur(status,.TRUE.,"def_var_nbju_ID")
   status = NF90_DEF_VAR(fidoutU,"nbiu",NF90_INT,(/dimID_xbu,dimID_yb/),nbiu_ID) ; call erreur(status,.TRUE.,"def_var_nbiu_ID")
   
   status = NF90_PUT_ATT(fidoutU,u2_ID,"_FillValue",0.)                         ; call erreur(status,.TRUE.,"put_att_u2_ID")
   status = NF90_PUT_ATT(fidoutU,u2_ID,"long_name","(sin) tidal east velocity") ; call erreur(status,.TRUE.,"put_att_u2_ID")
   status = NF90_PUT_ATT(fidoutU,u2_ID,"units","m/s")                           ; call erreur(status,.TRUE.,"put_att_u2_ID")
   status = NF90_PUT_ATT(fidoutU,u1_ID,"_FillValue",0.)                         ; call erreur(status,.TRUE.,"put_att_u1_ID")
   status = NF90_PUT_ATT(fidoutU,u1_ID,"long_name","(cos) tidal east velocity") ; call erreur(status,.TRUE.,"put_att_u1_ID")
   status = NF90_PUT_ATT(fidoutU,u1_ID,"units","m/s")                           ; call erreur(status,.TRUE.,"put_att_u1_ID")
   status = NF90_PUT_ATT(fidoutU,nbru_ID,"long_name","bdy discrete distance")   ; call erreur(status,.TRUE.,"put_att_nbru_ID")
   status = NF90_PUT_ATT(fidoutU,nbru_ID,"units","unitless")                    ; call erreur(status,.TRUE.,"put_att_nbru_ID")
   status = NF90_PUT_ATT(fidoutU,nbju_ID,"long_name","bdy j index")             ; call erreur(status,.TRUE.,"put_att_nbju_ID")
   status = NF90_PUT_ATT(fidoutU,nbju_ID,"units","unitless")                    ; call erreur(status,.TRUE.,"put_att_nbju_ID")
   status = NF90_PUT_ATT(fidoutU,nbiu_ID,"long_name","bdy i index")             ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
   status = NF90_PUT_ATT(fidoutU,nbiu_ID,"units","unitless")                    ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
   
   status = NF90_PUT_ATT(fidoutU,NF90_GLOBAL,"history","Created using extract_bdy_tides.f90")
   status = NF90_PUT_ATT(fidoutU,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
   call erreur(status,.TRUE.,"put_att_GLOBAL")
   
   status = NF90_ENDDEF(fidoutU) ; call erreur(status,.TRUE.,"fin_definition") 
   
   status = NF90_PUT_VAR(fidoutU,u2_ID,u2)     ; call erreur(status,.TRUE.,"var_u2_ID")
   status = NF90_PUT_VAR(fidoutU,u1_ID,u1)     ; call erreur(status,.TRUE.,"var_u1_ID")
   status = NF90_PUT_VAR(fidoutU,nbru_ID,nbru) ; call erreur(status,.TRUE.,"var_nbru_ID")
   status = NF90_PUT_VAR(fidoutU,nbju_ID,nbju) ; call erreur(status,.TRUE.,"var_nbju_ID")
   status = NF90_PUT_VAR(fidoutU,nbiu_ID,nbiu) ; call erreur(status,.TRUE.,"var_nbiu_ID")
   
   status = NF90_CLOSE(fidoutU) ; call erreur(status,.TRUE.,"final")         
   
   !------------------------
   ! Create gridV bdy file :
   
   write(file_out_V,403) TRIM(config_dir), TRIM(config), TRIM(harmstr)
   
   write(*,*) 'Creating ', TRIM(file_out_V)
   
   status = NF90_CREATE(TRIM(file_out_V),NF90_NOCLOBBER,fidoutV) ; call erreur(status,.TRUE.,'create file_out_V')                     
   
   status = NF90_DEF_DIM(fidoutV,"yb",myb,dimID_yb)   ; call erreur(status,.TRUE.,"def_dimID_yb")
   status = NF90_DEF_DIM(fidoutV,"xb",mxbv,dimID_xbv) ; call erreur(status,.TRUE.,"def_dimID_xb")
   
   status = NF90_DEF_VAR(fidoutV,"v2",NF90_FLOAT,(/dimID_xbv,dimID_yb/),v2_ID)   ; call erreur(status,.TRUE.,"def_var_v2_ID")
   status = NF90_DEF_VAR(fidoutV,"v1",NF90_FLOAT,(/dimID_xbv,dimID_yb/),v1_ID)   ; call erreur(status,.TRUE.,"def_var_v1_ID")
   status = NF90_DEF_VAR(fidoutV,"nbrv",NF90_INT,(/dimID_xbv,dimID_yb/),nbrv_ID) ; call erreur(status,.TRUE.,"def_var_nbrv_ID")
   status = NF90_DEF_VAR(fidoutV,"nbjv",NF90_INT,(/dimID_xbv,dimID_yb/),nbjv_ID) ; call erreur(status,.TRUE.,"def_var_nbjv_ID")
   status = NF90_DEF_VAR(fidoutV,"nbiv",NF90_INT,(/dimID_xbv,dimID_yb/),nbiv_ID) ; call erreur(status,.TRUE.,"def_var_nbiv_ID")
   
   status = NF90_PUT_ATT(fidoutV,v2_ID,"_FillValue",0.)                         ; call erreur(status,.TRUE.,"put_att_v2_ID")
   status = NF90_PUT_ATT(fidoutV,v2_ID,"long_name","(sin) tidal north velocity"); call erreur(status,.TRUE.,"put_att_v2_ID")
   status = NF90_PUT_ATT(fidoutV,v2_ID,"units","m/s")                           ; call erreur(status,.TRUE.,"put_att_v2_ID")
   status = NF90_PUT_ATT(fidoutV,v1_ID,"_FillValue",0.)                         ; call erreur(status,.TRUE.,"put_att_v1_ID")
   status = NF90_PUT_ATT(fidoutV,v1_ID,"long_name","(cos) tidal north velocity"); call erreur(status,.TRUE.,"put_att_v1_ID")
   status = NF90_PUT_ATT(fidoutV,v1_ID,"units","m/s")                           ; call erreur(status,.TRUE.,"put_att_v1_ID")
   status = NF90_PUT_ATT(fidoutV,nbrv_ID,"long_name","bdy discrete distance")   ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
   status = NF90_PUT_ATT(fidoutV,nbrv_ID,"units","unitless")                    ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
   status = NF90_PUT_ATT(fidoutV,nbjv_ID,"long_name","bdy j index")             ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
   status = NF90_PUT_ATT(fidoutV,nbjv_ID,"units","unitless")                    ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
   status = NF90_PUT_ATT(fidoutV,nbiv_ID,"long_name","bdy i index")             ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
   status = NF90_PUT_ATT(fidoutV,nbiv_ID,"units","unitless")                    ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
   
   status = NF90_PUT_ATT(fidoutV,NF90_GLOBAL,"history","Created using extract_bdy_tides.f90")
   status = NF90_PUT_ATT(fidoutV,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
   call erreur(status,.TRUE.,"put_att_GLOBAL")
   
   status = NF90_ENDDEF(fidoutV) ; call erreur(status,.TRUE.,"fin_definition") 
   
   status = NF90_PUT_VAR(fidoutV,v2_ID,v2)     ; call erreur(status,.TRUE.,"var_v2_ID")
   status = NF90_PUT_VAR(fidoutV,v1_ID,v1)     ; call erreur(status,.TRUE.,"var_v1_ID")
   status = NF90_PUT_VAR(fidoutV,nbrv_ID,nbrv) ; call erreur(status,.TRUE.,"var_nbrv_ID")
   status = NF90_PUT_VAR(fidoutV,nbjv_ID,nbjv) ; call erreur(status,.TRUE.,"var_nbjv_ID")
   status = NF90_PUT_VAR(fidoutV,nbiv_ID,nbiv) ; call erreur(status,.TRUE.,"var_nbiv_ID")
   
   status = NF90_CLOSE(fidoutV) ; call erreur(status,.TRUE.,"final")         
  
   !-----
 
   DEALLOCATE( z1, z2, u1, u2, v1, v2, Ha, Hg, Ua, Ug, Va, Vg, mskfes, lat, lon, zlon )

ENDDO  !-- kharm

end program modif

!========================================

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
