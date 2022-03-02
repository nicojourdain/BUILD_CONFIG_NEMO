program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Used to build netcdf coordinate file for BDY
!
! 0- Initializations
! 1a- Read coordinates of the entire domain
! 1b- Read mask of the entire domain
! 2- Calculate BDY dimensions
! 3- Exctract coordinates along BDY
! 4- Write BDY coordinates in a netcdf file
!
! hisotry: - Mar. 2015: initial verison (N. Jourdain, CNRS-LGGE)
!          - Feb. 2017: use of namelist (N. Jourdain, CNRS-IGE)
!          - May  2021: remove 3 consecutive land points (useful for land-processor masking in NEMO4)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /bdy/ nn_bdy_east, nn_bdy_west, nn_bdy_north, nn_bdy_south, nn_harm
namelist /bdy_east/  ii_bdy_east,  j1_bdy_east,  j2_bdy_east
namelist /bdy_west/  ii_bdy_west,  j1_bdy_west,  j2_bdy_west
namelist /bdy_north/ jj_bdy_north, i1_bdy_north, i2_bdy_north
namelist /bdy_south/ jj_bdy_south, i1_bdy_south, i2_bdy_south

CHARACTER(LEN=50)                    :: config
CHARACTER(LEN=150)                   :: config_dir
INTEGER                              :: nn_bdy_east, nn_bdy_west, nn_bdy_north, nn_bdy_south, nn_harm
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: ii_bdy_east,  j1_bdy_east,  j2_bdy_east,               &
&                                       ii_bdy_west,  j1_bdy_west,  j2_bdy_west,               &
&                                       jj_bdy_north, i1_bdy_north, i2_bdy_north,              &
&                                       jj_bdy_south, i1_bdy_south, i2_bdy_south

INTEGER                              :: fidA, status, dimID_xbv, dimID_xbu, dimID_xbt, dimID_yb, mxbv, mxbu,    &
&                                       mxbt, myb, e2v_ID, e1v_ID, gphiv_ID, glamv_ID, e2u_ID, e1u_ID, gphiu_ID,&
&                                       glamu_ID, e2t_ID, e1t_ID, gphit_ID, glamt_ID, nbrv_ID, nbjv_ID, nbiv_ID,&
&                                       nbru_ID, nbju_ID, nbiu_ID, nbrt_ID, nbjt_ID, nbit_ID, fidM, dimID_x,    &
&                                       dimID_y, k, kt, ku, kv, kt0, ku0, kv0, mx, my, kkbdy, k1, k2, ndouble,  &
&                                       fidJ, vmaskutil_ID, umaskutil_ID, tmaskutil_ID, nn_perio
CHARACTER(LEN=150)                   :: file_in_coord_CHLD, file_in_mask_CHLD, file_out                     
INTEGER*1,ALLOCATABLE,DIMENSION(:,:) :: vmaskutil, umaskutil, tmaskutil
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbrv, nbjv, nbiv, nbru, nbju, nbiu, nbrt, nbjt, nbit    
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbjv_tmp, nbiv_tmp, nbju_tmp, nbiu_tmp, nbjt_tmp, nbit_tmp    
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: e1t_bdy, e2t_bdy, e1u_bdy, e2u_bdy, e1v_bdy, e2v_bdy,                   &
&                                       gphit_bdy, glamt_bdy, gphiu_bdy, glamu_bdy, gphiv_bdy, glamv_bdy
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: e1t, e2t, e1u, e2u, e1v, e2v, gphit, glamt, gphiu, glamu, gphiv, glamv

!=================================================================================
! 0- Initializations 
!=============================================================================mask
!- default values :
nn_bdy_east  = 0
nn_bdy_west  = 0
nn_bdy_north = 0
nn_bdy_south = 0
nn_harm      = 0

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=bdy)
CLOSE(1)

write(*,*) 'nb of segments constituting the eastern  bdy = ', nn_bdy_east
write(*,*) 'nb of segments constituting the western  bdy = ', nn_bdy_west
write(*,*) 'nb of segments constituting the northern bdy = ', nn_bdy_north
write(*,*) 'nb of segments constituting the southern bdy = ', nn_bdy_south

IF (nn_bdy_east  .gt. 0) ALLOCATE( ii_bdy_east(nn_bdy_east),    j1_bdy_east(nn_bdy_east),   j2_bdy_east(nn_bdy_east)   )
IF (nn_bdy_west  .gt. 0) ALLOCATE( ii_bdy_west(nn_bdy_west),    j1_bdy_west(nn_bdy_west),   j2_bdy_west(nn_bdy_west)   )
IF (nn_bdy_north .gt. 0) ALLOCATE( jj_bdy_north(nn_bdy_north),  i1_bdy_north(nn_bdy_north), i2_bdy_north(nn_bdy_north) )
IF (nn_bdy_south .gt. 0) ALLOCATE( jj_bdy_south(nn_bdy_south),  i1_bdy_south(nn_bdy_south), i2_bdy_south(nn_bdy_south) )

OPEN (UNIT=1, FILE='namelist_pre' )
IF (nn_bdy_east  .gt. 0) READ (UNIT=1, NML=bdy_east)
IF (nn_bdy_west  .gt. 0) READ (UNIT=1, NML=bdy_west)
IF (nn_bdy_north .gt. 0) READ (UNIT=1, NML=bdy_north)
IF (nn_bdy_south .gt. 0) READ (UNIT=1, NML=bdy_south)
CLOSE(1)

write(file_in_coord_CHLD,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/coordinates_',a,'.nc')
write(file_in_mask_CHLD,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/mesh_mask_',a,'.nc')
write(file_out,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coordinates_bdy_',a,'.nc')

!=================================================================================
! 1a- Read coordinates of the entire domain
!=================================================================================

write(*,*) 'Reading ', TRIM(file_in_coord_CHLD)

status = NF90_OPEN(TRIM(file_in_coord_CHLD),0,fidA) ; call erreur(status,.TRUE.,"open coordinates") 

status = NF90_INQ_DIMID(fidA,"x",dimID_x)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"X",dimID_x) 
call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidA,"y",dimID_y)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidA,"Y",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")

status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx) ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my) ; call erreur(status,.TRUE.,"inq_dim_y")

ALLOCATE(  e1t(mx,my), e2t(mx,my)  )
ALLOCATE(  e1u(mx,my), e2u(mx,my)  )
ALLOCATE(  e1v(mx,my), e2v(mx,my)  )
ALLOCATE(  gphit(mx,my), glamt(mx,my)  )
ALLOCATE(  gphiu(mx,my), glamu(mx,my)  )
ALLOCATE(  gphiv(mx,my), glamv(mx,my)  )

status = NF90_INQ_VARID(fidA,"e1t",e1t_ID)     ; call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidA,"e2t",e2t_ID)     ; call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidA,"e1u",e1u_ID)     ; call erreur(status,.TRUE.,"inq_e1u_ID")
status = NF90_INQ_VARID(fidA,"e2u",e2u_ID)     ; call erreur(status,.TRUE.,"inq_e2u_ID")
status = NF90_INQ_VARID(fidA,"e1v",e1v_ID)     ; call erreur(status,.TRUE.,"inq_e1v_ID")
status = NF90_INQ_VARID(fidA,"e2v",e2v_ID)     ; call erreur(status,.TRUE.,"inq_e2v_ID")
status = NF90_INQ_VARID(fidA,"gphit",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidA,"glamt",glamt_ID) ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidA,"gphiu",gphiu_ID) ; call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidA,"glamu",glamu_ID) ; call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidA,"gphiv",gphiv_ID) ; call erreur(status,.TRUE.,"inq_gphiv_ID")
status = NF90_INQ_VARID(fidA,"glamv",glamv_ID) ; call erreur(status,.TRUE.,"inq_glamv_ID")

status = NF90_GET_VAR(fidA,e1t_ID,e1t)     ; call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidA,e2t_ID,e2t)     ; call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidA,e1u_ID,e1u)     ; call erreur(status,.TRUE.,"getvar_e1u")
status = NF90_GET_VAR(fidA,e2u_ID,e2u)     ; call erreur(status,.TRUE.,"getvar_e2u")
status = NF90_GET_VAR(fidA,e1v_ID,e1v)     ; call erreur(status,.TRUE.,"getvar_e1v")
status = NF90_GET_VAR(fidA,e2v_ID,e2v)     ; call erreur(status,.TRUE.,"getvar_e2v")
status = NF90_GET_VAR(fidA,gphit_ID,gphit) ; call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidA,glamt_ID,glamt) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidA,gphiu_ID,gphiu) ; call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidA,glamu_ID,glamu) ; call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidA,gphiv_ID,gphiv) ; call erreur(status,.TRUE.,"getvar_gphiv")
status = NF90_GET_VAR(fidA,glamv_ID,glamv) ; call erreur(status,.TRUE.,"getvar_glamv")

status = NF90_CLOSE(fidA) ; call erreur(status,.TRUE.,"close file")     

!=================================================================================
! 1b- Read mask of the entire domain
!=================================================================================

write(*,*) 'Reading ', TRIM(file_in_mask_CHLD)

status = NF90_OPEN(TRIM(file_in_mask_CHLD),0,fidJ) ; call erreur(status,.TRUE.,"open coordinates") 

ALLOCATE(  tmaskutil(mx,my)  )
ALLOCATE(  umaskutil(mx,my)  )
ALLOCATE(  vmaskutil(mx,my)  )

status = NF90_INQ_VARID(fidJ,"tmaskutil",tmaskutil_ID)     ; call erreur(status,.TRUE.,"inq_tmaskutil_ID")
status = NF90_INQ_VARID(fidJ,"umaskutil",umaskutil_ID)     ; call erreur(status,.TRUE.,"inq_umaskutil_ID")
status = NF90_INQ_VARID(fidJ,"vmaskutil",vmaskutil_ID)     ; call erreur(status,.TRUE.,"inq_vmaskutil_ID")

status = NF90_GET_VAR(fidJ,tmaskutil_ID,tmaskutil)     ; call erreur(status,.TRUE.,"getvar_tmaskutil")
status = NF90_GET_VAR(fidJ,umaskutil_ID,umaskutil)     ; call erreur(status,.TRUE.,"getvar_umaskutil")
status = NF90_GET_VAR(fidJ,vmaskutil_ID,vmaskutil)     ; call erreur(status,.TRUE.,"getvar_vmaskutil")

status = NF90_CLOSE(fidJ) ; call erreur(status,.TRUE.,"close file")     

!=================================================================================
! 2- Calculate BDY dimensions
!=================================================================================

! Check if grid is east-west periodic :
nn_perio = 0
if (nn_bdy_north .gt. 0 .and. i1_bdy_north(1) .eq. 1 .and. i2_bdy_north(1) .eq. mx ) nn_perio = 1
if (nn_bdy_south .gt. 0 .and. i1_bdy_south(1) .eq. 1 .and. i2_bdy_south(1) .eq. mx ) nn_perio = 1
if ( nn_perio .eq. 1 ) write(*,*) 'NB: ASSUMING EAST-WEST PERIODIC DOMAIN !!'

mxbt=0
mxbu=0
mxbv=0

do kkbdy=1,nn_bdy_east
  mxbt=mxbt+j2_bdy_east(kkbdy)-j1_bdy_east(kkbdy)+1
  mxbu=mxbu+j2_bdy_east(kkbdy)-j1_bdy_east(kkbdy)+1
  mxbv=mxbv+j2_bdy_east(kkbdy)-j1_bdy_east(kkbdy)
  !WED12: mxbv=mxbv+j2_bdy_east(kkbdy)-j1_bdy_east(kkbdy)+MIN(2,j1_bdy_east(kkbdy))
enddo

do kkbdy=1,nn_bdy_west
  mxbt=mxbt+j2_bdy_west(kkbdy)-j1_bdy_west(kkbdy)+1
  mxbu=mxbu+j2_bdy_west(kkbdy)-j1_bdy_west(kkbdy)+1
  mxbv=mxbv+j2_bdy_west(kkbdy)-j1_bdy_west(kkbdy)
  ! WED12: mxbv=mxbv+j2_bdy_west(kkbdy)-j1_bdy_west(kkbdy)+MIN(2,j1_bdy_west(kkbdy))
enddo

do kkbdy=1,nn_bdy_north
  mxbt=mxbt+i2_bdy_north(kkbdy)-i1_bdy_north(kkbdy)+1
  mxbu=mxbu+i2_bdy_north(kkbdy)-i1_bdy_north(kkbdy)+nn_perio
  !WED12: mxbu=mxbu+i2_bdy_north(kkbdy)-i1_bdy_north(kkbdy)+MIN(2,i1_bdy_north(kkbdy))
  mxbv=mxbv+i2_bdy_north(kkbdy)-i1_bdy_north(kkbdy)+1
enddo

do kkbdy=1,nn_bdy_south
  mxbt=mxbt+i2_bdy_south(kkbdy)-i1_bdy_south(kkbdy)+1
  mxbu=mxbu+i2_bdy_south(kkbdy)-i1_bdy_south(kkbdy)+nn_perio
  !WED12: mxbu=mxbu+i2_bdy_south(kkbdy)-i1_bdy_south(kkbdy)+MIN(2,i1_bdy_south(kkbdy))
  mxbv=mxbv+i2_bdy_south(kkbdy)-i1_bdy_south(kkbdy)+1
enddo

write(*,*) 'mxbt = ', mxbt
write(*,*) 'mxbu = ', mxbu
write(*,*) 'mxbv = ', mxbv

myb=1 !- degenerated dimension

ALLOCATE(  nbit_tmp(mxbt,myb), nbjt_tmp(mxbt,myb)  )
ALLOCATE(  nbiu_tmp(mxbu,myb), nbju_tmp(mxbu,myb)  )
ALLOCATE(  nbiv_tmp(mxbv,myb), nbjv_tmp(mxbv,myb)  )

!=================================================================================
! 3- Exctract coordinates along BDY:
!    NB: from NEMO4, removing any point with 3 consecutive land values
!        (t/u/vmaskutil=0)
!=================================================================================

write(*,*) 'Exctracting coordinates along BDY...'

kt=0
ku=0
kv=0

do kkbdy=1,nn_bdy_east
  do k=j1_bdy_east(kkbdy),j2_bdy_east(kkbdy)
    kt=kt+1
    if ( SUM(tmaskutil(ii_bdy_east(kkbdy),MAX(k-1,1):MIN(k+1,my))).ne.0 ) then
      nbit_tmp(kt,1)=ii_bdy_east(kkbdy)
      nbjt_tmp(kt,1)=k
    else
      nbit_tmp(kt,1)=999999
      nbjt_tmp(kt,1)=999999
    endif
  enddo
  do k=j1_bdy_east(kkbdy),j2_bdy_east(kkbdy)
    ku=ku+1
    if ( SUM(umaskutil(ii_bdy_east(kkbdy)-1,MAX(k-1,1):MIN(k+1,my))).ne.0 ) then
      nbiu_tmp(ku,1)=ii_bdy_east(kkbdy)-1
      nbju_tmp(ku,1)=k
    else
      nbiu_tmp(ku,1)=999999
      nbju_tmp(ku,1)=999999
    endif
  enddo
  do k=j1_bdy_east(kkbdy),j2_bdy_east(kkbdy)-1
  !WED12: do k=MAX(1,j1_bdy_east(kkbdy)-1),j2_bdy_east(kkbdy)
    kv=kv+1
    if ( SUM(vmaskutil(ii_bdy_east(kkbdy),MAX(k-1,1):MIN(k+1,my))).ne.0 ) then
      nbiv_tmp(kv,1)=ii_bdy_east(kkbdy)
      nbjv_tmp(kv,1)=k
    else
      nbiv_tmp(kv,1)=999999
      nbjv_tmp(kv,1)=999999
    endif
  enddo
enddo

do kkbdy=1,nn_bdy_west
  do k=j1_bdy_west(kkbdy),j2_bdy_west(kkbdy)
    kt=kt+1
    if ( SUM(tmaskutil(ii_bdy_west(kkbdy),MAX(k-1,1):MIN(k+1,my))).ne.0 ) then
      nbit_tmp(kt,1)=ii_bdy_west(kkbdy)
      nbjt_tmp(kt,1)=k
    else
      nbit_tmp(kt,1)=999999
      nbjt_tmp(kt,1)=999999
    endif
  enddo
  do k=j1_bdy_west(kkbdy),j2_bdy_west(kkbdy)
    ku=ku+1
    if ( SUM(umaskutil(ii_bdy_west(kkbdy),MAX(k-1,1):MIN(k+1,my))).ne.0 ) then
      nbiu_tmp(ku,1)=ii_bdy_west(kkbdy)
      nbju_tmp(ku,1)=k
    else
      nbiu_tmp(ku,1)=999999
      nbju_tmp(ku,1)=999999
    endif
  enddo
  do k=j1_bdy_west(kkbdy),j2_bdy_west(kkbdy)-1
  !WED12: do k=MAX(1,j1_bdy_west(kkbdy)-1),j2_bdy_west(kkbdy)
    kv=kv+1
    if ( SUM(vmaskutil(ii_bdy_west(kkbdy),MAX(k-1,1):MIN(k+1,my))).ne.0 ) then
      nbiv_tmp(kv,1)=ii_bdy_west(kkbdy)
      nbjv_tmp(kv,1)=k
    else
      nbiv_tmp(kv,1)=999999
      nbjv_tmp(kv,1)=999999
    endif
  enddo
enddo

do kkbdy=1,nn_bdy_north
  do k=i1_bdy_north(kkbdy),i2_bdy_north(kkbdy)
    kt=kt+1
    if ( SUM(tmaskutil(MAX(k-1,1):MIN(k+1,mx),jj_bdy_north(kkbdy))).ne.0 ) then
      nbit_tmp(kt,1)=k
      nbjt_tmp(kt,1)=jj_bdy_north(kkbdy)
    else
      nbit_tmp(kt,1)=999999
      nbjt_tmp(kt,1)=999999
    endif
  enddo
  do k=i1_bdy_north(kkbdy),i2_bdy_north(kkbdy)-1+nn_perio
  !WED12: do k=MAX(1,i1_bdy_north(kkbdy)-1),i2_bdy_north(kkbdy)
    ku=ku+1
    if ( SUM(umaskutil(MAX(k-1,1):MIN(k+1,mx),jj_bdy_north(kkbdy))).ne.0 ) then
      nbiu_tmp(ku,1)=k
      nbju_tmp(ku,1)=jj_bdy_north(kkbdy)
    else
      nbiu_tmp(ku,1)=999999
      nbju_tmp(ku,1)=999999
    endif
  enddo
  do k=i1_bdy_north(kkbdy),i2_bdy_north(kkbdy)
    kv=kv+1
    if ( SUM(vmaskutil(MAX(k-1,1):MIN(k+1,mx),jj_bdy_north(kkbdy)-1)).ne.0 ) then
      nbiv_tmp(kv,1)=k
      nbjv_tmp(kv,1)=jj_bdy_north(kkbdy)-1
    else
      nbiv_tmp(kv,1)=999999
      nbjv_tmp(kv,1)=999999
    endif
  enddo
enddo

do kkbdy=1,nn_bdy_south
  do k=i1_bdy_south(kkbdy),i2_bdy_south(kkbdy)
    kt=kt+1
    if ( SUM(tmaskutil(MAX(k-1,1):MIN(k+1,mx),jj_bdy_south(kkbdy))).ne.0 ) then
      nbit_tmp(kt,1)=k
      nbjt_tmp(kt,1)=jj_bdy_south(kkbdy)
    else
      nbit_tmp(kt,1)=999999
      nbjt_tmp(kt,1)=999999
    endif
  enddo
  do k=i1_bdy_south(kkbdy),i2_bdy_south(kkbdy)-1+nn_perio
  !WED12: do k=MAX(1,i1_bdy_south(kkbdy)-1),i2_bdy_south(kkbdy)
    ku=ku+1
    if ( SUM(umaskutil(MAX(k-1,1):MIN(k+1,mx),jj_bdy_south(kkbdy))).ne.0 ) then
      nbiu_tmp(ku,1)=k
      nbju_tmp(ku,1)=jj_bdy_south(kkbdy)
    else
      nbiu_tmp(ku,1)=999999
      nbju_tmp(ku,1)=999999
    endif
  enddo
  do k=i1_bdy_south(kkbdy),i2_bdy_south(kkbdy)
    kv=kv+1
    if ( SUM(vmaskutil(MAX(k-1,1):MIN(k+1,mx),jj_bdy_south(kkbdy))).ne.0 ) then
      nbiv_tmp(kv,1)=k
      nbjv_tmp(kv,1)=jj_bdy_south(kkbdy)
    else
      nbiv_tmp(kv,1)=999999
      nbjv_tmp(kv,1)=999999
    endif
  enddo
enddo

!--
! Check for locations defined two times (useful for complex boundaries)
! And remove them and adjust dimensions :

write(*,*) 'Check for double values or triple land points on grid T'

ndouble=0
do k1=1,mxbt-1
  if ( nbit_tmp(k1,1) .eq. 999999 .and. nbjt_tmp(k1,1) .eq. 999999 ) then ! triple land point
    ndouble=ndouble+1
  else
    do k2=k1+1,mxbt
      if ( nbit_tmp(k1,1).eq.nbit_tmp(k2,1) .and. nbjt_tmp(k1,1).eq.nbjt_tmp(k2,1) ) then
        ndouble=ndouble+1
        nbit_tmp(k1,1) = 999999
        nbjt_tmp(k1,1) = 999999
        write(*,*) 'adjusting ', ndouble, k1, k2, nbit_tmp(k2,1), nbjt_tmp(k2,1)
      endif
    enddo
  endif
enddo
ALLOCATE( nbit(mxbt-ndouble,myb), nbjt(mxbt-ndouble,myb) )
kt0=0
do k1=1,mxbt
  if ( nbit_tmp(k1,1).ne.999999 .and. nbjt_tmp(k1,1).ne.999999 ) then
    kt0=kt0+1
    nbit(kt0,:)=nbit_tmp(k1,:)
    nbjt(kt0,:)=nbjt_tmp(k1,:)
  endif
enddo
!DEALLOCATE( nbit_tmp, nbjt_tmp )
mxbt=mxbt-ndouble

write(*,*) 'Check for double values or triple land points on grid U'
ndouble=0
do k1=1,mxbu-1
  if ( nbiu_tmp(k1,1) .eq. 999999 .and. nbju_tmp(k1,1) .eq. 999999 ) then ! triple land point
    ndouble=ndouble+1
  else
    do k2=k1+1,mxbu
      if ( nbiu_tmp(k1,1).eq.nbiu_tmp(k2,1) .and. nbju_tmp(k1,1).eq.nbju_tmp(k2,1) ) then
        ndouble=ndouble+1
        nbiu_tmp(k1,1) = 999999
        nbju_tmp(k1,1) = 999999
        write(*,*) 'adjusting ', ndouble, k1, k2, nbiu_tmp(k2,1), nbju_tmp(k2,1)
      endif
    enddo
  endif
enddo
ALLOCATE( nbiu(mxbu-ndouble,myb), nbju(mxbu-ndouble,myb) )
ku0=0
do k1=1,mxbu
  if ( nbiu_tmp(k1,1).ne.999999 .and. nbju_tmp(k1,1).ne.999999 ) then
    ku0=ku0+1
    nbiu(ku0,:)=nbiu_tmp(k1,:)
    nbju(ku0,:)=nbju_tmp(k1,:)
  endif
enddo
!DEALLOCATE( nbiu_tmp, nbju_tmp )
mxbu=mxbu-ndouble

write(*,*) 'Check for double values or triple land points on grid V'
ndouble=0
do k1=1,mxbv-1
  if ( nbiv_tmp(k1,1) .eq. 999999 .and. nbjv_tmp(k1,1) .eq. 999999 ) then ! triple land point
    ndouble=ndouble+1
  else
    do k2=k1+1,mxbv
      if ( nbiv_tmp(k1,1).eq.nbiv_tmp(k2,1) .and. nbjv_tmp(k1,1).eq.nbjv_tmp(k2,1) ) then
        ndouble=ndouble+1
        nbiv_tmp(k1,1) = 999999
        nbjv_tmp(k1,1) = 999999
        write(*,*) 'adjusting ', ndouble, k1, k2, nbiv_tmp(k2,1), nbjv_tmp(k2,1)
      endif
    enddo
  endif
enddo
ALLOCATE( nbiv(mxbv-ndouble,myb), nbjv(mxbv-ndouble,myb) )
kv0=0
do k1=1,mxbv
  if ( nbiv_tmp(k1,1).ne.999999 .and. nbjv_tmp(k1,1).ne.999999 ) then
    kv0=kv0+1
    nbiv(kv0,:)=nbiv_tmp(k1,:)
    nbjv(kv0,:)=nbjv_tmp(k1,:)
  endif
enddo
!DEALLOCATE( nbiv_tmp, nbjv_tmp )
mxbv=mxbv-ndouble

write(*,*) 'Adjusted dimensions are:'
write(*,*) 'mxbt = ', mxbt
write(*,*) 'mxbu = ', mxbu
write(*,*) 'mxbv = ', mxbv

!---

ALLOCATE(  e1t_bdy(mxbt,myb), e2t_bdy(mxbt,myb)  )
ALLOCATE(  e1u_bdy(mxbu,myb), e2u_bdy(mxbu,myb)  )
ALLOCATE(  e1v_bdy(mxbv,myb), e2v_bdy(mxbv,myb)  )
ALLOCATE(  gphit_bdy(mxbt,myb), glamt_bdy(mxbt,myb)  )
ALLOCATE(  gphiu_bdy(mxbu,myb), glamu_bdy(mxbu,myb)  )
ALLOCATE(  gphiv_bdy(mxbv,myb), glamv_bdy(mxbv,myb)  )
ALLOCATE(  nbrt(mxbt,myb)  )
ALLOCATE(  nbru(mxbu,myb)  )
ALLOCATE(  nbrv(mxbv,myb)  )

!---

do k=1,mxbt
  write(*,*) '%% T ', k, nbit(k,1), nbjt(k,1)
  e1t_bdy  (k,1) = e1t  ( nbit(k,1), nbjt(k,1) )
  e2t_bdy  (k,1) = e2t  ( nbit(k,1), nbjt(k,1) )
  gphit_bdy(k,1) = gphit( nbit(k,1), nbjt(k,1) )
  glamt_bdy(k,1) = glamt( nbit(k,1), nbjt(k,1) )
enddo

do k=1,mxbu
  write(*,*) '%% U ', k, nbiu(k,1), nbju(k,1)
  e1u_bdy  (k,1) = e1u  ( nbiu(k,1), nbju(k,1) )
  e2u_bdy  (k,1) = e2u  ( nbiu(k,1), nbju(k,1) )
  gphiu_bdy(k,1) = gphiu( nbiu(k,1), nbju(k,1) )
  glamu_bdy(k,1) = glamu( nbiu(k,1), nbju(k,1) )
enddo

do k=1,mxbv
  write(*,*) '%% V ', k, nbiv(k,1), nbjv(k,1)
  e1v_bdy  (k,1) = e1v  ( nbiv(k,1), nbjv(k,1) )
  e2v_bdy  (k,1) = e2v  ( nbiv(k,1), nbjv(k,1) )
  gphiv_bdy(k,1) = gphiv( nbiv(k,1), nbjv(k,1) )
  glamv_bdy(k,1) = glamv( nbiv(k,1), nbjv(k,1) )
enddo

!- only one point along each bdy
nbrt(:,:)=1
nbru(:,:)=1
nbrv(:,:)=1

!=================================================================================
! 4- Write BDY coordinates in a netcdf file :
!=================================================================================

write(*,*) 'Writing ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY coordinates')

status = NF90_DEF_DIM(fidM,"xbv",mxbv,dimID_xbv)  ; call erreur(status,.TRUE.,"def_dimID_xbv")
status = NF90_DEF_DIM(fidM,"xbu",mxbu,dimID_xbu)  ; call erreur(status,.TRUE.,"def_dimID_xbu")
status = NF90_DEF_DIM(fidM,"xbt",mxbt,dimID_xbt)  ; call erreur(status,.TRUE.,"def_dimID_xbt")
status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)     ; call erreur(status,.TRUE.,"def_dimID_yb")

status = NF90_DEF_VAR(fidM,"e2v",NF90_FLOAT,(/dimID_xbv,dimID_yb/),e2v_ID)     ; call erreur(status,.TRUE.,"def_var_e2v_ID")
status = NF90_DEF_VAR(fidM,"e1v",NF90_FLOAT,(/dimID_xbv,dimID_yb/),e1v_ID)     ; call erreur(status,.TRUE.,"def_var_e1v_ID")
status = NF90_DEF_VAR(fidM,"gphiv",NF90_FLOAT,(/dimID_xbv,dimID_yb/),gphiv_ID) ; call erreur(status,.TRUE.,"def_var_gphiv_ID")
status = NF90_DEF_VAR(fidM,"glamv",NF90_FLOAT,(/dimID_xbv,dimID_yb/),glamv_ID) ; call erreur(status,.TRUE.,"def_var_glamv_ID")
status = NF90_DEF_VAR(fidM,"e2u",NF90_FLOAT,(/dimID_xbu,dimID_yb/),e2u_ID)     ; call erreur(status,.TRUE.,"def_var_e2u_ID")
status = NF90_DEF_VAR(fidM,"e1u",NF90_FLOAT,(/dimID_xbu,dimID_yb/),e1u_ID)     ; call erreur(status,.TRUE.,"def_var_e1u_ID")
status = NF90_DEF_VAR(fidM,"gphiu",NF90_FLOAT,(/dimID_xbu,dimID_yb/),gphiu_ID) ; call erreur(status,.TRUE.,"def_var_gphiu_ID")
status = NF90_DEF_VAR(fidM,"glamu",NF90_FLOAT,(/dimID_xbu,dimID_yb/),glamu_ID) ; call erreur(status,.TRUE.,"def_var_glamu_ID")
status = NF90_DEF_VAR(fidM,"e2t",NF90_FLOAT,(/dimID_xbt,dimID_yb/),e2t_ID)     ; call erreur(status,.TRUE.,"def_var_e2t_ID")
status = NF90_DEF_VAR(fidM,"e1t",NF90_FLOAT,(/dimID_xbt,dimID_yb/),e1t_ID)     ; call erreur(status,.TRUE.,"def_var_e1t_ID")
status = NF90_DEF_VAR(fidM,"gphit",NF90_FLOAT,(/dimID_xbt,dimID_yb/),gphit_ID) ; call erreur(status,.TRUE.,"def_var_gphit_ID")
status = NF90_DEF_VAR(fidM,"glamt",NF90_FLOAT,(/dimID_xbt,dimID_yb/),glamt_ID) ; call erreur(status,.TRUE.,"def_var_glamt_ID")
status = NF90_DEF_VAR(fidM,"nbrv",NF90_INT,(/dimID_xbv,dimID_yb/),nbrv_ID)     ; call erreur(status,.TRUE.,"def_var_nbrv_ID")
status = NF90_DEF_VAR(fidM,"nbjv",NF90_INT,(/dimID_xbv,dimID_yb/),nbjv_ID)     ; call erreur(status,.TRUE.,"def_var_nbjv_ID")
status = NF90_DEF_VAR(fidM,"nbiv",NF90_INT,(/dimID_xbv,dimID_yb/),nbiv_ID)     ; call erreur(status,.TRUE.,"def_var_nbiv_ID")
status = NF90_DEF_VAR(fidM,"nbru",NF90_INT,(/dimID_xbu,dimID_yb/),nbru_ID)     ; call erreur(status,.TRUE.,"def_var_nbru_ID")
status = NF90_DEF_VAR(fidM,"nbju",NF90_INT,(/dimID_xbu,dimID_yb/),nbju_ID)     ; call erreur(status,.TRUE.,"def_var_nbju_ID")
status = NF90_DEF_VAR(fidM,"nbiu",NF90_INT,(/dimID_xbu,dimID_yb/),nbiu_ID)     ; call erreur(status,.TRUE.,"def_var_nbiu_ID")
status = NF90_DEF_VAR(fidM,"nbrt",NF90_INT,(/dimID_xbt,dimID_yb/),nbrt_ID)     ; call erreur(status,.TRUE.,"def_var_nbrt_ID")
status = NF90_DEF_VAR(fidM,"nbjt",NF90_INT,(/dimID_xbt,dimID_yb/),nbjt_ID)     ; call erreur(status,.TRUE.,"def_var_nbjt_ID")
status = NF90_DEF_VAR(fidM,"nbit",NF90_INT,(/dimID_xbt,dimID_yb/),nbit_ID)     ; call erreur(status,.TRUE.,"def_var_nbit_ID")

status = NF90_PUT_ATT(fidM,e2v_ID,"units","meters")          ; call erreur(status,.TRUE.,"put_att_e2v_ID")
status = NF90_PUT_ATT(fidM,e1v_ID,"units","meters")          ; call erreur(status,.TRUE.,"put_att_e1v_ID")
status = NF90_PUT_ATT(fidM,gphiv_ID,"units","degrees_north") ; call erreur(status,.TRUE.,"put_att_gphiv_ID")
status = NF90_PUT_ATT(fidM,glamv_ID,"units","degrees_east")  ; call erreur(status,.TRUE.,"put_att_glamv_ID")
status = NF90_PUT_ATT(fidM,e2u_ID,"units","meters")          ; call erreur(status,.TRUE.,"put_att_e2u_ID")
status = NF90_PUT_ATT(fidM,e1u_ID,"units","meters")          ; call erreur(status,.TRUE.,"put_att_e1u_ID")
status = NF90_PUT_ATT(fidM,gphiu_ID,"units","degrees_north") ; call erreur(status,.TRUE.,"put_att_gphiu_ID")
status = NF90_PUT_ATT(fidM,glamu_ID,"units","degrees_east")  ; call erreur(status,.TRUE.,"put_att_glamu_ID")
status = NF90_PUT_ATT(fidM,e2t_ID,"units","meters")          ; call erreur(status,.TRUE.,"put_att_e2t_ID")
status = NF90_PUT_ATT(fidM,e1t_ID,"units","meters")          ; call erreur(status,.TRUE.,"put_att_e1t_ID")
status = NF90_PUT_ATT(fidM,gphit_ID,"units","degrees_north") ; call erreur(status,.TRUE.,"put_att_gphit_ID")
status = NF90_PUT_ATT(fidM,glamt_ID,"units","degrees_east")  ; call erreur(status,.TRUE.,"put_att_glamt_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using build_coordinates_bdy.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")

status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"fin_definition") 

status = NF90_PUT_VAR(fidM,e2v_ID,e2v_bdy)     ; call erreur(status,.TRUE.,"var_e2v_ID")
status = NF90_PUT_VAR(fidM,e1v_ID,e1v_bdy)     ; call erreur(status,.TRUE.,"var_e1v_ID")
status = NF90_PUT_VAR(fidM,gphiv_ID,gphiv_bdy) ; call erreur(status,.TRUE.,"var_gphiv_ID")
status = NF90_PUT_VAR(fidM,glamv_ID,glamv_bdy) ; call erreur(status,.TRUE.,"var_glamv_ID")
status = NF90_PUT_VAR(fidM,e2u_ID,e2u_bdy)     ; call erreur(status,.TRUE.,"var_e2u_ID")
status = NF90_PUT_VAR(fidM,e1u_ID,e1u_bdy)     ; call erreur(status,.TRUE.,"var_e1u_ID")
status = NF90_PUT_VAR(fidM,gphiu_ID,gphiu_bdy) ; call erreur(status,.TRUE.,"var_gphiu_ID")
status = NF90_PUT_VAR(fidM,glamu_ID,glamu_bdy) ; call erreur(status,.TRUE.,"var_glamu_ID")
status = NF90_PUT_VAR(fidM,e2t_ID,e2t_bdy)     ; call erreur(status,.TRUE.,"var_e2t_ID")
status = NF90_PUT_VAR(fidM,e1t_ID,e1t_bdy)     ; call erreur(status,.TRUE.,"var_e1t_ID")
status = NF90_PUT_VAR(fidM,gphit_ID,gphit_bdy) ; call erreur(status,.TRUE.,"var_gphit_ID")
status = NF90_PUT_VAR(fidM,glamt_ID,glamt_bdy) ; call erreur(status,.TRUE.,"var_glamt_ID")
status = NF90_PUT_VAR(fidM,nbrv_ID,nbrv)       ; call erreur(status,.TRUE.,"var_nbrv_ID")
status = NF90_PUT_VAR(fidM,nbjv_ID,nbjv)       ; call erreur(status,.TRUE.,"var_nbjv_ID")
status = NF90_PUT_VAR(fidM,nbiv_ID,nbiv)       ; call erreur(status,.TRUE.,"var_nbiv_ID")
status = NF90_PUT_VAR(fidM,nbru_ID,nbru)       ; call erreur(status,.TRUE.,"var_nbru_ID")
status = NF90_PUT_VAR(fidM,nbju_ID,nbju)       ; call erreur(status,.TRUE.,"var_nbju_ID")
status = NF90_PUT_VAR(fidM,nbiu_ID,nbiu)       ; call erreur(status,.TRUE.,"var_nbiu_ID")
status = NF90_PUT_VAR(fidM,nbrt_ID,nbrt)       ; call erreur(status,.TRUE.,"var_nbrt_ID")
status = NF90_PUT_VAR(fidM,nbjt_ID,nbjt)       ; call erreur(status,.TRUE.,"var_nbjt_ID")
status = NF90_PUT_VAR(fidM,nbit_ID,nbit)       ; call erreur(status,.TRUE.,"var_nbit_ID")

status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close new coordinate file")         

end program modif

!==============================================
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
    WRITE(*,*) 'ERROR: ', iret
    message=NF90_STRERROR(iret)
    WRITE(*,*) 'WHICH MEANS:',TRIM(message)
    IF ( lstop ) STOP
  ENDIF
  !
END SUBROUTINE erreur
