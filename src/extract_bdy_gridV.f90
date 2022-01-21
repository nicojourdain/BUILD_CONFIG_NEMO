program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Used to build netcdf gridV file for BDY
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read input file dimensions in first existing file for specified time window
! 3- Process all gridV files over specified period
!
! history : - Jan. 2015: initial version (N. Jourdain, CNRS-LGGE)
!           - Feb. 2017: version with namelist (N. Jourdain, CNRS-IGE)
!           - Jan. 2022: new convention for variable names (PAR/CHLD/CHLD)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /bdy_data/ nn_yeari, nn_yearf, data_dir, data_prefix, nn_bdy_eosmatch, &
&                   data_suffix_T, data_suffix_S, data_suffix_U, data_suffix_V, &
&                   data_suffix_ssh, data_suffix_ice, file_data_mask,           &
&                   file_data_zgr, file_data_hgr, sep1, sep2

CHARACTER(LEN=50)                    :: config, sep1, sep2
CHARACTER(LEN=150)                   :: config_dir, data_dir, data_prefix, data_suffix_T, data_suffix_S, &
&                                       data_suffix_U, data_suffix_V, data_suffix_ssh, data_suffix_ice,  &
&                                       file_data_mask, file_data_zgr, file_data_hgr
INTEGER                              :: nn_yeari, nn_yearf, nn_bdy_eosmatch

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbv,                       &
&                                       myb, mxbv, glamv_ID, gphiv_ID, e1v_ID, e2v_ID, e1v_PAR_ID,   &
&                                       nbiv_ID, nbjv_ID, nbrv_ID, mtime, dimID_x, dimID_y,          &
&                                       mlon, mlat, mdepthv, kday, kmonth, kyear, kbdy, nfmt, fidN,  &
&                                       kt, kz,     lon_ID, lat_ID, depthv_ID, vobtcrty_ID, fidM,    &
&                                       vomecrty_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidV,  &
&                                       dimID_time_counter, dimID_depthv, time_ID, dimID_time, fidS, &
&                                       i, j, k, l, fidC, imin_EXT, jmin_EXT, iPAR, jPAR,      &
&                                       depth_ID, ai, aj, bi, bj, kfmt, fidMSKIN, vmask_PAR_ID,      &
&                                       e3v_PAR_ID, fidZGRIN, fidHGRIN, e3v_ID
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_gridV, file_bdy_gridV2d,  &
&                                       file_in_coord_CHLD, command_str, file_bdy_gridV3d
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: list_fmt
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbiv, nbjv, nbrv
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: vmask_PAR
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: e1v_PAR, glamv, gphiv, e1v, e2v, nav_lon, nav_lat, nav_lon_bdy, nav_lat_bdy
REAL*4,ALLOCATABLE,DIMENSION(:,:,:)  :: e3v_PAR, vobtcrty, vobtcrty_bdy
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:):: vomecrty, vomecrty_bdy, e3v
REAL*4,ALLOCATABLE,DIMENSION(:)      :: depthv
REAL*8,ALLOCATABLE,DIMENSION(:)      :: time
REAL*4                               :: thic
LOGICAL                              :: existfile, ln_vvl

!=================================================================================
!- 0- Initialiartions
!=================================================================================

write(*,*) 'Reading namelist parameters'

! Default values (replaced with namelist values if specified):
config_dir        = '.'
nn_bdy_eosmatch   =   1

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=bdy_data)
CLOSE(1)

!- bdy coordinates :
write(file_coord,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/coordinates_bdy_',a,'.nc')

!- name of child domain coordinates (input) :
write(file_in_coord_CHLD,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coordinates_',a,'.nc')

!=================================================================================
! 1- Read information on grids
!=================================================================================

!- Read global attributes of coordinate file to get grid correspondance :
!       i_EXT = ai * i_PAR + bi
!       j_EXT = aj * j_PAR + bj

write(*,*) 'Reading parameters for grid correspondence in the attributes of ', TRIM(file_in_coord_CHLD)
status = NF90_OPEN(TRIM(file_in_coord_CHLD),0,fidC); call erreur(status,.TRUE.,"read global grid coord input")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "ai", ai); call erreur(status,.TRUE.,"read att1")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "bi", bi); call erreur(status,.TRUE.,"read att2")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "aj", aj); call erreur(status,.TRUE.,"read att3")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "bj", bj); call erreur(status,.TRUE.,"read att4")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "imin_extraction", imin_EXT); call erreur(status,.TRUE.,"read att5")
status = NF90_GET_ATT(fidC, NF90_GLOBAL, "jmin_extraction", jmin_EXT); call erreur(status,.TRUE.,"read att6")
status = NF90_CLOSE(fidC)                         ; call erreur(status,.TRUE.,"end read fidC")

!- Read BDY coordinates 

write(*,*) 'Reading BDY coordinates in ', TRIM(file_coord)
status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read bdy coordinate file") 

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbv",dimID_xbv) ; call erreur(status,.TRUE.,"inq_dimID_xbv")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbv,len=mxbv) ; call erreur(status,.TRUE.,"inq_dim_xbv")

ALLOCATE(  glamv(mxbv,myb)  ) 
ALLOCATE(  gphiv(mxbv,myb)  ) 
ALLOCATE(  e1v(mxbv,myb)  ) 
ALLOCATE(  e2v(mxbv,myb)  ) 
ALLOCATE(  nbiv(mxbv,myb)  ) 
ALLOCATE(  nbjv(mxbv,myb)  ) 
ALLOCATE(  nbrv(mxbv,myb)  ) 

status = NF90_INQ_VARID(fidCOORD,"glamv",glamv_ID) ; call erreur(status,.TRUE.,"inq_glamv_ID")
status = NF90_INQ_VARID(fidCOORD,"gphiv",gphiv_ID) ; call erreur(status,.TRUE.,"inq_gphiv_ID")
status = NF90_INQ_VARID(fidCOORD,"e1v",e1v_ID)     ; call erreur(status,.TRUE.,"inq_e1v_ID")
status = NF90_INQ_VARID(fidCOORD,"e2v",e2v_ID)     ; call erreur(status,.TRUE.,"inq_e2v_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiv",nbiv_ID)   ; call erreur(status,.TRUE.,"inq_nbiv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjv",nbjv_ID)   ; call erreur(status,.TRUE.,"inq_nbjv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrv",nbrv_ID)   ; call erreur(status,.TRUE.,"inq_nbrv_ID")

status = NF90_GET_VAR(fidCOORD,glamv_ID,glamv) ; call erreur(status,.TRUE.,"getvar_glamv")
status = NF90_GET_VAR(fidCOORD,gphiv_ID,gphiv) ; call erreur(status,.TRUE.,"getvar_gphiv")
status = NF90_GET_VAR(fidCOORD,e1v_ID,e1v)     ; call erreur(status,.TRUE.,"getvar_e1v")
status = NF90_GET_VAR(fidCOORD,e2v_ID,e2v)     ; call erreur(status,.TRUE.,"getvar_e2v")
status = NF90_GET_VAR(fidCOORD,nbiv_ID,nbiv)   ; call erreur(status,.TRUE.,"getvar_nbiv")
status = NF90_GET_VAR(fidCOORD,nbjv_ID,nbjv)   ; call erreur(status,.TRUE.,"getvar_nbjv")
status = NF90_GET_VAR(fidCOORD,nbrv_ID,nbrv)   ; call erreur(status,.TRUE.,"getvar_nbrv")

status = NF90_CLOSE(fidCOORD) ; call erreur(status,.TRUE.,"close coordinate file")     


!=================================================================================
! 2- Read input file dimensions in first existing file for specified time window
!=================================================================================

!- accepted input format :
191 FORMAT(a,'/',i4.4,'/',a,i4.4,a,i2.2,a,i2.2,a,'.nc')  ! <data_dir>/YYYY/<data_prefix>YYYY<sep1>MM<sep2>DD<data_suffix>.nc  
192 FORMAT(a,'/',i4.4,'/',a,i4.4,a,i2.2,a,'.nc')         ! <data_dir>/YYYY/<data_prefix>YYYY<sep1>MM<data_suffix>.nc  
193 FORMAT(a,'/',a,i4.4,a,i2.2,a,i2.2,a,'.nc')           ! <data_dir>/<data_prefix>YYYY<sep1>MM<sep2>DD<data_suffix>.nc  
194 FORMAT(a,'/',a,i4.4,a,i2.2,a,'.nc')                  ! <data_dir>/<data_prefix>YYYY<sep1>MM<data_suffix>.nc 

ALLOCATE(list_fmt(4))
list_fmt=(/191,192,193,194/)

kyear=nn_yeari
kmonth=1
DO kday=1,31

  do kfmt=1,size(list_fmt)
     nfmt=list_fmt(kfmt)
     SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridV,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_V)
        CASE(192)
          write(file_in_gridV,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_V) 
        CASE(193)
          write(file_in_gridV,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_V)
        CASE(194)
          write(file_in_gridV,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_V)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     write(*,*) 'Looking for existence of ', TRIM(file_in_gridV)
     inquire(file=file_in_gridV, exist=existfile)
     if ( existfile ) then; write(*,*) 'BINGO !'; exit ; endif
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading T,S input dimensions in ', TRIM(file_in_gridV)
    status = NF90_OPEN(TRIM(file_in_gridV),0,fidV)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidV,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidV,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidV,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")
    status = NF90_INQ_DIMID(fidV,"z",dimID_depthv)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidV,"depth",dimID_depthv)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidV,"depthv",dimID_depthv)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidV,"deptht",dimID_depthv)
    call erreur(status,.TRUE.,"inq_dimID_depthv")

    status = NF90_INQUIRE_DIMENSION(fidV,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidV,dimID_x,len=mlon)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidV,dimID_y,len=mlat)         ; call erreur(status,.TRUE.,"inq_dim_y")
    status = NF90_INQUIRE_DIMENSION(fidV,dimID_depthv,len=mdepthv) ; call erreur(status,.TRUE.,"inq_dim_depthv")

    ALLOCATE( nav_lon(mlon,mlat), nav_lat(mlon,mlat) )
    ALLOCATE( depthv(mdepthv) )

    status = NF90_INQ_VARID(fidV,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidV,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
    status = NF90_INQ_VARID(fidV,"depth",depthv_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"depthv",depthv_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"deptht",depthv_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"nav_lev",depthv_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"z",depthv_ID)
    call erreur(status,.TRUE.,"inq_depthv_ID")
        
    status = NF90_GET_VAR(fidV,lon_ID,nav_lon)    ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidV,lat_ID,nav_lat)    ; call erreur(status,.TRUE.,"getvar_lat")
    status = NF90_GET_VAR(fidV,depthv_ID,depthv)  ; call erreur(status,.TRUE.,"getvar_depthv")

    status = NF90_CLOSE(fidV)                     ; call erreur(status,.TRUE.,"fin_lecture")

    exit

  ELSEIF ( kday .eq. 31 ) THEN

    write(*,*) 'No velocity file found for first month of ', kyear
    write(*,*) '        >>>>>>>>>>>>>>>>>> stop !!'
    stop

  ENDIF

ENDDO

!- Read vmask and e3v in parent grid file:
status = NF90_OPEN(TRIM(file_data_mask),0,fidMSKIN);    call erreur(status,.TRUE.,"read mask input") 
ALLOCATE(  vmask_PAR(mlon,mlat,mdepthv)  ) 
status = NF90_INQ_VARID(fidMSKIN,"vmask",vmask_PAR_ID); call erreur(status,.TRUE.,"inq_vmask_PAR_ID")
status = NF90_GET_VAR(fidMSKIN,vmask_PAR_ID,vmask_PAR); call erreur(status,.TRUE.,"getvar_vmask_PAR")
status = NF90_CLOSE(fidMSKIN);                          call erreur(status,.TRUE.,"end read mask_PAR")

!- Read e1v and e3v in parent grid file:
status = NF90_OPEN(TRIM(file_data_zgr),0,fidZGRIN);     call erreur(status,.TRUE.,"read mask input")
ALLOCATE(  e3v_PAR(mlon,mlat,mdepthv)  )
status = NF90_INQ_VARID(fidZGRIN,"e3v",e3v_PAR_ID) 
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidZGRIN,"e3v_0",e3v_PAR_ID)
call erreur(status,.TRUE.,"inq_e3v_PAR_ID")
status = NF90_GET_VAR(fidZGRIN,e3v_PAR_ID,e3v_PAR);     call erreur(status,.TRUE.,"getvar_e3v_PAR")
status = NF90_CLOSE(fidZGRIN);                          call erreur(status,.TRUE.,"end read mask_PAR")
!-
status = NF90_OPEN(TRIM(file_data_hgr),0,fidHGRIN);     call erreur(status,.TRUE.,"read mask input")
ALLOCATE(  e1v_PAR(mlon,mlat)  )
status = NF90_INQ_VARID(fidHGRIN,"e1v",e1v_PAR_ID);     call erreur(status,.TRUE.,"inq_e1v_PAR_ID")
status = NF90_GET_VAR(fidHGRIN,e1v_PAR_ID,e1v_PAR);     call erreur(status,.TRUE.,"getvar_e1v_PAR")
status = NF90_CLOSE(fidHGRIN);                          call erreur(status,.TRUE.,"end read mask_PAR")

!--

ALLOCATE( nav_lon_bdy(mxbv,1), nav_lat_bdy(mxbv,1) )

do kbdy=1,mxbv
  iPAR=NINT(FLOAT(nbiv(kbdy,1)+imin_EXT-1-bi)/ai)
  jPAR=NINT(FLOAT(nbjv(kbdy,1)+jmin_EXT-1-bj)/aj)
  nav_lon_bdy(kbdy,1) = nav_lon(iPAR,jPAR)
  nav_lat_bdy(kbdy,1) = nav_lat(iPAR,jPAR)
enddo

!--

write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir -pv ',a,'/BDY')
CALL system(TRIM(command_str))

!=================================================================================
! 3- Process all gridV files over specified period
!=================================================================================

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridV,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_V)
        CASE(192)
          write(file_in_gridV,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_V) 
        CASE(193)
          write(file_in_gridV,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_V)
        CASE(194)
          write(file_in_gridV,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_V)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_gridV, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          401 FORMAT(a,'/BDY/bdyV_u2d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridV2d,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
          501 FORMAT(a,'/BDY/bdyV_u3d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridV3d,501) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 ) then
          402 FORMAT(a,'/BDY/bdyV_u2d_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridV2d,402) TRIM(config_dir), kyear, kmonth, TRIM(config)
          502 FORMAT(a,'/BDY/bdyV_u3d_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridV3d,502) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_bdy_gridV3d and file_bdy_gridV3d  >>>> stop'
          stop
        endif

        ALLOCATE( vomecrty(mlon,mlat,mdepthv,mtime)  )
        ALLOCATE( vobtcrty(mlon,mlat,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read input velocity : 

        write(*,*) 'Reading velocities in ', TRIM(file_in_gridV)
        
        status = NF90_OPEN(TRIM(file_in_gridV),0,fidV)                ; call erreur(status,.TRUE.,"read EXT TS") 

        status = NF90_INQ_VARID(fidV,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_GET_VAR(fidV,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_ATT(fidV,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidV,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")

        status = NF90_INQ_VARID(fidV,"vomecrty",vomecrty_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"voce",vomecrty_ID)
        call erreur(status,.TRUE.,"None existing velocity, or unknown variable name")
        status = NF90_GET_VAR(fidV,vomecrty_ID,vomecrty)              ; call erreur(status,.TRUE.,"getvar_vomecrty")

        status = NF90_INQ_VARID(fidV,"e3v",e3v_ID)
        if ( status .ne. 0 ) then
          write(*,*) 'No e3v found in grid_V file => assuming not in vvl mode (i.e. constant e3v)'
          ln_vvl = .false.
        else
          ALLOCATE( e3v(mlon,mlat,mdepthv,mtime)  )
          write(*,*) 'Barotropic component calculated from time-dependent e3v'
          status = NF90_GET_VAR(fidV,e3v_ID,e3v) ; call erreur(status,.TRUE.,"getvar_e3v")
          ln_vvl = .true.
        endif

        status = NF90_CLOSE(fidV)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Remove possible NaNs :
 
        do i=1,mlon
        do j=1,mlat
        do k=1,mdepthv
        do l=1,mtime
          if ( .not. abs(vomecrty(i,j,k,l)) .lt. 99.0 ) then
            vomecrty(i,j,k,l) = 0.0
          endif
        enddo
        enddo
        enddo
        enddo

        !---------------------------------------
        ! Calculate the barotropic component :
  
        do l=1,mtime
          if ( ln_vvl ) e3v_PAR(:,:,:) = e3v(:,:,:,l)
          do i=1,mlon
          do j=1,mlat
            vobtcrty(i,j,l) = 0.0
            thic=0.0
            do k=1,mdepthv
              if ( vmask_PAR(i,j,k) .eq. 1 ) then
                vobtcrty(i,j,l) = vobtcrty(i,j,l) + vomecrty(i,j,k,l) * e3v_PAR(i,j,k) * vmask_PAR(i,j,k)
                thic            = thic            +                     e3v_PAR(i,j,k) * vmask_PAR(i,j,k)
              endif
            enddo
            if ( thic .gt. 1.e-3 ) then
              ! barotropic component :
              vobtcrty(i,j,l) = vobtcrty(i,j,l) / thic
            else
              vobtcrty(i,j,l) = 0.0
            endif
            ! baroclinic component :
            do k=1,mdepthv
              vomecrty(i,j,k,l) = ( vomecrty(i,j,k,l) - vobtcrty(i,j,l) ) * vmask_PAR(i,j,k)
            enddo
          enddo
          enddo
        enddo

        !---------------------------------------
        ! Fill values on bdyV :
        ! (we ensure transport conservation in case grids are 
        ! slightly stretched with respect to each other)     
 
        ALLOCATE( vomecrty_bdy(mxbv,1,mdepthv,mtime)  )
        ALLOCATE( vobtcrty_bdy(mxbv,1,mtime)  )

        do kbdy=1,mxbv
          iPAR=NINT(FLOAT(nbiv(kbdy,1)+imin_EXT-1-bi)/ai)
          jPAR=NINT(FLOAT(nbjv(kbdy,1)+jmin_EXT-1-bj)/aj)
          write(*,*) kbdy, iPAR, jPAR, nbiv(kbdy,1), nbjv(kbdy,1)
          do kt=1,mtime
          do kz=1,mdepthv
            vomecrty_bdy(kbdy,1,kz,kt) = vomecrty( iPAR, jPAR, kz, kt ) * e1v_PAR( iPAR, jPAR ) / (e1v(kbdy,1)*ai)
            vobtcrty_bdy(kbdy,1,   kt) = vobtcrty( iPAR, jPAR,     kt ) * e1v_PAR( iPAR, jPAR ) / (e1v(kbdy,1)*ai)
          enddo
          enddo
        enddo

        !--------------------------------------
        ! Write BDY netcdf file for barotropic velocity :

        write(*,*) 'Creating ', TRIM(file_bdy_gridV2d)
        status = NF90_CREATE(TRIM(file_bdy_gridV2d),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidM,"xbV",mxbV,dimID_xbV)                             ; call erreur(status,.TRUE.,"def_dimID_xbV")

        status = NF90_DEF_VAR(fidM,"vobtcrty",NF90_FLOAT,(/dimID_xbV,dimID_yb,dimID_time_counter/),vobtcrty_ID)
        call erreur(status,.TRUE.,"def_var_vobtcrty_ID")
        status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_xbV,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_xbV,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"nbrdta",NF90_INT,(/dimID_xbV,dimID_yb/),nbrv_ID)
        call erreur(status,.TRUE.,"def_var_nbrv_ID")
        status = NF90_DEF_VAR(fidM,"nbjdta",NF90_INT,(/dimID_xbV,dimID_yb/),nbjv_ID)
        call erreur(status,.TRUE.,"def_var_nbjv_ID")
        status = NF90_DEF_VAR(fidM,"nbidta",NF90_INT,(/dimID_xbV,dimID_yb/),nbiv_ID)
        call erreur(status,.TRUE.,"def_var_nbiv_ID")
       
        status = NF90_PUT_ATT(fidM,vobtcrty_ID,"long_name","Barotropic velocity")   ; call erreur(status,.TRUE.,"put_att_vobtcrty_ID")
        status = NF90_PUT_ATT(fidM,vobtcrty_ID,"units","m/s")                       ; call erreur(status,.TRUE.,"put_att_vobtcrty_ID")
        status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,nbrv_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidM,nbrv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidM,nbjv_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidM,nbjv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidM,nbiv_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        status = NF90_PUT_ATT(fidM,nbiv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bdy_gridV.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidM,vobtcrty_ID,vobtcrty_bdy) ; call erreur(status,.TRUE.,"var_vobtcrty_ID")
        status = NF90_PUT_VAR(fidM,nav_lat_ID,nav_lat_bdy)   ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidM,nav_lon_ID,nav_lon_bdy)   ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidM,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidM,nbrv_ID,nbrv)             ; call erreur(status,.TRUE.,"var_nbrv_ID")
        status = NF90_PUT_VAR(fidM,nbjv_ID,nbjv)             ; call erreur(status,.TRUE.,"var_nbjv_ID")
        status = NF90_PUT_VAR(fidM,nbiv_ID,nbiv)             ; call erreur(status,.TRUE.,"var_nbiv_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close BDY file")

        !--------------------------------------
        ! Write BDY netcdf file for baroclinic velocity :

        write(*,*) 'Creating ', TRIM(file_bdy_gridV3d)
        status = NF90_CREATE(TRIM(file_bdy_gridV3d),NF90_NOCLOBBER,fidN) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidN,"depthv",mdepthv,dimID_depthv)                    ; call erreur(status,.TRUE.,"def_dimID_depthv")
        status = NF90_DEF_DIM(fidN,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidN,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidN,"xbV",mxbV,dimID_xbV)                             ; call erreur(status,.TRUE.,"def_dimID_xbV")

        status = NF90_DEF_VAR(fidN,"vomecrty",NF90_FLOAT,(/dimID_xbV,dimID_yb,dimID_depthv,dimID_time_counter/),vomecrty_ID)
        call erreur(status,.TRUE.,"def_var_vomecrty_ID")
        status = NF90_DEF_VAR(fidN,"depthv",NF90_FLOAT,(/dimID_depthv/),depthv_ID)
        call erreur(status,.TRUE.,"def_var_depthv_ID")
        status = NF90_DEF_VAR(fidN,"nav_lat",NF90_FLOAT,(/dimID_xbV,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidN,"nav_lon",NF90_FLOAT,(/dimID_xbV,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidN,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidN,"nbrdta",NF90_INT,(/dimID_xbV,dimID_yb/),nbrv_ID)
        call erreur(status,.TRUE.,"def_var_nbrv_ID")
        status = NF90_DEF_VAR(fidN,"nbjdta",NF90_INT,(/dimID_xbV,dimID_yb/),nbjv_ID)
        call erreur(status,.TRUE.,"def_var_nbjv_ID")
        status = NF90_DEF_VAR(fidN,"nbidta",NF90_INT,(/dimID_xbV,dimID_yb/),nbiv_ID)
        call erreur(status,.TRUE.,"def_var_nbiv_ID")
       
        status = NF90_PUT_ATT(fidN,vomecrty_ID,"long_name","Baroclinic Velocity")   ; call erreur(status,.TRUE.,"put_att_vomecrty_ID")
        status = NF90_PUT_ATT(fidN,vomecrty_ID,"units","m/s")                       ; call erreur(status,.TRUE.,"put_att_vomecrty_ID")
        status = NF90_PUT_ATT(fidN,depthv_ID,"long_name","Vertical T levels")       ; call erreur(status,.TRUE.,"put_att_depthv_ID")
        status = NF90_PUT_ATT(fidN,depthv_ID,"units","m")                           ; call erreur(status,.TRUE.,"put_att_depthv_ID")
        status = NF90_PUT_ATT(fidN,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidN,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidN,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidN,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidN,nbrv_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidN,nbrv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidN,nbjv_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidN,nbjv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidN,nbiv_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        status = NF90_PUT_ATT(fidN,nbiv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        
        status = NF90_PUT_ATT(fidN,NF90_GLOBAL,"history","Created using extract_bdy_gridV.f90")
        status = NF90_PUT_ATT(fidN,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidN) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidN,vomecrty_ID,vomecrty_bdy) ; call erreur(status,.TRUE.,"var_vomecrty_ID")
        status = NF90_PUT_VAR(fidN,depthv_ID,depthv)         ; call erreur(status,.TRUE.,"var_depthv_ID")
        status = NF90_PUT_VAR(fidN,nav_lat_ID,nav_lat_bdy)   ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidN,nav_lon_ID,nav_lon_bdy)   ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidN,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidN,nbrv_ID,nbrv)             ; call erreur(status,.TRUE.,"var_nbrv_ID")
        status = NF90_PUT_VAR(fidN,nbjv_ID,nbjv)             ; call erreur(status,.TRUE.,"var_nbjv_ID")
        status = NF90_PUT_VAR(fidN,nbiv_ID,nbiv)             ; call erreur(status,.TRUE.,"var_nbiv_ID")
        
        status = NF90_CLOSE(fidN) ; call erreur(status,.TRUE.,"close BDY file")
        
        !--       
        
        !--       
        DEALLOCATE( vomecrty, vobtcrty, time )
        DEALLOCATE( vomecrty_bdy, vobtcrty_bdy )
        IF (ln_vvl) DEALLOCATE(e3v)

        !--
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          write(*,*) 'Looking for next existing day in this month/year'
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 ) then
          write(*,*) 'Only one file per month => switching to next month'
          exit
        else
          write(*,*) 'Keep in mind to decide how to end the main loop according to new file format  >>>> stop'
          stop
        endif

      ENDIF

    ENDDO !-kday
  ENDDO !- kmonth
ENDDO !- kyear

end program modif

!=============================================

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
