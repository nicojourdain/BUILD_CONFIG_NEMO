program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to build netcdf gridU file for BDY
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read input file dimensions in first existing file for specified time window
! 3- Process all gridU files over specified period
!
! history : - Feb. 2017: version with namelist (N. Jourdain)
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

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbt, dimID_xbu, dimID_xbv, &
&                                       myb, mxbu, glamu_ID, gphiu_ID, e1u_ID, e2u_ID, e2u_GLO_ID,   &
&                                       nbiu_ID, nbju_ID, nbru_ID, mtime, dimID_x, dimID_y,          &
&                                       mlon, mlat, mdepthu, kday, kmonth, kyear, kbdy, nfmt, fidN,  &
&                                       kt, kz,     lon_ID, lat_ID, depthu_ID, vobtcrtx_ID, fidM,    &
&                                       vozocrtx_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidU,  &
&                                       dimID_time_counter, dimID_depthu, time_ID, dimID_time, fidS, &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO,      &
&                                       depth_ID, ai, aj, bi, bj, kfmt, fidMSKIN, umask_GLO_ID,      &
&                                       e3u_GLO_ID, fidZGRIN, fidHGRIN, e3u_ID
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_gridU, file_bdy_gridU2d,  &
&                                       file_in_coord_REG, command_str, file_bdy_gridU3d
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: list_fmt
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbiu, nbju, nbru
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: umask_GLO
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: e2u_GLO, glamu, gphiu, e1u, e2u, nav_lon, nav_lat, nav_lon_bdy, nav_lat_bdy
REAL*4,ALLOCATABLE,DIMENSION(:,:,:)  :: e3u_GLO, vobtcrtx, vobtcrtx_bdy
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:):: e3u, vozocrtx, vozocrtx_bdy
REAL*4,ALLOCATABLE,DIMENSION(:)      :: depthu
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

!- name of regional coordinates (input) :
write(file_in_coord_REG,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coordinates_',a,'.nc')

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

!- Read BDY coordinates 

write(*,*) 'Reading BDY coordinates in ', TRIM(file_coord)
status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read bdy coordinate file") 

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbu",dimID_xbu) ; call erreur(status,.TRUE.,"inq_dimID_xbu")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbu,len=mxbu) ; call erreur(status,.TRUE.,"inq_dim_xbu")

ALLOCATE(  glamu(mxbu,myb)  ) 
ALLOCATE(  gphiu(mxbu,myb)  ) 
ALLOCATE(  e1u(mxbu,myb)  ) 
ALLOCATE(  e2u(mxbu,myb)  ) 
ALLOCATE(  nbiu(mxbu,myb)  ) 
ALLOCATE(  nbju(mxbu,myb)  ) 
ALLOCATE(  nbru(mxbu,myb)  ) 

status = NF90_INQ_VARID(fidCOORD,"glamu",glamu_ID) ; call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidCOORD,"gphiu",gphiu_ID) ; call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidCOORD,"e1u",e1u_ID)     ; call erreur(status,.TRUE.,"inq_e1u_ID")
status = NF90_INQ_VARID(fidCOORD,"e2u",e2u_ID)     ; call erreur(status,.TRUE.,"inq_e2u_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiu",nbiu_ID)   ; call erreur(status,.TRUE.,"inq_nbiu_ID")
status = NF90_INQ_VARID(fidCOORD,"nbju",nbju_ID)   ; call erreur(status,.TRUE.,"inq_nbju_ID")
status = NF90_INQ_VARID(fidCOORD,"nbru",nbru_ID)   ; call erreur(status,.TRUE.,"inq_nbru_ID")

status = NF90_GET_VAR(fidCOORD,glamu_ID,glamu) ; call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidCOORD,gphiu_ID,gphiu) ; call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidCOORD,e1u_ID,e1u)     ; call erreur(status,.TRUE.,"getvar_e1u")
status = NF90_GET_VAR(fidCOORD,e2u_ID,e2u)     ; call erreur(status,.TRUE.,"getvar_e2u")
status = NF90_GET_VAR(fidCOORD,nbiu_ID,nbiu)   ; call erreur(status,.TRUE.,"getvar_nbiu")
status = NF90_GET_VAR(fidCOORD,nbju_ID,nbju)   ; call erreur(status,.TRUE.,"getvar_nbju")
status = NF90_GET_VAR(fidCOORD,nbru_ID,nbru)   ; call erreur(status,.TRUE.,"getvar_nbru")

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
          write(file_in_gridU,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_U)
        CASE(192)
          write(file_in_gridU,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_U) 
        CASE(193)
          write(file_in_gridU,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_U)
        CASE(194)
          write(file_in_gridU,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_U)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     write(*,*) 'Looking for existence of ', TRIM(file_in_gridU)
     inquire(file=file_in_gridU, exist=existfile)
     if ( existfile ) then; write(*,*) 'BINGO !'; exit ; endif
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading T,S input dimensions in ', TRIM(file_in_gridU)
    status = NF90_OPEN(TRIM(file_in_gridU),0,fidU)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidU,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidU,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidU,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")
    status = NF90_INQ_DIMID(fidU,"z",dimID_depthu)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"depth",dimID_depthu)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"depthu",dimID_depthu)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"deptht",dimID_depthu)
    call erreur(status,.TRUE.,"inq_dimID_depthu")

    status = NF90_INQUIRE_DIMENSION(fidU,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_x,len=mlon)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_y,len=mlat)         ; call erreur(status,.TRUE.,"inq_dim_y")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_depthu,len=mdepthu) ; call erreur(status,.TRUE.,"inq_dim_depthu")

    ALLOCATE( nav_lon(mlon,mlat), nav_lat(mlon,mlat) )
    ALLOCATE( depthu(mdepthu) )

    status = NF90_INQ_VARID(fidU,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidU,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
    status = NF90_INQ_VARID(fidU,"depth",depthu_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"depthu",depthu_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"deptht",depthu_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"nav_lev",depthu_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"z",depthu_ID)
    call erreur(status,.TRUE.,"inq_depthu_ID")
        
    status = NF90_GET_VAR(fidU,lon_ID,nav_lon)    ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidU,lat_ID,nav_lat)    ; call erreur(status,.TRUE.,"getvar_lat")
    status = NF90_GET_VAR(fidU,depthu_ID,depthu)  ; call erreur(status,.TRUE.,"getvar_depthu")

    status = NF90_CLOSE(fidU)                     ; call erreur(status,.TRUE.,"fin_lecture")

    exit

  ELSEIF ( kday .eq. 31 ) THEN

    write(*,*) 'No velocity file found for first month of ', kyear
    write(*,*) '        >>>>>>>>>>>>>>>>>> stop !!'
    stop

  ENDIF

ENDDO

!- Read umask and e3u in large-scale/global file:
status = NF90_OPEN(TRIM(file_data_mask),0,fidMSKIN);    call erreur(status,.TRUE.,"read mask input") 
ALLOCATE(  umask_GLO(mlon,mlat,mdepthu)  ) 
status = NF90_INQ_VARID(fidMSKIN,"umask",umask_GLO_ID); call erreur(status,.TRUE.,"inq_umask_GLO_ID")
status = NF90_GET_VAR(fidMSKIN,umask_GLO_ID,umask_GLO); call erreur(status,.TRUE.,"getvar_umask_GLO")
status = NF90_CLOSE(fidMSKIN);                          call erreur(status,.TRUE.,"end read mask_GLO")


!- Read e2u and e3u in large-scale/global file:
status = NF90_OPEN(TRIM(file_data_zgr),0,fidZGRIN);     call erreur(status,.TRUE.,"read mask input")
ALLOCATE(  e3u_GLO(mlon,mlat,mdepthu)  )
status = NF90_INQ_VARID(fidZGRIN,"e3u",e3u_GLO_ID) 
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidZGRIN,"e3u_0",e3u_GLO_ID)
call erreur(status,.TRUE.,"inq_e3u_GLO_ID")
status = NF90_GET_VAR(fidZGRIN,e3u_GLO_ID,e3u_GLO);     call erreur(status,.TRUE.,"getvar_e3u_GLO")
status = NF90_CLOSE(fidZGRIN);                          call erreur(status,.TRUE.,"end read mask_GLO")
!-
status = NF90_OPEN(TRIM(file_data_hgr),0,fidHGRIN);     call erreur(status,.TRUE.,"read mask input")
ALLOCATE(  e2u_GLO(mlon,mlat)  )
status = NF90_INQ_VARID(fidHGRIN,"e2u",e2u_GLO_ID);     call erreur(status,.TRUE.,"inq_e2u_GLO_ID")
status = NF90_GET_VAR(fidHGRIN,e2u_GLO_ID,e2u_GLO);     call erreur(status,.TRUE.,"getvar_e2u_GLO")
status = NF90_CLOSE(fidHGRIN);                          call erreur(status,.TRUE.,"end read mask_GLO")


!--

ALLOCATE( nav_lon_bdy(mxbu,1), nav_lat_bdy(mxbu,1) )

do kbdy=1,mxbu
  iGLO=NINT(FLOAT(nbiu(kbdy,1)+imin_ORCA12-1-bi)/ai)
  jGLO=NINT(FLOAT(nbju(kbdy,1)+jmin_ORCA12-1-bj)/aj)
  nav_lon_bdy(kbdy,1) = nav_lon(iGLO,jGLO)
  nav_lat_bdy(kbdy,1) = nav_lat(iGLO,jGLO)
enddo

!--

write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir ',a,'/BDY')
CALL system(TRIM(command_str))

!=================================================================================
! 3- Process all gridU files over specified period
!=================================================================================

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridU,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_U)
        CASE(192)
          write(file_in_gridU,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_U) 
        CASE(193)
          write(file_in_gridU,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_U)
        CASE(194)
          write(file_in_gridU,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_U)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_gridU, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          401 FORMAT(a,'/BDY/bdyU_u2d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridU2d,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
          501 FORMAT(a,'/BDY/bdyU_u3d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridU3d,501) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 ) then
          402 FORMAT(a,'/BDY/bdyU_u2d_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridU2d,402) TRIM(config_dir), kyear, kmonth, TRIM(config)
          502 FORMAT(a,'/BDY/bdyU_u3d_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridU3d,502) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_bdy_gridU3d and file_bdy_gridU3d  >>>> stop'
          stop
        endif

        ALLOCATE( vozocrtx(mlon,mlat,mdepthu,mtime)  )
        ALLOCATE( vobtcrtx(mlon,mlat,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read input velocity : 

        write(*,*) 'Reading velocities in ', TRIM(file_in_gridU)
        
        status = NF90_OPEN(TRIM(file_in_gridU),0,fidU)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidU,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_GET_VAR(fidU,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_ATT(fidU,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidU,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")

        status = NF90_INQ_VARID(fidU,"vozocrtx",vozocrtx_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"uoce",vozocrtx_ID)
        call erreur(status,.TRUE.,"None existing velocity, or unknown variable name")
        status = NF90_GET_VAR(fidU,vozocrtx_ID,vozocrtx)              ; call erreur(status,.TRUE.,"getvar_vozocrtx")

        status = NF90_INQ_VARID(fidU,"e3u",e3u_ID)
        if ( status .ne. 0 ) then
          write(*,*) 'No e3u found in grid_U file => assuming not in vvl mode (i.e. constant e3u)'
          ln_vvl = .false.
        else
          ALLOCATE( e3u(mlon,mlat,mdepthu,mtime)  )
          write(*,*) 'Barotropic component calculated from time-dependent e3u'
          status = NF90_GET_VAR(fidU,e3u_ID,e3u) ; call erreur(status,.TRUE.,"getvar_e3u")
          ln_vvl = .true.
        endif
 
        status = NF90_CLOSE(fidU)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Remove possible NaNs :
 
        do i=1,mlon
        do j=1,mlat
        do k=1,mdepthu
        do l=1,mtime
          if ( .not. abs(vozocrtx(i,j,k,l)) .lt. 99.0 ) then
            vozocrtx(i,j,k,l) = 0.0
          endif
        enddo
        enddo
        enddo
        enddo

        !---------------------------------------
        ! Calculate the barotropic component :
  
        do l=1,mtime
          if ( ln_vvl ) e3u_GLO(:,:,:) = e3u(:,:,:,l)
          do i=1,mlon
          do j=1,mlat
            vobtcrtx(i,j,l) = 0.0
            thic=0.0
            do k=1,mdepthu
              if ( umask_GLO(i,j,k) .eq. 1 ) then
                vobtcrtx(i,j,l) = vobtcrtx(i,j,l) + vozocrtx(i,j,k,l) * e3u_GLO(i,j,k) * umask_GLO(i,j,k)
                thic            = thic            +                     e3u_GLO(i,j,k) * umask_GLO(i,j,k)
              endif
            enddo
            if ( thic .gt. 1.e-3 ) then
              ! barotropic component :
              vobtcrtx(i,j,l) = vobtcrtx(i,j,l) / thic
            else
              vobtcrtx(i,j,l) = 0.0
            endif
            ! baroclinic component :
            do k=1,mdepthu
              vozocrtx(i,j,k,l) = ( vozocrtx(i,j,k,l) - vobtcrtx(i,j,l) ) * umask_GLO(i,j,k)
            enddo
          enddo
          enddo
        enddo

        !---------------------------------------
        ! Fill values on bdyT :
        ! (we ensure transport conservation in case grids are 
        ! slightly stretched with respect to each other)    
        ! NB: no mask used here, i.e. we assume that velocities = 0 at masked points... 
      
        ALLOCATE( vozocrtx_bdy(mxbu,1,mdepthu,mtime)  )
        ALLOCATE( vobtcrtx_bdy(mxbu,1,mtime)  )

        do kbdy=1,mxbu
          iGLO=NINT(FLOAT(nbiu(kbdy,1)+imin_ORCA12-1-bi)/ai)
          jGLO=NINT(FLOAT(nbju(kbdy,1)+jmin_ORCA12-1-bj)/aj)
          write(*,*) kbdy, iGLO, jGLO, nbiu(kbdy,1), nbju(kbdy,1)
          do kt=1,mtime
           do kz=1,mdepthu
             vozocrtx_bdy(kbdy,1,kz,kt) = vozocrtx( iGLO, jGLO, kz, kt ) * e2u_GLO( iGLO, jGLO ) / (e2u(kbdy,1)*aj)
             vobtcrtx_bdy(kbdy,1,   kt) = vobtcrtx( iGLO, jGLO,     kt ) * e2u_GLO( iGLO, jGLO ) / (e2u(kbdy,1)*aj)
           enddo
          enddo
        enddo

        !--------------------------------------
        ! Write BDY netcdf file for barotropic velocity :

        write(*,*) 'Creating ', TRIM(file_bdy_gridU2d)
        status = NF90_CREATE(TRIM(file_bdy_gridU2d),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidM,"xbU",mxbU,dimID_xbU)                             ; call erreur(status,.TRUE.,"def_dimID_xbU")

        status = NF90_DEF_VAR(fidM,"vobtcrtx",NF90_FLOAT,(/dimID_xbU,dimID_yb,dimID_time_counter/),vobtcrtx_ID)
        call erreur(status,.TRUE.,"def_var_vobtcrtx_ID")
        status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_xbU,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_xbU,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"nbrdta",NF90_INT,(/dimID_xbU,dimID_yb/),nbru_ID)
        call erreur(status,.TRUE.,"def_var_nbru_ID")
        status = NF90_DEF_VAR(fidM,"nbjdta",NF90_INT,(/dimID_xbU,dimID_yb/),nbju_ID)
        call erreur(status,.TRUE.,"def_var_nbju_ID")
        status = NF90_DEF_VAR(fidM,"nbidta",NF90_INT,(/dimID_xbU,dimID_yb/),nbiu_ID)
        call erreur(status,.TRUE.,"def_var_nbiu_ID")
       
        status = NF90_PUT_ATT(fidM,vobtcrtx_ID,"long_name","Barotropic velocity")   ; call erreur(status,.TRUE.,"put_att_vobtcrtx_ID")
        status = NF90_PUT_ATT(fidM,vobtcrtx_ID,"units","m/s")                       ; call erreur(status,.TRUE.,"put_att_vobtcrtx_ID")
        status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,nbru_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidM,nbru_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidM,nbju_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidM,nbju_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidM,nbiu_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        status = NF90_PUT_ATT(fidM,nbiu_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bdy_gridU.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidM,vobtcrtx_ID,vobtcrtx_bdy) ; call erreur(status,.TRUE.,"var_vobtcrtx_ID")
        status = NF90_PUT_VAR(fidM,nav_lat_ID,nav_lat_bdy)   ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidM,nav_lon_ID,nav_lon_bdy)   ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidM,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidM,nbru_ID,nbru)             ; call erreur(status,.TRUE.,"var_nbru_ID")
        status = NF90_PUT_VAR(fidM,nbju_ID,nbju)             ; call erreur(status,.TRUE.,"var_nbju_ID")
        status = NF90_PUT_VAR(fidM,nbiu_ID,nbiu)             ; call erreur(status,.TRUE.,"var_nbiu_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close BDY file")

        !--------------------------------------
        ! Write BDY netcdf file for baroclinic velocity :

        write(*,*) 'Creating ', TRIM(file_bdy_gridU3d)
        status = NF90_CREATE(TRIM(file_bdy_gridU3d),NF90_NOCLOBBER,fidN) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidN,"depthu",mdepthu,dimID_depthu)                    ; call erreur(status,.TRUE.,"def_dimID_depthu")
        status = NF90_DEF_DIM(fidN,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidN,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidN,"xbU",mxbU,dimID_xbU)                             ; call erreur(status,.TRUE.,"def_dimID_xbU")

        status = NF90_DEF_VAR(fidN,"vozocrtx",NF90_FLOAT,(/dimID_xbU,dimID_yb,dimID_depthu,dimID_time_counter/),vozocrtx_ID)
        call erreur(status,.TRUE.,"def_var_vozocrtx_ID")
        status = NF90_DEF_VAR(fidN,"depthu",NF90_FLOAT,(/dimID_depthu/),depthu_ID)
        call erreur(status,.TRUE.,"def_var_depthu_ID")
        status = NF90_DEF_VAR(fidN,"nav_lat",NF90_FLOAT,(/dimID_xbU,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidN,"nav_lon",NF90_FLOAT,(/dimID_xbU,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidN,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidN,"nbrdta",NF90_INT,(/dimID_xbU,dimID_yb/),nbru_ID)
        call erreur(status,.TRUE.,"def_var_nbru_ID")
        status = NF90_DEF_VAR(fidN,"nbjdta",NF90_INT,(/dimID_xbU,dimID_yb/),nbju_ID)
        call erreur(status,.TRUE.,"def_var_nbju_ID")
        status = NF90_DEF_VAR(fidN,"nbidta",NF90_INT,(/dimID_xbU,dimID_yb/),nbiu_ID)
        call erreur(status,.TRUE.,"def_var_nbiu_ID")
       
        status = NF90_PUT_ATT(fidN,vozocrtx_ID,"long_name","Baroclinic Velocity")   ; call erreur(status,.TRUE.,"put_att_vozocrtx_ID")
        status = NF90_PUT_ATT(fidN,vozocrtx_ID,"units","m/s")                       ; call erreur(status,.TRUE.,"put_att_vozocrtx_ID")
        status = NF90_PUT_ATT(fidN,depthu_ID,"long_name","Vertical T levels")       ; call erreur(status,.TRUE.,"put_att_depthu_ID")
        status = NF90_PUT_ATT(fidN,depthu_ID,"units","m")                           ; call erreur(status,.TRUE.,"put_att_depthu_ID")
        status = NF90_PUT_ATT(fidN,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidN,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidN,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidN,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidN,nbru_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidN,nbru_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidN,nbju_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidN,nbju_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidN,nbiu_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        status = NF90_PUT_ATT(fidN,nbiu_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        
        status = NF90_PUT_ATT(fidN,NF90_GLOBAL,"history","Created using extract_bdy_gridU.f90")
        status = NF90_PUT_ATT(fidN,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidN) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidN,vozocrtx_ID,vozocrtx_bdy) ; call erreur(status,.TRUE.,"var_vozocrtx_ID")
        status = NF90_PUT_VAR(fidN,depthu_ID,depthu)         ; call erreur(status,.TRUE.,"var_depthu_ID")
        status = NF90_PUT_VAR(fidN,nav_lat_ID,nav_lat_bdy)   ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidN,nav_lon_ID,nav_lon_bdy)   ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidN,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidN,nbru_ID,nbru)             ; call erreur(status,.TRUE.,"var_nbru_ID")
        status = NF90_PUT_VAR(fidN,nbju_ID,nbju)             ; call erreur(status,.TRUE.,"var_nbju_ID")
        status = NF90_PUT_VAR(fidN,nbiu_ID,nbiu)             ; call erreur(status,.TRUE.,"var_nbiu_ID")
        
        status = NF90_CLOSE(fidN) ; call erreur(status,.TRUE.,"close BDY file")
        
        !--       
        
        !--       
        DEALLOCATE( vozocrtx, vobtcrtx, time )
        DEALLOCATE( vozocrtx_bdy, vobtcrtx_bdy )
        IF (ln_vvl) DEALLOCATE(e3u)

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
