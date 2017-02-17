program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to build netcdf coordinate file for BDY
!
! history : - Feb. 2017: version with namelist (N. Jourdain)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

use gsw_mod_kinds
use gsw_mod_netcdf
use gsw_mod_toolbox

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /bdy_data/ nn_yeari, nn_yearf, data_dir, data_prefix, nn_bdy_eosmatch, &
&                   data_suffix_T, data_suffix_S, data_suffix_U, data_suffix_V,
&                   data_suffix_ssh, data_suffix_ice

CHARACTER(LEN=50)                    :: config
CHARACTER(LEN=150)                   :: config_dir, data_dir, data_dir, data_prefix, data_suffix_T, data_suffix_S, &
&                                       data_suffix_U, data_suffix_V, data_suffix_ssh, data_suffix_ice
INTEGER                              :: nn_yeari, nn_yearf, nn_bdy_eosmatch

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbt, dimID_xbu, dimID_xbv, &
&                                       myb, mxbt, mxbu, mxbv, glamt_ID, gphit_ID, e1t_ID, e2t_ID,   &
&                                       glamu_ID, gphiu_ID, e1u_ID, e2u_ID, glamv_ID, gphiv_ID,      &
&                                       e1v_ID, e2v_ID, nbit_ID, nbjt_ID, nbrt_ID, nbiu_ID, nbju_ID, &
&                                       nbru_ID, nbiv_ID, nbjv_ID, nbrv_ID, mtime, dimID_x, dimID_y, &
&                                       mlon, mlat, mdeptht, kday, kmonth, kyear, kbdy, nfmt,        &
&                                       kt, kz,     lon_ID, lat_ID, deptht_ID, vosaline_ID, fidM,    &
&                                       votemper_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidA,  &
&                                       dimID_time_counter, dimID_deptht, time_ID, dimID_time,       &
&                                       i, j, k, l
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_gridT, file_bdy_gridT, file_mesh_mask
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbit, nbjt, nbrt, nbiu, nbju, nbru, nbiv, nbjv, nbrv    
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: glamt, gphit, e1t, e2t, glamu, gphiu, e1u, e2u, glamv, gphiv,&
&                                       e1v, e2v, nav_lon, nav_lat, nav_lon_bdy, nav_lat_bdy
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:):: votemper, vosaline, votemper_bdy, vosaline_bdy
REAL*4,ALLOCATABLE,DIMENSION(:)      :: deptht
REAL*8,ALLOCATABLE,DIMENSION(:)      :: time
LOGICAL                              :: existfile

!---------------------------------------

write(file_coord,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/coordinates_bdy_',a,'.nc')

write(file_mesh_mask,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/mesh_mask_',a,'.nc')

!- structure of bdyT file name (see file_bdy_gridT) :
401 FORMAT('bdyT_tra_y',i4.4,'d',i3.3,'_AMU12x.nc')

!---------------------------------------                   
! Read BDY coordinates 

status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read coordinate file") 

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbt",dimID_xbt) ; call erreur(status,.TRUE.,"inq_dimID_xbt")
status = NF90_INQ_DIMID(fidCOORD,"xbu",dimID_xbu) ; call erreur(status,.TRUE.,"inq_dimID_xbu")
status = NF90_INQ_DIMID(fidCOORD,"xbv",dimID_xbv) ; call erreur(status,.TRUE.,"inq_dimID_xbv")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbt,len=mxbt) ; call erreur(status,.TRUE.,"inq_dim_xbt")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbu,len=mxbu) ; call erreur(status,.TRUE.,"inq_dim_xbu")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbv,len=mxbv) ; call erreur(status,.TRUE.,"inq_dim_xbv")

ALLOCATE(  glamt(mxbt,myb)  ) 
ALLOCATE(  gphit(mxbt,myb)  ) 
ALLOCATE(  e1t(mxbt,myb)  ) 
ALLOCATE(  e2t(mxbt,myb)  ) 
ALLOCATE(  nbit(mxbt,myb)  ) 
ALLOCATE(  nbjt(mxbt,myb)  ) 
ALLOCATE(  nbrt(mxbt,myb)  ) 

status = NF90_INQ_VARID(fidCOORD,"glamt",glamt_ID) ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidCOORD,"gphit",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidCOORD,"e1t",e1t_ID)     ; call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidCOORD,"e2t",e2t_ID)     ; call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidCOORD,"nbit",nbit_ID)   ; call erreur(status,.TRUE.,"inq_nbit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjt",nbjt_ID)   ; call erreur(status,.TRUE.,"inq_nbjt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrt",nbrt_ID)   ; call erreur(status,.TRUE.,"inq_nbrt_ID")

status = NF90_GET_VAR(fidCOORD,glamt_ID,glamt) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidCOORD,gphit_ID,gphit) ; call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidCOORD,e1t_ID,e1t)     ; call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidCOORD,e2t_ID,e2t)     ; call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidCOORD,nbit_ID,nbit)   ; call erreur(status,.TRUE.,"getvar_nbit")
status = NF90_GET_VAR(fidCOORD,nbjt_ID,nbjt)   ; call erreur(status,.TRUE.,"getvar_nbjt")
status = NF90_GET_VAR(fidCOORD,nbrt_ID,nbrt)   ; call erreur(status,.TRUE.,"getvar_nbrt")

status = NF90_CLOSE(fidCOORD) ; call erreur(status,.TRUE.,"close coordinate file")     

!-----------------------------------------------
! Search input file dimensions :

191 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')  ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD_<data_suffix_T>.nc
192 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')           ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_<data_suffix_T>.nc
193 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')        ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD.nc
194 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'.nc')                 ! <data_dir>/YYYY/<data_prefix>_YYYY_MM.nc
195 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')           ! <data_dir>/<data_prefix>_YYYY_MM_DD_<data_suffix_T>.nc
196 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')                    ! <data_dir>/<data_prefix>_YYYY_MM_<data_suffix_T>.nc
197 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')                 ! <data_dir>/<data_prefix>_YYYY_MM_DD.nc
198 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'.nc')                          ! <data_dir>/<data_prefix>_YYYY_MM.nc

kyear=nn_yeari
kmonth=1
DO kday=1,31

  write(file_in_gridT,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
  inquire(file=file_in_gridT, exist=existfile)
  if ( existfile ) then
    nfmt=191
  else
    write(file_in_gridT,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T) 
    inquire(file=file_in_gridT, exist=existfile)
    if ( existfile ) then
      nfmt=192
    else
      write(file_in_gridT,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
      inquire(file=file_in_gridT, exist=existfile)
      if ( existfile ) then
        nfmt=193
      else
        write(file_in_gridT,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        inquire(file=file_in_gridT, exist=existfile)
        if ( existfile ) then
          nfmt=194
        else
          write(file_in_gridT,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
          inquire(file=file_in_gridT, exist=existfile)
          if ( existfile ) then
            nfmt=195
          else
            write(file_in_gridT,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T)
            inquire(file=file_in_gridT, exist=existfile)
            if ( existfile ) then
              nfmt=196
            else
              write(file_in_gridT,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
              inquire(file=file_in_gridT, exist=existfile)
              if ( existfile ) then
                nfmt=197
              else
                write(file_in_gridT,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
                inquire(file=file_in_gridT, exist=existfile)
                if ( existfile ) nfmt=198
              endif
            endif
          endif
        endif
      endif
    endif
  endif

  IF ( existfile ) THEN

    status = NF90_OPEN(TRIM(file_in_gridT),0,fidA)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidA,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidA,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidA,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")
    status = NF90_INQ_DIMID(fidA,"z",dimID_deptht)          ; call erreur(status,.TRUE.,"inq_dimID_deptht")

    status = NF90_INQUIRE_DIMENSION(fidA,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mlon)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=mlat)         ; call erreur(status,.TRUE.,"inq_dim_y")
    status = NF90_INQUIRE_DIMENSION(fidA,dimID_deptht,len=mdeptht) ; call erreur(status,.TRUE.,"inq_dim_deptht")

    ALLOCATE( nav_lon(mlon,mlat), nav_lat(mlon,mlat) )
    ALLOCATE( deptht(mdeptht) )

    status = NF90_INQ_VARID(fidA,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidA,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
    status = NF90_INQ_VARID(fidA,"deptht",deptht_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"depth"depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"nav_lev",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidA,"z",depth_ID)
    call erreur(status,.TRUE.,"inq_deptht_ID")
        
    status = NF90_GET_VAR(fidA,lon_ID,nav_lon)                    ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidA,lat_ID,nav_lat)                    ; call erreur(status,.TRUE.,"getvar_lat")
    status = NF90_GET_VAR(fidA,deptht_ID,deptht)                  ; call erreur(status,.TRUE.,"getvar_deptht")

    status = NF90_CLOSE(fidA)                                     ; call erreur(status,.TRUE.,"fin_lecture")

    exit

  ENDIF

ENDDO

!--

ALLOCATE( nav_lon_bdy(mxbt,1), nav_lat_bdy(mxbt,1) )

do kbdy=1,mxbt
  nav_lon_bdy(kbdy,1) = nav_lon( nbit(kbdy,1) , nbjt(kbdy,1) )
  nav_lat_bdy(kbdy,1) = nav_lat( nbit(kbdy,1) , nbjt(kbdy,1) )
enddo

!-----------------------------------------------
! Process all gridT files over specified period: 

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridT,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
          write(file_in_gridS,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_S)
        CASE(192)
          write(file_in_gridT,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T) 
          write(file_in_gridS,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_S) 
        CASE(193)
          write(file_in_gridT,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridS,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_gridT,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
          write(file_in_gridS,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_gridT,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
          write(file_in_gridS,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_S)
        CASE(196) 
          write(file_in_gridT,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T)
          write(file_in_gridS,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_S)
        CASE(197) 
          write(file_in_gridT,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridS,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_gridT,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
          write(file_in_gridS,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_gridT, exist=existfile)

      IF ( existfile ) THEN

        write(file_bdy_gridT,401) kyear, kday

        ALLOCATE( votemper(mlon,mlat,mdeptht,mtime)  )
        ALLOCATE( vosaline(mlon,mlat,mdeptht,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read input temperature :

        write(*,*) 'Reading ', TRIM(file_in_gridT)
        
        status = NF90_OPEN(TRIM(file_in_gridT),0,fidA)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidA,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidA,"votemper",votemper_ID)          ; call erreur(status,.TRUE.,"inq_votemper_ID")
        status = NF90_INQ_VARID(fidA,"vosaline",vosaline_ID)          ; call erreur(status,.TRUE.,"inq_vosaline_ID")
        
        status = NF90_GET_VAR(fidA,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_VAR(fidA,votemper_ID,votemper)              ; call erreur(status,.TRUE.,"getvar_votemper")
        status = NF90_GET_VAR(fidA,vosaline_ID,vosaline)              ; call erreur(status,.TRUE.,"getvar_vosaline")

        status = NF90_GET_ATT(fidA,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidA,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidA)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Remove NaNs :
 
        do i=1,mlon
        do j=1,mlat
        do k=1,mdeptht
        do l=1,mtime
          if ( .not. vosaline(i,j,k,l) .gt. 0.5 ) then
            vosaline(i,j,k,l) = 0.0
            votemper(i,j,k,l) = 0.0
          endif
        enddo
        enddo
        enddo
        enddo

        !---------------------------------------
        ! Fill values on bdyT :
      
        ALLOCATE( votemper_bdy(mxbt,1,mdeptht,mtime)  )
        ALLOCATE( vosaline_bdy(mxbt,1,mdeptht,mtime)  )
 
        do kbdy=1,mxbt
          do kt=1,mtime
          do kz=1,mdeptht
            votemper_bdy(kbdy,1,kz,kt) = votemper( nbit(kbdy,1), nbjt(kbdy,1), kz, kt )
            vosaline_bdy(kbdy,1,kz,kt) = vosaline( nbit(kbdy,1), nbjt(kbdy,1), kz, kt )
          enddo
          enddo
        enddo

        !------------------------------------------------
        ! Convert to conservative temperature if needed :

        if ( nn_bdy_eosmatch .eq. 0 ) then
          write(*,*) 'Converting from EOS80 to TEOS10 ...'
          do kbdy=1,mxbt
          do kt=1,mtime
          do kz=1,mdeptht
            if ( vosaline_bdy(kbdy,1,kz,kt) .gt. 0.5 ) then
              vosaline_bdy(kbdy,1,kz,kt) = gsw_sa_from_sp( DBLE(vosaline_bdy(kbdy,1,kz,kt)), DBLE(deptht(kz)), DBLE(nav_lon_bdy(kbdy,1)), DBLE(nav_lat_bdy(kbdy,1)) )
              votemper_bdy(kbdy,1,kz,kt) = gsw_ct_from_pt( DBLE(vosaline_bdy(kbdy,1,kz,kt)), DBLE(votemper_bdy(kbdy,1,kz,kt)) )
            endif
          enddo
          enddo
          enddo
        elseif ( nn_bdy_eosmatch .ne. 1 ) then
          write(*,*) '~!@#$%^* Error: nn_bdy_eosmatch should be 0 or 1 >>>>> stop !!'
          stop
        endif

        !--------------------------------------
        ! Write BDY netcdf file

        write(*,*) 'Creating ', TRIM(file_bdy_gridT)
        status = NF90_CREATE(TRIM(file_bdy_gridT),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidM,"deptht",mdeptht,dimID_deptht)                    ; call erreur(status,.TRUE.,"def_dimID_deptht")
        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidM,"xbT",mxbT,dimID_xbT)                             ; call erreur(status,.TRUE.,"def_dimID_xbT")

        status = NF90_DEF_VAR(fidM,"votemper",NF90_FLOAT,(/dimID_xbT,dimID_yb,dimID_deptht,dimID_time_counter/),votemper_ID)
        call erreur(status,.TRUE.,"def_var_votemper_ID")
        status = NF90_DEF_VAR(fidM,"vosaline",NF90_FLOAT,(/dimID_xbT,dimID_yb,dimID_deptht,dimID_time_counter/),vosaline_ID)
        call erreur(status,.TRUE.,"def_var_vosaline_ID")
        status = NF90_DEF_VAR(fidM,"deptht",NF90_FLOAT,(/dimID_deptht/),deptht_ID)
        call erreur(status,.TRUE.,"def_var_deptht_ID")
        status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_xbT,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_xbT,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"nbrdta",NF90_INT,(/dimID_xbT,dimID_yb/),nbrt_ID)
        call erreur(status,.TRUE.,"def_var_nbrt_ID")
        status = NF90_DEF_VAR(fidM,"nbjdta",NF90_INT,(/dimID_xbT,dimID_yb/),nbjt_ID)
        call erreur(status,.TRUE.,"def_var_nbjt_ID")
        status = NF90_DEF_VAR(fidM,"nbidta",NF90_INT,(/dimID_xbT,dimID_yb/),nbit_ID)
        call erreur(status,.TRUE.,"def_var_nbit_ID")
       
        if ( nn_bdy_eosmatch .eq. 0 ) then 
          status = NF90_PUT_ATT(fidM,votemper_ID,"long_name","Conservative Temperature") ; call erreur(status,.TRUE.,"put_att_votemper_ID")
          status = NF90_PUT_ATT(fidM,votemper_ID,"units","degC")                         ; call erreur(status,.TRUE.,"put_att_votemper_ID")
          status = NF90_PUT_ATT(fidM,vosaline_ID,"long_name","Absolute Salinity")        ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
          status = NF90_PUT_ATT(fidM,vosaline_ID,"units","g/kg")                         ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        else
          status = NF90_PUT_ATT(fidM,votemper_ID,"long_name","Temperature")           ; call erreur(status,.TRUE.,"put_att_votemper_ID")
          status = NF90_PUT_ATT(fidM,votemper_ID,"units","degC")                         ; call erreur(status,.TRUE.,"put_att_votemper_ID")
          status = NF90_PUT_ATT(fidM,vosaline_ID,"long_name","Salinity")              ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
          status = NF90_PUT_ATT(fidM,vosaline_ID,"units","psu")                       ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        endif
        status = NF90_PUT_ATT(fidM,deptht_ID,"long_name","Vertical T levels")       ; call erreur(status,.TRUE.,"put_att_deptht_ID")
        status = NF90_PUT_ATT(fidM,deptht_ID,"units","m")                           ; call erreur(status,.TRUE.,"put_att_deptht_ID")
        status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,nbrt_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
        status = NF90_PUT_ATT(fidM,nbrt_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
        status = NF90_PUT_ATT(fidM,nbjt_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
        status = NF90_PUT_ATT(fidM,nbjt_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
        status = NF90_PUT_ATT(fidM,nbit_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbit_ID")
        status = NF90_PUT_ATT(fidM,nbit_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbit_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bdy_gridT.f90")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidM,votemper_ID,votemper_bdy) ; call erreur(status,.TRUE.,"var_votemper_ID")
        status = NF90_PUT_VAR(fidM,vosaline_ID,vosaline_bdy) ; call erreur(status,.TRUE.,"var_vosaline_ID")
        status = NF90_PUT_VAR(fidM,deptht_ID,deptht)         ; call erreur(status,.TRUE.,"var_deptht_ID")
        status = NF90_PUT_VAR(fidM,nav_lat_ID,nav_lat_bdy)   ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidM,nav_lon_ID,nav_lon_bdy)   ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidM,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidM,nbrt_ID,nbrt)             ; call erreur(status,.TRUE.,"var_nbrt_ID")
        status = NF90_PUT_VAR(fidM,nbjt_ID,nbjt)             ; call erreur(status,.TRUE.,"var_nbjt_ID")
        status = NF90_PUT_VAR(fidM,nbit_ID,nbit)             ; call erreur(status,.TRUE.,"var_nbit_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close BDY file")
        
        !--       
        DEALLOCATE( votemper, vosaline, time )
        DEALLOCATE( votemper_bdy, vosaline_bdy )

      ENDIF

    ENDDO
ENDDO

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
