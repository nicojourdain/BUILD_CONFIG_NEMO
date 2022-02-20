program modif                                         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To extract the bathymetry from a global/large grid ("EXT", e.g. eORCA12) onto a child/regional ("CHLD") grid.
!
! If ln_coarse_bdy=.true. the bathymetry is coarser (which is refered to as the parent grid, i.e. "PAR", 
! e.g. from ORCA025) along the boundaries.
!
! 0- Initializations
! 1- Read bathymetry in which to extract (e.g. eORCA12)
! 2- Read parent bathymetry used for consistent bathymetry along boundaries  [if ln_coarse_bdy]
! 3- Find relationship between the two grids (at least where they overlap)   [if ln_coarse_bdy]
! 4- Extract variables on the CHILD grid
! 5- Manual corrections for WED12      [ if congig == WED12 ]
! 6- Writing new regional bathymetry file
! 7- Reading coordinates on global domain
! 8- Writing coordinates file for the regional domain
!
! History: - Jan. 2017: initial version (N. Jourdain, CNRS-IGE)
!          -
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /griddata/ inputdir, file_in_coord_extract, file_in_bathy_extract, file_bathy_out, ln_coarse_bdy, file_in_bathy_bdy, ln_isfcav, &
& nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract, rn_latref, rn_lonref
INTEGER                               :: nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: inputdir, file_in_bathy_extract, file_in_coord_extract, file_in_bathy_bdy, config_dir
LOGICAL                               :: ln_coarse_bdy, ln_isfcav
REAL(KIND=4)                          :: rn_latref, rn_lonref

!-- local variables :
INTEGER :: fidEXT, fidPAR, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, isf_draft_ID, Bathymetry_isf_ID, Bathymetry_ID, &
&          my_EXT, mx_EXT,  my_PAR, mx_PAR,  my_CHLD, mx_CHLD, imin_EXT, imax_EXT, jmin_EXT, jmax_EXT,            &
&          ai, aj, bi, bj, iCHLD, jCHLD, jtmp, npts, kk, ki, kj, ni1, ni2, nj1, nj2, pi, pj, kiref, kjref, mx_tmp, my_tmp, e2f_ID, e2v_ID,  &
&          e2u_ID, e2t_ID, e1f_ID, e1v_ID, e1u_ID, e1t_ID, gphif_ID, gphiv_ID, gphiu_ID, gphit_ID, glamf_ID, glamv_ID, glamu_ID, glamt_ID,&
&          fidCOORDreg, fidCOORDpar, i0, j0, Nbox 

CHARACTER(LEN=150) :: aaa, file_bathy_out, file_coord_out

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: nav_lat_EXT, nav_lon_EXT, isf_draft_EXT, Bathymetry_isf_EXT, Bathymetry_EXT,     &
&                                          nav_lat_PAR, nav_lon_PAR, Bathymetry_PAR, Bathymetry_isf_PAR, isf_draft_PAR,&
&                                          nav_lat_CHLD, nav_lon_CHLD, isf_draft_CHLD, Bathymetry_isf_CHLD, Bathymetry_CHLD

REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: e2f, e2v, e2u, e2t, e1f, e1v, e1u, e1t, gphif, gphiv, gphiu, gphit, glamf, glamv, glamu, glamt, &
&                                          e2f_CHLD, e2v_CHLD, e2u_CHLD, e2t_CHLD, e1f_CHLD, e1v_CHLD, e1u_CHLD, e1t_CHLD,                         &
&                                          gphif_CHLD, gphiv_CHLD, gphiu_CHLD, gphit_CHLD, glamf_CHLD, glamv_CHLD, glamu_CHLD, glamt_CHLD

INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: imin, imax, jmin, jmax

REAL(KIND=4) :: eps, dist, distmin
 
!=================================================================================
! 0- Initializations 
!=================================================================================

! Default values (replaced with namelist values if specified):
config_dir        = '.'
file_in_bathy_bdy = 'not_used'
ln_coarse_bdy     = .false.
ln_isfcav         = .false.
rn_latref         = -45.0
rn_lonref         = 170.0

!- read namelist values
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=griddata)
CLOSE(1)

! name of regional bathymetry file (output file) :
write(file_bathy_out,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/bathy_meter_',a,'.nc')

! name of regional coordinates file (output file) :
write(file_coord_out,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

!-
imin_EXT = nn_imin_extract
imax_EXT = nn_imax_extract
jmin_EXT = nn_jmin_extract
jmax_EXT = nn_jmax_extract

eps = 1.e-5

!=================================================================================
! 1- Read bathymetry in which to extract (EXT)
!=================================================================================

write(*,*) 'Extracting regional bathymetry from : ', TRIM(file_in_bathy_extract)

status = NF90_OPEN(TRIM(file_in_bathy_extract),0,fidEXT) ; call erreur(status,.TRUE.,"read_bathy_to_extract") 

status = NF90_INQ_DIMID(fidEXT,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_EXT")
status = NF90_INQ_DIMID(fidEXT,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_EXT")

status = NF90_INQUIRE_DIMENSION(fidEXT,dimID_y,len=my_EXT); call erreur(status,.TRUE.,"inq_dim_y_EXT")
status = NF90_INQUIRE_DIMENSION(fidEXT,dimID_x,len=mx_EXT); call erreur(status,.TRUE.,"inq_dim_x_EXT")

ALLOCATE(  nav_lat_EXT        (mx_EXT,my_EXT)  ) 
ALLOCATE(  nav_lon_EXT        (mx_EXT,my_EXT)  ) 
ALLOCATE(  isf_draft_EXT      (mx_EXT,my_EXT)  ) 
ALLOCATE(  Bathymetry_isf_EXT (mx_EXT,my_EXT)  ) 
ALLOCATE(  Bathymetry_EXT     (mx_EXT,my_EXT)  ) 

status = NF90_INQ_VARID(fidEXT,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidEXT,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidEXT,"latitude",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID_EXT")
status = NF90_INQ_VARID(fidEXT,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidEXT,"lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidEXT,"longitude",nav_lon_ID)
call erreur(status,.TRUE.,"inq_nav_lon_ID_EXT")
status = NF90_INQ_VARID(fidEXT,"Bathymetry_isf",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidEXT,"Bathymetry",Bathymetry_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidEXT,"bathy_metry",Bathymetry_ID)
call erreur(status,.TRUE.,"inq_Bathymetry_ID_EXT")

status = NF90_GET_VAR(fidEXT,nav_lat_ID,nav_lat_EXT); call erreur(status,.TRUE.,"getvar_nav_lat_EXT")
status = NF90_GET_VAR(fidEXT,nav_lon_ID,nav_lon_EXT); call erreur(status,.TRUE.,"getvar_nav_lon_EXT")
status = NF90_GET_VAR(fidEXT,Bathymetry_ID,Bathymetry_EXT); call erreur(status,.TRUE.,"getvar_Bathymetry_EXT")

if ( ln_isfcav ) then
  status = NF90_INQ_VARID(fidEXT,"isf_draft",isf_draft_ID)
  if ( status .eq. 0 ) then
    status = NF90_GET_VAR(fidEXT,isf_draft_ID,isf_draft_EXT)
    call erreur(status,.TRUE.,"getvar_isf_draft_EXT")
  else
    write(*,*) 'WARNING : no isf_draft in global bathymetry file !! isf_draft will be set to zero !!!!!!!!!'
    isf_draft_EXT(:,:) = 0.e0
  endif
  status = NF90_INQ_VARID(fidEXT,"Bathymetry_isf",Bathymetry_isf_ID)
  if ( status .eq. 0 ) then
    status = NF90_GET_VAR(fidEXT,Bathymetry_isf_ID,Bathymetry_isf_EXT)
    call erreur(status,.TRUE.,"getvar_Bathymetry_isf_EXT")
  else
    write(*,*) 'WARNING : no Bathymetry_isf in global bathymetry file !! Bathymetry_isf will be set to standard Bathymetry !!!!!!!!!'
    Bathymetry_isf_EXT(:,:) = Bathymetry_EXT(:,:)
  endif
endif

status = NF90_CLOSE(fidEXT); call erreur(status,.TRUE.,"close_grid_to_extract")     

if ( ln_coarse_bdy ) then

    !=================================================================================
    ! 2- Read parent bathymetry used for consistent bathymetry along boundaries (e.g. eORCA025)
    !=================================================================================
    
    write(*,*) 'Reading coarse bathymetry for consistent boundaries: ', TRIM(file_in_bathy_bdy)
    
    status = NF90_OPEN(TRIM(file_in_bathy_bdy),0,fidPAR); call erreur(status,.TRUE.,"read_coarse_bathymetry") 
    
    status = NF90_INQ_DIMID(fidPAR,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_PAR")
    status = NF90_INQ_DIMID(fidPAR,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_PAR")
    
    status = NF90_INQUIRE_DIMENSION(fidPAR,dimID_y,len=my_PAR); call erreur(status,.TRUE.,"inq_dim_y_PAR")
    status = NF90_INQUIRE_DIMENSION(fidPAR,dimID_x,len=mx_PAR); call erreur(status,.TRUE.,"inq_dim_x_PAR")
    
    ALLOCATE(  Bathymetry_PAR     (mx_PAR,my_PAR)  ) 
    ALLOCATE(  nav_lat_PAR        (mx_PAR,my_PAR)  ) 
    ALLOCATE(  nav_lon_PAR        (mx_PAR,my_PAR)  ) 
    
    status = NF90_INQ_VARID(fidPAR,"Bathymetry",Bathymetry_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"bathy_metry",Bathymetry_ID)
    call erreur(status,.TRUE.,"inq_Bathymetry_ID_PAR")
    status = NF90_GET_VAR(fidPAR,Bathymetry_ID,Bathymetry_PAR); call erreur(status,.TRUE.,"getvar_Bathymetry_PAR")

    if ( ln_isfcav ) then
      ALLOCATE(  Bathymetry_isf_PAR (mx_PAR,my_PAR)  )
      ALLOCATE(  isf_draft_PAR      (mx_PAR,my_PAR)  )
      status = NF90_INQ_VARID(fidPAR,"Bathymetry_isf",Bathymetry_isf_ID)
      if ( status .ne. 0 ) then
        Bathymetry_isf_PAR(:,:) = Bathymetry_PAR(:,:)
      else
        write(*,*) 'Taking Bathymetry_isf from the coarse dataset used for boundaries'
        status = NF90_GET_VAR(fidPAR,Bathymetry_isf_ID,Bathymetry_isf_PAR)
        call erreur(status,.TRUE.,"getvar_Bathymetry_isf_PAR")
      endif
      status = NF90_INQ_VARID(fidPAR,"isf_draft",isf_draft_ID)
      if ( status .ne. 0 ) then
        isf_draft_PAR(:,:) = 0.0
      else
        write(*,*) 'Taking isf_draft from the coarse dataset used for boundaries'
        status = NF90_GET_VAR(fidPAR,isf_draft_ID,isf_draft_PAR)
        call erreur(status,.TRUE.,"getvar_isf_draft_PAR")
      endif          
    endif

    status = NF90_INQ_VARID(fidPAR,"nav_lat",nav_lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"lat",nav_lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"latitude",nav_lat_ID)
    call erreur(status,.TRUE.,"inq_nav_lat_ID_PAR")
    status = NF90_INQ_VARID(fidPAR,"nav_lon",nav_lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"lon",nav_lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidPAR,"longitude",nav_lon_ID)
    call erreur(status,.TRUE.,"inq_nav_lon_ID_PAR")
    status = NF90_GET_VAR(fidPAR,nav_lat_ID,nav_lat_PAR);       call erreur(status,.TRUE.,"getvar_nav_lat_PAR")
    status = NF90_GET_VAR(fidPAR,nav_lon_ID,nav_lon_PAR);       call erreur(status,.TRUE.,"getvar_nav_lon_PAR")
       
    status = NF90_CLOSE(fidPAR); call erreur(status,.TRUE.,"close_coarse_bathy_file")
        
    !=================================================================================
    ! 3- Find relationship between the two grids (at least where they overlap)
    !=================================================================================
    
    write(*,*) 'Finding relationship between the two grids (at least where they overlap)'
  
    write(*,*) 'Reference point for grid match :', rn_lonref, rn_latref
 
    distmin=1000.0
    do ki=1,mx_PAR
    do kj=1,my_PAR
      dist = sqrt( ( nav_lat_PAR(ki,kj) - rn_latref )**2  + ( nav_lon_PAR(ki,kj) - rn_lonref )**2 )
      if ( dist .lt. distmin ) then
        distmin=dist
        kiref=ki
        kjref=kj
      endif
    enddo
    enddo
    write(*,*) kiref, kjref, mx_PAR, my_PAR 

    do pi=1,mx_EXT
    do pj=1,my_EXT
      if     (       abs( nav_lon_EXT(pi,pj) - nav_lon_PAR(kiref,kjref) ) .lt. eps &
      &        .and. abs( nav_lat_EXT(pi,pj) - nav_lat_PAR(kiref,kjref) ) .lt. eps ) then
        ni1=pi
        nj1=pj
      elseif (       abs( nav_lon_EXT(pi,pj) - nav_lon_PAR(kiref+1,kjref+1) ) .lt. eps &
      &        .and. abs( nav_lat_EXT(pi,pj) - nav_lat_PAR(kiref+1,kjref+1) ) .lt. eps ) then
        ni2=pi
        nj2=pj
      endif
    enddo
    enddo

    ! i_EXT = ai * i_PAR + bi
    ! j_EXT = aj * j_PAR + bj
    ai = ni2 - ni1
    bi = ni1 - ai*kiref
    aj = nj2 - nj1
    bj = nj1 - aj*kjref
 
    if ( ai .ne. aj ) then
      write(*,*) '~!@#$ ERROR : grid ratio different along x and y   >>>>>>> STOP !!'
      write(*,*) ai, aj, bi, bj
      stop
    else
      write(*,*) 'Grid resolution ratio is : ', ai
      write(*,*) ai, aj, bi, bj
    endif    

    ! width of the halo where the regional bathy is exactly the coarse bathy :
    ! (there is a second halo of same width for the transition)
    npts=CEILING(ai*1.5)  ! in nb of points on the regional grid
    
    !----------------------------------------------------------------------
    ! Adjust domain bounds to match PAR on the points on which BDY
    ! conditions will be applied (i.e. 2nd pt given the 1-pt masked halo):
    
    imin_EXT = ai * FLOOR(   FLOAT(imin_EXT-bi)/ai ) + bi - 1
    imax_EXT = ai * CEILING( FLOAT(imax_EXT-bi)/ai ) + bi + 1
    jmin_EXT = aj * FLOOR(   FLOAT(jmin_EXT-bj)/aj ) + bj - 1
    jmax_EXT = aj * CEILING( FLOAT(jmax_EXT-bj)/aj ) + bj + 1
    
endif
    
write(aaa,112) imin_EXT, imax_EXT, jmin_EXT, jmax_EXT
112 FORMAT('The area that is actually extracted is (',i4,':',i4,',',i4,':',i4,')')
    
mx_CHLD = imax_EXT - imin_EXT + 1
my_CHLD = jmax_EXT - jmin_EXT + 1
    
write(*,*) ' '
write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
write(*,*) ' Note for later (overwritten namelist parameters):'
write(*,113) imin_EXT, imax_EXT, jmin_EXT, jmax_EXT
113 FORMAT('To crop the global files, use:   ncks -F -d x,',i4,',',i4,' -d y,',i4,',',i4,' file_EXT.nc file_CHLD.nc')
write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
write(*,*) ' Do not forget to adapt NEMO s namelist with:'
write(*,*) '    jpidta  =  ', mx_CHLD
write(*,*) '    jpjdta  =  ', my_CHLD
write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
write(*,*) ' '
    
!=================================================================================
! 4- Extract variables on the CHLD grid :
!=================================================================================
    
ALLOCATE(  nav_lat_CHLD        (mx_CHLD,my_CHLD)  )
ALLOCATE(  nav_lon_CHLD        (mx_CHLD,my_CHLD)  )
if ( ln_isfcav) then
  ALLOCATE(  isf_draft_CHLD      (mx_CHLD,my_CHLD)  )
  ALLOCATE(  Bathymetry_isf_CHLD (mx_CHLD,my_CHLD)  )
endif
ALLOCATE(  Bathymetry_CHLD     (mx_CHLD,my_CHLD)  )
    
nav_lat_CHLD        (1:mx_CHLD,1:my_CHLD) = nav_lat_EXT        (imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
nav_lon_CHLD        (1:mx_CHLD,1:my_CHLD) = nav_lon_EXT        (imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
if ( ln_isfcav) then
  isf_draft_CHLD      (1:mx_CHLD,1:my_CHLD) = isf_draft_EXT      (imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
  Bathymetry_isf_CHLD (1:mx_CHLD,1:my_CHLD) = Bathymetry_isf_EXT (imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
endif
Bathymetry_CHLD     (1:mx_CHLD,1:my_CHLD) = Bathymetry_EXT     (imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
    
if ( ln_coarse_bdy ) then
    
    !---------------------------------------
    ! Adjust bathy along the edges of the CHLD grid :
    
    write(*,*) 'Halo with bathy from coarse resolution...'
    
    !=== put exactly PAR over a npts-point halo :
    do jCHLD=1,my_CHLD
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of PAR's grid
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
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of PAR's grid
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
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of PAR's grid
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

    !=== smooth transition from PAR to EXT (over npts points again) :
    do jCHLD=npts+1,my_CHLD-npts
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      !-- Western BDY :
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of PAR's grid
        do iCHLD=npts+1,2*npts
          if ( ln_isfcav) then
            if (       isf_draft_EXT(iCHLD+imin_EXT-1,jCHLD+jmin_EXT-1) .gt. 1.0       &
            &    .and. isf_draft_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_CHLD   (iCHLD,jCHLD) = (  (2*npts+1-iCHLD)   * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                                    + (iCHLD-npts) * isf_draft_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
            !else
            !  isf_draft_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            if (       Bathymetry_isf_EXT(iCHLD+imin_EXT-1,jCHLD+jmin_EXT-1) .gt. 1.0       &
            &    .and. Bathymetry_isf_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_CHLD   (iCHLD,jCHLD) = (  (2*npts+1-iCHLD)   * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                                    + (iCHLD-npts) * Bathymetry_isf_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
            !else
            !  Bathymetry_isf_CHLD   (iCHLD,jCHLD) = 0.0
            endif
          endif
          Bathymetry_CHLD    (iCHLD,jCHLD) = (  (2*npts+1-iCHLD)   * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
          &                                      + (iCHLD-npts) * Bathymetry_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
        enddo
        !-- Eastern BDY
        do iCHLD=mx_CHLD-2*npts+1,mx_CHLD-npts
          if ( ln_isfcav) then
            if (       isf_draft_EXT(iCHLD+imin_EXT-1,jCHLD+jmin_EXT-1) .gt. 1.0        &
            &    .and. isf_draft_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp) .gt. 1.0 ) then
              isf_draft_CHLD   (iCHLD,jCHLD) = (   (iCHLD-mx_CHLD+2*npts) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                               + (mx_CHLD-npts+1-iCHLD) * isf_draft_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
            !else
            !  isf_draft_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            if (       Bathymetry_isf_EXT(iCHLD+imin_EXT-1,jCHLD+jmin_EXT-1) .gt. 1.0        &
            &    .and. Bathymetry_isf_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp) .gt. 1.0 ) then
              Bathymetry_isf_CHLD   (iCHLD,jCHLD) = (   (iCHLD-mx_CHLD+2*npts) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                               + (mx_CHLD-npts+1-iCHLD) * Bathymetry_isf_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
            !else
            !  Bathymetry_isf_CHLD   (iCHLD,jCHLD) = 0.0
            endif
          endif
          Bathymetry_CHLD    (iCHLD,jCHLD) = (   (iCHLD-mx_CHLD+2*npts) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
          &                                 + (mx_CHLD-npts+1-iCHLD) * Bathymetry_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
        enddo
      endif
    enddo
    !-- Southern BDY :
    do jCHLD=npts+1,2*npts
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of PAR's grid
        do iCHLD=2*npts+1,mx_CHLD-2*npts
          if ( ln_isfcav) then
            if (       isf_draft_EXT(iCHLD+imin_EXT-1,jCHLD+jmin_EXT-1) .gt. 1.0       &
            &    .and. isf_draft_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_CHLD   (iCHLD,jCHLD) = (   (2*npts+1-jCHLD) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                               +   (jCHLD-npts)   * isf_draft_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
            !else
            !  isf_draft_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            if (       Bathymetry_isf_EXT(iCHLD+imin_EXT-1,jCHLD+jmin_EXT-1) .gt. 1.0       &
            &    .and. Bathymetry_isf_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_CHLD   (iCHLD,jCHLD) = (   (2*npts+1-jCHLD) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                                    +   (jCHLD-npts)   * Bathymetry_isf_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
            !else
            !  Bathymetry_isf_CHLD   (iCHLD,jCHLD) = 0.0
            endif
          endif
          Bathymetry_CHLD    (iCHLD,jCHLD) = (   (2*npts+1-jCHLD) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
          &                                 + (jCHLD-npts)     * Bathymetry_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
        enddo
      endif
    enddo
    !-- Northern BDY :
    do jCHLD=my_CHLD-2*npts+1,my_CHLD-npts
      jtmp = NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of PAR's grid
        do iCHLD=2*npts+1,mx_CHLD-2*npts
          if ( ln_isfcav) then
            if (       isf_draft_EXT(iCHLD+imin_EXT-1,jCHLD+jmin_EXT-1) .gt. 1.0       &
            &    .and. isf_draft_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_CHLD   (iCHLD,jCHLD) = (   (jCHLD-my_CHLD+2*npts) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                               + (my_CHLD-npts+1-jCHLD) * isf_draft_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
            !else
            !  isf_draft_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            if (       Bathymetry_isf_EXT(iCHLD+imin_EXT-1,jCHLD+jmin_EXT-1) .gt. 1.0       &
            &    .and. Bathymetry_isf_PAR(NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_CHLD   (iCHLD,jCHLD) = (   (jCHLD-my_CHLD+2*npts) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
              &                               + (my_CHLD-npts+1-jCHLD) * Bathymetry_isf_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
            !else
            !  Bathymetry_isf_CHLD   (iCHLD,jCHLD) = 0.0
            endif
            Bathymetry_isf_CHLD(iCHLD,jCHLD) = (   (jCHLD-my_CHLD+2*npts) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
            &                                 + (my_CHLD-npts+1-jCHLD) * Bathymetry_isf_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
          endif
          Bathymetry_CHLD    (iCHLD,jCHLD) = (   (jCHLD-my_CHLD+2*npts) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), jtmp ) &
          &                                 + (my_CHLD-npts+1-jCHLD) * Bathymetry_EXT  ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
        enddo
      endif
    enddo
   
endif

!=================================================================================
! 5- Manual corrections for some configurations :
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
    &           NINT(FLOAT(325+i0+imin_EXT-1-bi)/ai)*ai+bi-imin_EXT+1-1 /) ! WARNING: last number must match with BDY (i.e. next unmasked point neads to be on BDY) !!
    jmin = (/ 494+j0 , 807+j0 , 835+j0 , 862+j0 , 876+j0 , 894+j0 ,            &
    &           NINT(FLOAT(899+j0+jmin_EXT-1-bj)/aj)*aj+bj-jmin_EXT+1+1 ,&
    &           NINT(FLOAT(903+j0+jmin_EXT-1-bj)/aj)*aj+bj-jmin_EXT+1+1 /) ! WARNING: last two numbers must match with BDY (i.e. next unmasked point neads to be on BDY) !!
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

    !----- put coarse data (PAR) along modified North-Western corner :
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

    !---- smooth transitions between PAR and EXT :
    !- bdy_west(1) :
    do jCHLD=jmin(8),my_CHLD
      do iCHLD=imax(8)+npts+1,imax(8)+2*npts
        isf_draft_CHLD     (iCHLD,jCHLD) = (   (imax(8)+2*npts-iCHLD+1) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                 +     (iCHLD-imax(8)-npts) * isf_draft_EXT ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
        Bathymetry_isf_CHLD(iCHLD,jCHLD) = (   (imax(8)+2*npts-iCHLD+1) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                 +     (iCHLD-imax(8)-npts) * Bathymetry_isf_EXT ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
        Bathymetry_CHLD    (iCHLD,jCHLD) = (   (imax(8)+2*npts-iCHLD+1) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                 +     (iCHLD-imax(8)-npts) * Bathymetry_EXT ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)
      enddo
    enddo
    !- bdy_north(1) :
    do jCHLD=jmin(8)-2*npts,jmin(8)-npts-1
      do iCHLD=imin(8),imax(8)
        isf_draft_CHLD     (iCHLD,jCHLD) = ( (jCHLD-jmin(8)+2*npts+1) * isf_draft_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jCHLD) * isf_draft_EXT ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)                     
        Bathymetry_isf_CHLD(iCHLD,jCHLD) = ( (jCHLD-jmin(8)+2*npts+1) * Bathymetry_isf_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jCHLD) * Bathymetry_isf_EXT ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)                     
        Bathymetry_CHLD    (iCHLD,jCHLD) = ( (jCHLD-jmin(8)+2*npts+1) * Bathymetry_PAR ( NINT(FLOAT(iCHLD+imin_EXT-1-bi)/ai), NINT(FLOAT(jCHLD+jmin_EXT-1-bj)/aj) ) &
        &                                   + (jmin(8)-npts-jCHLD) * Bathymetry_EXT ( iCHLD+imin_EXT-1, jCHLD+jmin_EXT-1 ) ) / (npts+1)                     
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

    ! To keep the boxes at the same position:
    i0 = 1771 - imin_EXT
    j0 =   30 - jmin_EXT

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

    ! no isf over a safety zone (2*npts wide halo) from the eastern and western BDY :
    isf_draft_CHLD     (1:2*npts,:) = 0.0
    Bathymetry_isf_CHLD(1:2*npts,:) = Bathymetry_CHLD(1:2*npts,:)
    isf_draft_CHLD     (mx_CHLD-2*npts+1:mx_CHLD,:) = 0.0
    Bathymetry_isf_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:) = Bathymetry_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:)

elseif ( TRIM(config) == 'eAMUXL12.L121' ) then

    write(*,*) 'Special correction for config ', TRIM(config)

    ! to avoid hole in an ice shelf:
    isf_draft_CHLD(557,199:200) = Bathymetry_isf_CHLD(557,199:200)

    ! no isf over a safety zone (2*npts wide halo) from the eastern and western BDY :
    isf_draft_CHLD     (1:2*npts,:) = 0.0
    Bathymetry_isf_CHLD(1:2*npts,:) = Bathymetry_CHLD(1:2*npts,:)
    isf_draft_CHLD     (mx_CHLD-2*npts+1:mx_CHLD,:) = 0.0
    Bathymetry_isf_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:) = Bathymetry_CHLD(mx_CHLD-2*npts+1:mx_CHLD,:)

endif

!=================================================================================
! 6- Writing new CHLD bathymetry file :
!=================================================================================

write(*,*) 'Creating ', TRIM(file_bathy_out)

status = NF90_CREATE(TRIM(file_bathy_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')                     

status = NF90_DEF_DIM(fidM,"y",my_CHLD,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx_CHLD,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")

status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lat_ID); call erreur(status,.TRUE.,"def_var_nav_lat_ID")
status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lon_ID); call erreur(status,.TRUE.,"def_var_nav_lon_ID")
if ( ln_isfcav) then
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
if ( ln_isfcav) then
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

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bathy_coord.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL1")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"domain",TRIM(aaa))             ; call erreur(status,.TRUE.,"put_att_GLOBAL2")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imin_extraction",imin_EXT)  ; call erreur(status,.TRUE.,"put_att_GLOBAL3")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imax_extraction",imax_EXT)  ; call erreur(status,.TRUE.,"put_att_GLOBAL4")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmin_extraction",jmin_EXT)  ; call erreur(status,.TRUE.,"put_att_GLOBAL5")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmax_extraction",jmax_EXT)  ; call erreur(status,.TRUE.,"put_att_GLOBAL6")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidM,nav_lat_ID,nav_lat_CHLD); call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,nav_lon_CHLD); call erreur(status,.TRUE.,"var_nav_lon_ID")
if ( ln_isfcav) then
  status = NF90_PUT_VAR(fidM,isf_draft_ID,isf_draft_CHLD);           call erreur(status,.TRUE.,"var_isf_draft_ID")
  status = NF90_PUT_VAR(fidM,Bathymetry_isf_ID,Bathymetry_isf_CHLD); call erreur(status,.TRUE.,"var_Bathymetry_isf_ID")
endif
status = NF90_PUT_VAR(fidM,Bathymetry_ID,Bathymetry_CHLD); call erreur(status,.TRUE.,"var_Bathymetry_ID")

status = NF90_CLOSE(fidM)                    
call erreur(status,.TRUE.,"final")         

DEALLOCATE( nav_lat_EXT, nav_lon_EXT, isf_draft_EXT, Bathymetry_isf_EXT, Bathymetry_EXT)
DEALLOCATE( nav_lat_PAR, nav_lon_PAR, Bathymetry_PAR)
DEALLOCATE( nav_lat_CHLD, nav_lon_CHLD, isf_draft_CHLD, Bathymetry_isf_CHLD, Bathymetry_CHLD)

!=================================================================================
! 7- Reading coordinates on global domain
!=================================================================================

write(*,*) 'Extracting regional domain from : ', TRIM(file_in_coord_extract)

status = NF90_OPEN(TRIM(file_in_coord_extract),0,fidCOORDpar); call erreur(status,.TRUE.,"read input coordinates") 

status = NF90_INQ_DIMID(fidCOORDpar,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidCOORDpar,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
  
status = NF90_INQUIRE_DIMENSION(fidCOORDpar,dimID_y,len=my_tmp); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidCOORDpar,dimID_x,len=mx_tmp); call erreur(status,.TRUE.,"inq_dim_x")
if ( mx_tmp .ne. mx_EXT .or. my_tmp .ne. my_EXT ) then
  write(*,*) '~!@#$%^* ERROR : mismatch between bathymetry and coordinates dimensions >>>>>> stop !!'
  stop
endif

!ALLOCATE(  mask(mx_EXT,my_EXT)  ) 
ALLOCATE(  e2f(mx_EXT,my_EXT)  ) 
ALLOCATE(  e2v(mx_EXT,my_EXT)  ) 
ALLOCATE(  e2u(mx_EXT,my_EXT)  ) 
ALLOCATE(  e2t(mx_EXT,my_EXT)  ) 
ALLOCATE(  e1f(mx_EXT,my_EXT)  ) 
ALLOCATE(  e1v(mx_EXT,my_EXT)  ) 
ALLOCATE(  e1u(mx_EXT,my_EXT)  ) 
ALLOCATE(  e1t(mx_EXT,my_EXT)  ) 
ALLOCATE(  gphif(mx_EXT,my_EXT)  ) 
ALLOCATE(  gphiv(mx_EXT,my_EXT)  ) 
ALLOCATE(  gphiu(mx_EXT,my_EXT)  ) 
ALLOCATE(  gphit(mx_EXT,my_EXT)  ) 
ALLOCATE(  glamf(mx_EXT,my_EXT)  ) 
ALLOCATE(  glamv(mx_EXT,my_EXT)  ) 
ALLOCATE(  glamu(mx_EXT,my_EXT)  ) 
ALLOCATE(  glamt(mx_EXT,my_EXT)  ) 

!status = NF90_INQ_VARID(fidCOORDpar,"mask",mask_ID);   call erreur(status,.TRUE.,"inq_mask_ID")
status = NF90_INQ_VARID(fidCOORDpar,"e2f",e2f_ID);     call erreur(status,.TRUE.,"inq_e2f_ID")
status = NF90_INQ_VARID(fidCOORDpar,"e2v",e2v_ID);     call erreur(status,.TRUE.,"inq_e2v_ID")
status = NF90_INQ_VARID(fidCOORDpar,"e2u",e2u_ID);     call erreur(status,.TRUE.,"inq_e2u_ID")
status = NF90_INQ_VARID(fidCOORDpar,"e2t",e2t_ID);     call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidCOORDpar,"e1f",e1f_ID);     call erreur(status,.TRUE.,"inq_e1f_ID")
status = NF90_INQ_VARID(fidCOORDpar,"e1v",e1v_ID);     call erreur(status,.TRUE.,"inq_e1v_ID")
status = NF90_INQ_VARID(fidCOORDpar,"e1u",e1u_ID);     call erreur(status,.TRUE.,"inq_e1u_ID")
status = NF90_INQ_VARID(fidCOORDpar,"e1t",e1t_ID);     call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidCOORDpar,"gphif",gphif_ID); call erreur(status,.TRUE.,"inq_gphif_ID")
status = NF90_INQ_VARID(fidCOORDpar,"gphiv",gphiv_ID); call erreur(status,.TRUE.,"inq_gphiv_ID")
status = NF90_INQ_VARID(fidCOORDpar,"gphiu",gphiu_ID); call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidCOORDpar,"gphit",gphit_ID); call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidCOORDpar,"glamf",glamf_ID); call erreur(status,.TRUE.,"inq_glamf_ID")
status = NF90_INQ_VARID(fidCOORDpar,"glamv",glamv_ID); call erreur(status,.TRUE.,"inq_glamv_ID")
status = NF90_INQ_VARID(fidCOORDpar,"glamu",glamu_ID); call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidCOORDpar,"glamt",glamt_ID); call erreur(status,.TRUE.,"inq_glamt_ID")
   
!status = NF90_GET_VAR(fidCOORDpar,mask_ID,mask);   call erreur(status,.TRUE.,"getvar_mask")
status = NF90_GET_VAR(fidCOORDpar,e2f_ID,e2f);     call erreur(status,.TRUE.,"getvar_e2f")
status = NF90_GET_VAR(fidCOORDpar,e2v_ID,e2v);     call erreur(status,.TRUE.,"getvar_e2v")
status = NF90_GET_VAR(fidCOORDpar,e2u_ID,e2u);     call erreur(status,.TRUE.,"getvar_e2u")
status = NF90_GET_VAR(fidCOORDpar,e2t_ID,e2t);     call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidCOORDpar,e1f_ID,e1f);     call erreur(status,.TRUE.,"getvar_e1f")
status = NF90_GET_VAR(fidCOORDpar,e1v_ID,e1v);     call erreur(status,.TRUE.,"getvar_e1v")
status = NF90_GET_VAR(fidCOORDpar,e1u_ID,e1u);     call erreur(status,.TRUE.,"getvar_e1u")
status = NF90_GET_VAR(fidCOORDpar,e1t_ID,e1t);     call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidCOORDpar,gphif_ID,gphif); call erreur(status,.TRUE.,"getvar_gphif")
status = NF90_GET_VAR(fidCOORDpar,gphiv_ID,gphiv); call erreur(status,.TRUE.,"getvar_gphiv")
status = NF90_GET_VAR(fidCOORDpar,gphiu_ID,gphiu); call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidCOORDpar,gphit_ID,gphit); call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidCOORDpar,glamf_ID,glamf); call erreur(status,.TRUE.,"getvar_glamf")
status = NF90_GET_VAR(fidCOORDpar,glamv_ID,glamv); call erreur(status,.TRUE.,"getvar_glamv")
status = NF90_GET_VAR(fidCOORDpar,glamu_ID,glamu); call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidCOORDpar,glamt_ID,glamt); call erreur(status,.TRUE.,"getvar_glamt")

status = NF90_CLOSE(fidCOORDpar); call erreur(status,.TRUE.,"End read coordinates")     
   
!=================================================================================
! 8- Writing coordinates file for the regional domain :                                   
!=================================================================================

!ALLOCATE(  mask_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  e2f_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  e2v_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  e2u_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  e2t_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  e1f_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  e1v_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  e1u_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  e1t_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  gphif_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  gphiv_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  gphiu_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  gphit_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  glamf_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  glamv_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  glamu_CHLD(mx_CHLD,my_CHLD)  ) 
ALLOCATE(  glamt_CHLD(mx_CHLD,my_CHLD)  ) 

write(*,*) 'Creating ', TRIM(file_coord_out)
   
status = NF90_CREATE(TRIM(file_coord_out),NF90_NOCLOBBER,fidCOORDreg); call erreur(status,.TRUE.,'create file_coord_out') 

status = NF90_DEF_DIM(fidCOORDreg,"x",mx_CHLD,dimID_x); call erreur(status,.TRUE.,"def_dimID_coord_x")
status = NF90_DEF_DIM(fidCOORDreg,"y",my_CHLD,dimID_y); call erreur(status,.TRUE.,"def_dimID_coord_y")

!status = NF90_DEF_VAR(fidCOORDreg,"mask",NF90_DOUBLE,(/dimID_x,dimID_y/),mask_ID);   call erreur(status,.TRUE.,"def_var_mask_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e2f",NF90_DOUBLE,(/dimID_x,dimID_y/),e2f_ID);     call erreur(status,.TRUE.,"def_var_e2f_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e2v",NF90_DOUBLE,(/dimID_x,dimID_y/),e2v_ID);     call erreur(status,.TRUE.,"def_var_e2v_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e2u",NF90_DOUBLE,(/dimID_x,dimID_y/),e2u_ID);     call erreur(status,.TRUE.,"def_var_e2u_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e2t",NF90_DOUBLE,(/dimID_x,dimID_y/),e2t_ID);     call erreur(status,.TRUE.,"def_var_e2t_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e1f",NF90_DOUBLE,(/dimID_x,dimID_y/),e1f_ID);     call erreur(status,.TRUE.,"def_var_e1f_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e1v",NF90_DOUBLE,(/dimID_x,dimID_y/),e1v_ID);     call erreur(status,.TRUE.,"def_var_e1v_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e1u",NF90_DOUBLE,(/dimID_x,dimID_y/),e1u_ID);     call erreur(status,.TRUE.,"def_var_e1u_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e1t",NF90_DOUBLE,(/dimID_x,dimID_y/),e1t_ID);     call erreur(status,.TRUE.,"def_var_e1t_ID")
status = NF90_DEF_VAR(fidCOORDreg,"gphif",NF90_DOUBLE,(/dimID_x,dimID_y/),gphif_ID); call erreur(status,.TRUE.,"def_var_gphif_ID")
status = NF90_DEF_VAR(fidCOORDreg,"gphiv",NF90_DOUBLE,(/dimID_x,dimID_y/),gphiv_ID); call erreur(status,.TRUE.,"def_var_gphiv_ID")
status = NF90_DEF_VAR(fidCOORDreg,"gphiu",NF90_DOUBLE,(/dimID_x,dimID_y/),gphiu_ID); call erreur(status,.TRUE.,"def_var_gphiu_ID")
status = NF90_DEF_VAR(fidCOORDreg,"gphit",NF90_DOUBLE,(/dimID_x,dimID_y/),gphit_ID); call erreur(status,.TRUE.,"def_var_gphit_ID")
status = NF90_DEF_VAR(fidCOORDreg,"glamf",NF90_DOUBLE,(/dimID_x,dimID_y/),glamf_ID); call erreur(status,.TRUE.,"def_var_glamf_ID")
status = NF90_DEF_VAR(fidCOORDreg,"glamv",NF90_DOUBLE,(/dimID_x,dimID_y/),glamv_ID); call erreur(status,.TRUE.,"def_var_glamv_ID")
status = NF90_DEF_VAR(fidCOORDreg,"glamu",NF90_DOUBLE,(/dimID_x,dimID_y/),glamu_ID); call erreur(status,.TRUE.,"def_var_glamu_ID")
status = NF90_DEF_VAR(fidCOORDreg,"glamt",NF90_DOUBLE,(/dimID_x,dimID_y/),glamt_ID); call erreur(status,.TRUE.,"def_var_glamt_ID")

status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"history","Created using extract_bathy_coord.f90"); call erreur(status,.TRUE.,"put_att_EXTB1")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"domain",TRIM(aaa));                                call erreur(status,.TRUE.,"put_att_EXTB2")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"imin_extraction",imin_EXT);                        call erreur(status,.TRUE.,"put_att_EXTB3")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"imax_extraction",imax_EXT);                        call erreur(status,.TRUE.,"put_att_EXTB4")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"jmin_extraction",jmin_EXT);                        call erreur(status,.TRUE.,"put_att_EXTB5")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"jmax_extraction",jmax_EXT);                        call erreur(status,.TRUE.,"put_att_EXTB6")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"ai",ai);                                           call erreur(status,.TRUE.,"put_att_EXTB7")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"bi",bi);                                           call erreur(status,.TRUE.,"put_att_EXTB8")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"aj",aj);                                           call erreur(status,.TRUE.,"put_att_EXTB9")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"bj",bj);                                           call erreur(status,.TRUE.,"put_att_EXTB0")

status = NF90_ENDDEF(fidCOORDreg); call erreur(status,.TRUE.,"end_definition_coord") 

e2f_CHLD=e2f(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
e2v_CHLD=e2v(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
e2u_CHLD=e2u(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
e2t_CHLD=e2t(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
e1f_CHLD=e1f(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
e1v_CHLD=e1v(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
e1u_CHLD=e1u(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
e1t_CHLD=e1t(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
gphif_CHLD=gphif(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
gphiv_CHLD=gphiv(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
gphiu_CHLD=gphiu(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
gphit_CHLD=gphit(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
glamf_CHLD=glamf(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
glamv_CHLD=glamv(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
glamu_CHLD=glamu(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)
glamt_CHLD=glamt(imin_EXT:imax_EXT,jmin_EXT:jmax_EXT)

!status = NF90_PUT_VAR(fidCOORDreg,mask_ID,mask_CHLD); call erreur(status,.TRUE.,"var_mask_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2f_ID,e2f_CHLD);      call erreur(status,.TRUE.,"var_e2f_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2v_ID,e2v_CHLD);      call erreur(status,.TRUE.,"var_e2v_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2u_ID,e2u_CHLD);      call erreur(status,.TRUE.,"var_e2u_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2t_ID,e2t_CHLD);      call erreur(status,.TRUE.,"var_e2t_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1f_ID,e1f_CHLD);      call erreur(status,.TRUE.,"var_e1f_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1v_ID,e1v_CHLD);      call erreur(status,.TRUE.,"var_e1v_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1u_ID,e1u_CHLD);      call erreur(status,.TRUE.,"var_e1u_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1t_ID,e1t_CHLD);      call erreur(status,.TRUE.,"var_e1t_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphif_ID,gphif_CHLD);  call erreur(status,.TRUE.,"var_gphif_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphiv_ID,gphiv_CHLD);  call erreur(status,.TRUE.,"var_gphiv_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphiu_ID,gphiu_CHLD);  call erreur(status,.TRUE.,"var_gphiu_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphit_ID,gphit_CHLD);  call erreur(status,.TRUE.,"var_gphit_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamf_ID,glamf_CHLD);  call erreur(status,.TRUE.,"var_glamf_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamv_ID,glamv_CHLD);  call erreur(status,.TRUE.,"var_glamv_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamu_ID,glamu_CHLD);  call erreur(status,.TRUE.,"var_glamu_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamt_ID,glamt_CHLD);  call erreur(status,.TRUE.,"var_glamt_ID")

status = NF90_CLOSE(fidCOORDreg); call erreur(status,.TRUE.,"final")         

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
