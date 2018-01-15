program modif                                         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, IGE-CNRS, Jan. 2017
!
! Script to extract the bathymetry from a global grid (e.g. eORCA12).
!
! If ln_coarse_bdy=.true. the bathymetry is coarser (e.g. from ORCA025) along the boundaries.
!
! 0- Initializations
! 1- Read bathymetry in which to extract (e.g. eORCA12)
! 2- Read coarse bathymetry used for consistent bathymetry along boundaries  [if ln_coarse_bdy]
! 3- Find relationship between the two grids (at least where they overlap)   [if ln_coarse_bdy]
! 4- Extract variables on the REG grid
! 5- Manual corrections for WED12      [ if congig == WED12 ]
! 6- Writing new regional bathymetry file
! 7- Reading coordinates on global domain
! 8- Writing coordinates file for the regional domain
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
INTEGER :: fidORCA12, fidORCA025, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, isf_draft_ID, Bathymetry_isf_ID, Bathymetry_ID, &
&          my_ORCA12, mx_ORCA12,  my_ORCA025, mx_ORCA025,  my_REG, mx_REG, imin_ORCA12, imax_ORCA12, jmin_ORCA12, jmax_ORCA12,            &
&          ai, aj, bi, bj, iREG, jREG, jtmp, npts, kk, ki, kj, ni1, ni2, nj1, nj2, pi, pj, kiref, kjref, mx_tmp, my_tmp, e2f_ID, e2v_ID,  &
&          e2u_ID, e2t_ID, e1f_ID, e1v_ID, e1u_ID, e1t_ID, gphif_ID, gphiv_ID, gphiu_ID, gphit_ID, glamf_ID, glamv_ID, glamu_ID, glamt_ID,&
&          fidCOORDreg, fidCOORDpar, i0, j0 

CHARACTER(LEN=150) :: aaa, file_bathy_out, file_coord_out

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: nav_lat_ORCA12, nav_lon_ORCA12, isf_draft_ORCA12, Bathymetry_isf_ORCA12, Bathymetry_ORCA12,     &
&                                          nav_lat_ORCA025, nav_lon_ORCA025, Bathymetry_ORCA025, Bathymetry_isf_ORCA025, isf_draft_ORCA025,&
&                                          nav_lat_REG, nav_lon_REG, isf_draft_REG, Bathymetry_isf_REG, Bathymetry_REG

REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: e2f, e2v, e2u, e2t, e1f, e1v, e1u, e1t, gphif, gphiv, gphiu, gphit, glamf, glamv, glamu, glamt, &
&                                          e2f_REG, e2v_REG, e2u_REG, e2t_REG, e1f_REG, e1v_REG, e1u_REG, e1t_REG,                         &
&                                          gphif_REG, gphiv_REG, gphiu_REG, gphit_REG, glamf_REG, glamv_REG, glamu_REG, glamt_REG

INTEGER(KIND=4), DIMENSION(8) :: imin, imax, jmin

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
imin_ORCA12 = nn_imin_extract
imax_ORCA12 = nn_imax_extract
jmin_ORCA12 = nn_jmin_extract
jmax_ORCA12 = nn_jmax_extract

eps = 1.e-5

!=================================================================================
! 1- Read bathymetry in which to extract (e.g. eORCA12)
!=================================================================================

write(*,*) 'Extracting regional bathymetry from : ', TRIM(file_in_bathy_extract)

status = NF90_OPEN(TRIM(file_in_bathy_extract),0,fidORCA12) ; call erreur(status,.TRUE.,"read_bathy_to_extract") 

status = NF90_INQ_DIMID(fidORCA12,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_ORCA12")
status = NF90_INQ_DIMID(fidORCA12,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_ORCA12")

status = NF90_INQUIRE_DIMENSION(fidORCA12,dimID_y,len=my_ORCA12); call erreur(status,.TRUE.,"inq_dim_y_ORCA12")
status = NF90_INQUIRE_DIMENSION(fidORCA12,dimID_x,len=mx_ORCA12); call erreur(status,.TRUE.,"inq_dim_x_ORCA12")

ALLOCATE(  nav_lat_ORCA12        (mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  nav_lon_ORCA12        (mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  isf_draft_ORCA12      (mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  Bathymetry_isf_ORCA12 (mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  Bathymetry_ORCA12     (mx_ORCA12,my_ORCA12)  ) 

status = NF90_INQ_VARID(fidORCA12,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA12,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA12,"latitude",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID_ORCA12")
status = NF90_INQ_VARID(fidORCA12,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA12,"lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA12,"longitude",nav_lon_ID)
call erreur(status,.TRUE.,"inq_nav_lon_ID_ORCA12")
status = NF90_INQ_VARID(fidORCA12,"Bathymetry",Bathymetry_ID)
call erreur(status,.TRUE.,"inq_Bathymetry_ID_ORCA12")

status = NF90_GET_VAR(fidORCA12,nav_lat_ID,nav_lat_ORCA12); call erreur(status,.TRUE.,"getvar_nav_lat_ORCA12")
status = NF90_GET_VAR(fidORCA12,nav_lon_ID,nav_lon_ORCA12); call erreur(status,.TRUE.,"getvar_nav_lon_ORCA12")
status = NF90_GET_VAR(fidORCA12,Bathymetry_ID,Bathymetry_ORCA12); call erreur(status,.TRUE.,"getvar_Bathymetry_ORCA12")

if ( ln_isfcav ) then
  status = NF90_INQ_VARID(fidORCA12,"isf_draft",isf_draft_ID);              call erreur(status,.TRUE.,"inq_isf_draft_ID_ORCA12")
  status = NF90_INQ_VARID(fidORCA12,"Bathymetry_isf",Bathymetry_isf_ID);    call erreur(status,.TRUE.,"inq_Bathymetry_isf_ID_ORCA12")
  status = NF90_GET_VAR(fidORCA12,isf_draft_ID,isf_draft_ORCA12);           call erreur(status,.TRUE.,"getvar_isf_draft_ORCA12")
  status = NF90_GET_VAR(fidORCA12,Bathymetry_isf_ID,Bathymetry_isf_ORCA12); call erreur(status,.TRUE.,"getvar_Bathymetry_isf_ORCA12")
endif

status = NF90_CLOSE(fidORCA12); call erreur(status,.TRUE.,"close_grid_to_extract")     

if ( ln_coarse_bdy ) then

    !=================================================================================
    ! 2- Read coarse bathymetry used for consistent bathymetry along boundaries
    !=================================================================================
    
    write(*,*) 'Reading coarse bathymetry for consistent boundaries: ', TRIM(file_in_bathy_bdy)
    
    status = NF90_OPEN(TRIM(file_in_bathy_bdy),0,fidORCA025); call erreur(status,.TRUE.,"read_coarse_bathymetry") 
    
    status = NF90_INQ_DIMID(fidORCA025,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_ORCA025")
    status = NF90_INQ_DIMID(fidORCA025,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_ORCA025")
    
    status = NF90_INQUIRE_DIMENSION(fidORCA025,dimID_y,len=my_ORCA025); call erreur(status,.TRUE.,"inq_dim_y_ORCA025")
    status = NF90_INQUIRE_DIMENSION(fidORCA025,dimID_x,len=mx_ORCA025); call erreur(status,.TRUE.,"inq_dim_x_ORCA025")
    
    ALLOCATE(  Bathymetry_ORCA025     (mx_ORCA025,my_ORCA025)  ) 
    ALLOCATE(  nav_lat_ORCA025        (mx_ORCA025,my_ORCA025)  ) 
    ALLOCATE(  nav_lon_ORCA025        (mx_ORCA025,my_ORCA025)  ) 
    
    status = NF90_INQ_VARID(fidORCA025,"Bathymetry",Bathymetry_ID) ; call erreur(status,.TRUE.,"inq_Bathymetry_ID_ORCA025")
    status = NF90_GET_VAR(fidORCA025,Bathymetry_ID,Bathymetry_ORCA025); call erreur(status,.TRUE.,"getvar_Bathymetry_ORCA025")

    if ( ln_isfcav ) then
      ALLOCATE(  Bathymetry_isf_ORCA025 (mx_ORCA025,my_ORCA025)  )
      ALLOCATE(  isf_draft_ORCA025      (mx_ORCA025,my_ORCA025)  )
      status = NF90_INQ_VARID(fidORCA025,"Bathymetry_isf",Bathymetry_isf_ID)
      if ( status .ne. 0 ) then
        Bathymetry_isf_ORCA025(:,:) = Bathymetry_ORCA025(:,:)
      else
        write(*,*) 'Taking Bathymetry_isf from the coarse dataset used for boundaries'
        status = NF90_GET_VAR(fidORCA025,Bathymetry_isf_ID,Bathymetry_isf_ORCA025)
        call erreur(status,.TRUE.,"getvar_Bathymetry_isf_ORCA025")
      endif
      status = NF90_INQ_VARID(fidORCA025,"isf_draft",isf_draft_ID)
      if ( status .ne. 0 ) then
        isf_draft_ORCA025(:,:) = 0.0
      else
        write(*,*) 'Taking isf_draft from the coarse dataset used for boundaries'
        status = NF90_GET_VAR(fidORCA025,isf_draft_ID,isf_draft_ORCA025)
        call erreur(status,.TRUE.,"getvar_isf_draft_ORCA025")
      endif          
    endif

    status = NF90_INQ_VARID(fidORCA025,"nav_lat",nav_lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA025,"lat",nav_lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA025,"latitude",nav_lat_ID)
    call erreur(status,.TRUE.,"inq_nav_lat_ID_ORCA025")
    status = NF90_INQ_VARID(fidORCA025,"nav_lon",nav_lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA025,"lon",nav_lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA025,"longitude",nav_lon_ID)
    call erreur(status,.TRUE.,"inq_nav_lon_ID_ORCA025")
    status = NF90_GET_VAR(fidORCA025,nav_lat_ID,nav_lat_ORCA025);       call erreur(status,.TRUE.,"getvar_nav_lat_ORCA025")
    status = NF90_GET_VAR(fidORCA025,nav_lon_ID,nav_lon_ORCA025);       call erreur(status,.TRUE.,"getvar_nav_lon_ORCA025")
       
    status = NF90_CLOSE(fidORCA025); call erreur(status,.TRUE.,"close_coarse_bathy_file")
        
    !=================================================================================
    ! 3- Find relationship between the two grids (at least where they overlap)
    !=================================================================================
    
    write(*,*) 'Finding relationship between the two grids (at least where they overlap)'
  
    write(*,*) 'Reference point for grid match :', rn_lonref, rn_latref
 
    distmin=10.0
    do ki=1,mx_ORCA025
    do kj=1,my_ORCA025
      dist = sqrt( ( nav_lat_ORCA025(ki,kj) - rn_latref )**2  + ( nav_lon_ORCA025(ki,kj) - rn_lonref )**2 )
      if ( dist .lt. distmin ) then
        distmin=dist
        kiref=ki
        kjref=kj
      endif
    enddo
    enddo
    write(*,*) kiref, kjref, mx_ORCA025, my_ORCA025 

    do pi=1,mx_ORCA12
    do pj=1,my_ORCA12
      if     (       abs( nav_lon_ORCA12(pi,pj) - nav_lon_ORCA025(kiref,kjref) ) .lt. eps &
      &        .and. abs( nav_lat_ORCA12(pi,pj) - nav_lat_ORCA025(kiref,kjref) ) .lt. eps ) then
        ni1=pi
        nj1=pj
      elseif (       abs( nav_lon_ORCA12(pi,pj) - nav_lon_ORCA025(kiref+1,kjref+1) ) .lt. eps &
      &        .and. abs( nav_lat_ORCA12(pi,pj) - nav_lat_ORCA025(kiref+1,kjref+1) ) .lt. eps ) then
        ni2=pi
        nj2=pj
      endif
    enddo
    enddo

    ! i_ORCA12 = ai * i_ORCA025 + bi
    ! j_ORCA12 = aj * j_ORCA025 + bj
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
    ! Adjust domain bounds to match ORCA025 on the points on which BDY
    ! conditions will be applied (i.e. 2nd pt given the 1-pt masked halo):
    
    imin_ORCA12 = ai * FLOOR(   FLOAT(imin_ORCA12-bi)/ai ) + bi - 1
    imax_ORCA12 = ai * CEILING( FLOAT(imax_ORCA12-bi)/ai ) + bi + 1
    jmin_ORCA12 = aj * FLOOR(   FLOAT(jmin_ORCA12-bj)/aj ) + bj - 1
    jmax_ORCA12 = aj * CEILING( FLOAT(jmax_ORCA12-bj)/aj ) + bj + 1
    
endif
    
write(aaa,112) imin_ORCA12, imax_ORCA12, jmin_ORCA12, jmax_ORCA12
112 FORMAT('The area that is actually extracted is (',i4,':',i4,',',i4,':',i4,')')
    
mx_REG = imax_ORCA12 - imin_ORCA12 + 1
my_REG = jmax_ORCA12 - jmin_ORCA12 + 1
    
write(*,*) ' '
write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
write(*,*) ' Note for later:'
write(*,113) imin_ORCA12, imax_ORCA12, jmin_ORCA12, jmax_ORCA12
113 FORMAT('To crop the global files, use:   ncks -F -d x,',i4,',',i4,' -d y,',i4,',',i4,' file_eORCA12.nc file_REG.nc')
write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
write(*,*) ' Do not forget to adapt NEMO s namelist with:'
write(*,*) '    jpidta  =  ', mx_REG
write(*,*) '    jpjdta  =  ', my_REG
write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
write(*,*) ' '
    
!=================================================================================
! 4- Extract variables on the REG grid :
!=================================================================================
    
ALLOCATE(  nav_lat_REG        (mx_REG,my_REG)  )
ALLOCATE(  nav_lon_REG        (mx_REG,my_REG)  )
if ( ln_isfcav) then
  ALLOCATE(  isf_draft_REG      (mx_REG,my_REG)  )
  ALLOCATE(  Bathymetry_isf_REG (mx_REG,my_REG)  )
endif
ALLOCATE(  Bathymetry_REG     (mx_REG,my_REG)  )
    
nav_lat_REG        (1:mx_REG,1:my_REG) = nav_lat_ORCA12        (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
nav_lon_REG        (1:mx_REG,1:my_REG) = nav_lon_ORCA12        (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
if ( ln_isfcav) then
  isf_draft_REG      (1:mx_REG,1:my_REG) = isf_draft_ORCA12      (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
  Bathymetry_isf_REG (1:mx_REG,1:my_REG) = Bathymetry_isf_ORCA12 (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
endif
Bathymetry_REG     (1:mx_REG,1:my_REG) = Bathymetry_ORCA12     (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
    
if ( ln_coarse_bdy ) then
    
    !---------------------------------------
    ! Adjust bathy along the edges of the REG grid :
    
    write(*,*) 'Halo with bathy from coarse resolution...'
    
    !=== put exactly ORCA025 over a npts-point halo :
    do jREG=1,my_REG
      jtmp = NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of ORCA025's grid
        !-- Western BDY :
        do iREG=1,npts
          if ( ln_isfcav) then
            isf_draft_REG     (iREG,jREG) = isf_draft_ORCA025      ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
            Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
          endif
          Bathymetry_REG    (iREG,jREG) = Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
        enddo
        !--- Eastern BDY :
        do iREG=mx_REG-npts+1,mx_REG
          if ( ln_isfcav) then
            isf_draft_REG     (iREG,jREG) = isf_draft_ORCA025      ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
            Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
          endif
          Bathymetry_REG    (iREG,jREG) = Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    !-- Southern BDY :
    do jREG=1,npts
      jtmp = NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of ORCA025's grid
        do iREG=npts+1,mx_REG-npts
          if ( ln_isfcav) then
            isf_draft_REG     (iREG,jREG) = isf_draft_ORCA025      ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
            Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
          endif
          Bathymetry_REG    (iREG,jREG) = Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    !--- Northern BDY :
    do jREG=my_REG-npts+1,my_REG
      jtmp = NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of ORCA025's grid
        do iREG=npts+1,mx_REG-npts
          if ( ln_isfcav) then
            isf_draft_REG     (iREG,jREG) = isf_draft_ORCA025      ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
            Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
          endif
          Bathymetry_REG    (iREG,jREG) = Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp )
        enddo
      endif
    enddo
    
    write(*,*) 'Smooth transition...'

    !=== smooth transition from ORCA025 to ORCA12 (over npts points again) :
    do jREG=npts+1,my_REG-npts
      jtmp = NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
      !-- Western BDY :
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of ORCA025's grid
        do iREG=npts+1,2*npts
          if ( ln_isfcav) then
            if (       isf_draft_ORCA12(iREG+imin_ORCA12-1,jREG+jmin_ORCA12-1) .gt. 1.0       &
            &    .and. isf_draft_ORCA025(NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_REG   (iREG,jREG) = (  (2*npts+1-iREG)   * isf_draft_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
              &                                    + (iREG-npts) * isf_draft_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
            !else
            !  isf_draft_REG   (iREG,jREG) = 0.0
            endif
            if (       Bathymetry_isf_ORCA12(iREG+imin_ORCA12-1,jREG+jmin_ORCA12-1) .gt. 1.0       &
            &    .and. Bathymetry_isf_ORCA025(NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_REG   (iREG,jREG) = (  (2*npts+1-iREG)   * Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
              &                                    + (iREG-npts) * Bathymetry_isf_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
            !else
            !  Bathymetry_isf_REG   (iREG,jREG) = 0.0
            endif
          endif
          Bathymetry_REG    (iREG,jREG) = (  (2*npts+1-iREG)   * Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
          &                                      + (iREG-npts) * Bathymetry_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
        enddo
        !-- Eastern BDY
        do iREG=mx_REG-2*npts+1,mx_REG-npts
          if ( ln_isfcav) then
            if (       isf_draft_ORCA12(iREG+imin_ORCA12-1,jREG+jmin_ORCA12-1) .gt. 1.0        &
            &    .and. isf_draft_ORCA025(NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp) .gt. 1.0 ) then
              isf_draft_REG   (iREG,jREG) = (   (iREG-mx_REG+2*npts) * isf_draft_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
              &                               + (mx_REG-npts+1-iREG) * isf_draft_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
            !else
            !  isf_draft_REG   (iREG,jREG) = 0.0
            endif
            if (       Bathymetry_isf_ORCA12(iREG+imin_ORCA12-1,jREG+jmin_ORCA12-1) .gt. 1.0        &
            &    .and. Bathymetry_isf_ORCA025(NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp) .gt. 1.0 ) then
              Bathymetry_isf_REG   (iREG,jREG) = (   (iREG-mx_REG+2*npts) * Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
              &                               + (mx_REG-npts+1-iREG) * Bathymetry_isf_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
            !else
            !  Bathymetry_isf_REG   (iREG,jREG) = 0.0
            endif
          endif
          Bathymetry_REG    (iREG,jREG) = (   (iREG-mx_REG+2*npts) * Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
          &                                 + (mx_REG-npts+1-iREG) * Bathymetry_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
        enddo
      endif
    enddo
    !-- Southern BDY :
    do jREG=npts+1,2*npts
      jtmp = NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of ORCA025's grid
        do iREG=2*npts+1,mx_REG-2*npts
          if ( ln_isfcav) then
            if (       isf_draft_ORCA12(iREG+imin_ORCA12-1,jREG+jmin_ORCA12-1) .gt. 1.0       &
            &    .and. isf_draft_ORCA025(NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_REG   (iREG,jREG) = (   (2*npts+1-jREG) * isf_draft_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
              &                               +   (jREG-npts)   * isf_draft_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
            !else
            !  isf_draft_REG   (iREG,jREG) = 0.0
            endif
            if (       Bathymetry_isf_ORCA12(iREG+imin_ORCA12-1,jREG+jmin_ORCA12-1) .gt. 1.0       &
            &    .and. Bathymetry_isf_ORCA025(NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_REG   (iREG,jREG) = (   (2*npts+1-jREG) * Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
              &                                    +   (jREG-npts)   * Bathymetry_isf_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
            !else
            !  Bathymetry_isf_REG   (iREG,jREG) = 0.0
            endif
          endif
          Bathymetry_REG    (iREG,jREG) = (   (2*npts+1-jREG) * Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
          &                                 + (jREG-npts)     * Bathymetry_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
        enddo
      endif
    enddo
    !-- Northern BDY :
    do jREG=my_REG-2*npts+1,my_REG-npts
      jtmp = NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj)
      if ( jtmp .gt. 0 ) then !! we do not try to change the bathy southward of ORCA025's grid
        do iREG=2*npts+1,mx_REG-2*npts
          if ( ln_isfcav) then
            if (       isf_draft_ORCA12(iREG+imin_ORCA12-1,jREG+jmin_ORCA12-1) .gt. 1.0       &
            &    .and. isf_draft_ORCA025(NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai),jtmp) .gt. 1.0 ) then
              isf_draft_REG   (iREG,jREG) = (   (jREG-my_REG+2*npts) * isf_draft_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
              &                               + (my_REG-npts+1-jREG) * isf_draft_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
            !else
            !  isf_draft_REG   (iREG,jREG) = 0.0
            endif
            if (       Bathymetry_isf_ORCA12(iREG+imin_ORCA12-1,jREG+jmin_ORCA12-1) .gt. 1.0       &
            &    .and. Bathymetry_isf_ORCA025(NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai),jtmp) .gt. 1.0 ) then
              Bathymetry_isf_REG   (iREG,jREG) = (   (jREG-my_REG+2*npts) * Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
              &                               + (my_REG-npts+1-jREG) * Bathymetry_isf_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
            !else
            !  Bathymetry_isf_REG   (iREG,jREG) = 0.0
            endif
            Bathymetry_isf_REG(iREG,jREG) = (   (jREG-my_REG+2*npts) * Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
            &                                 + (my_REG-npts+1-jREG) * Bathymetry_isf_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
          endif
          Bathymetry_REG    (iREG,jREG) = (   (jREG-my_REG+2*npts) * Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), jtmp ) &
          &                                 + (my_REG-npts+1-jREG) * Bathymetry_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
        enddo
      endif
    enddo
   
endif

!=================================================================================
! 5- Manual corrections for WED12 :
!=================================================================================

if ( TRIM(config) == 'WED12' ) then
  
    ! Note that you DO NOT have to change the following in you change the domain
    ! size
    ! through modifications of nn_imin_extract, nn_imax_extract, ... in the
    ! namelist 

    ! To keep the boxes at the same position:
    i0 = imin_ORCA12 - 2464
    j0 = jmin_ORCA12 -  151

    !! correction to avoid a closed cavity of 2x2x2 pts (identified after first
    !! mesh_mask creation)
    !isf_draft_REG     (i0+241:i0+242,j0+667:j0+668) = 0.0
    !Bathymetry_isf_REG(i0+241:i0+242,j0+667:j0+668) = 0.0
    !
    !! no isf along eastern boundary (adapt manually to adjust more accurately) :
    !isf_draft_REG     (i0+1095:mx_REG,j0+668:j0+703) = 0.0
    !Bathymetry_isf_REG(i0+1095:mx_REG,j0+668:j0+703) = Bathymetry_REG(i0+1095:mx_REG,j0+668:j0+703)

    ! boxes to fill the Bellingshausen Sea : filled over [imin:imax,jmin:my]
    imin = (/   1+i0 , 192+i0 , 213+i0 , 237+i0 , 254+i0 , 275+i0 , 287+i0 , 299+i0 /)
    imax = (/ 191+i0 , 212+i0 , 236+i0 , 253+i0 , 274+i0 , 286+i0 , 298+i0 , &
    &           NINT(FLOAT(325+i0+imin_ORCA12-1-bi)/ai)*ai+bi-imin_ORCA12+1-1 /) ! WARNING: last number must match with ORCA025 (i.e. next unmasked point neads to be on ORCA025) !!
    jmin = (/ 494+j0 , 807+j0 , 835+j0 , 862+j0 , 876+j0 , 894+j0 ,            &
    &           NINT(FLOAT(899+j0+jmin_ORCA12-1-bj)/aj)*aj+bj-jmin_ORCA12+1+1 ,&
    &           NINT(FLOAT(903+j0+jmin_ORCA12-1-bj)/aj)*aj+bj-jmin_ORCA12+1+1 /) ! WARNING: last two numbers must match with ORCA025 (i.e. next unmasked point neads to be on ORCA025) !!

    do kk=1,size(imin)
      imin(kk) = MIN( MAX( 1, imin(kk) ), mx_REG )
      imax(kk) = MIN( MAX( 1, imax(kk) ), mx_REG )
      jmin(kk) = MIN( MAX( 1, jmin(kk) ), my_REG )
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

    !----- put ORCA025 along modified North-Western corner :
    !- bdy_west(1) :
    do jREG=jmin(8),my_REG
      do iREG=imax(8)+1,imax(8)+npts
        isf_draft_REG     (iREG,jREG) = isf_draft_ORCA025      ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) )
        Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) )
        Bathymetry_REG    (iREG,jREG) = Bathymetry_ORCA025     ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) )
      enddo
    enddo
    !- bdy_north(1) :
    do jREG=jmin(8)-npts, jmin(8)-1
      do iREG=imin(8),imax(8)
        isf_draft_REG     (iREG,jREG) = isf_draft_ORCA025      ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) )
        Bathymetry_isf_REG(iREG,jREG) = Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) )
        Bathymetry_REG    (iREG,jREG) = Bathymetry_ORCA025     ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) )
      enddo
    enddo

    !---- smooth transitions :
    !- bdy_west(1) :
    do jREG=jmin(8),my_REG
      do iREG=imax(8)+npts+1,imax(8)+2*npts
        isf_draft_REG     (iREG,jREG) = (   (imax(8)+2*npts-iREG+1) * isf_draft_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) ) &
        &                                 +     (iREG-imax(8)-npts) * isf_draft_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
        Bathymetry_isf_REG(iREG,jREG) = (   (imax(8)+2*npts-iREG+1) * Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) ) &
        &                                 +     (iREG-imax(8)-npts) * Bathymetry_isf_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
        Bathymetry_REG    (iREG,jREG) = (   (imax(8)+2*npts-iREG+1) * Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) ) &
        &                                 +     (iREG-imax(8)-npts) * Bathymetry_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)
      enddo
    enddo
    !- bdy_north(1) :
    do jREG=jmin(8)-2*npts,jmin(8)-npts-1
      do iREG=imin(8),imax(8)
        isf_draft_REG     (iREG,jREG) = ( (jREG-jmin(8)+2*npts+1)    * isf_draft_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) ) &
        &                                      + (jmin(8)-npts-jREG) * isf_draft_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)                     
        Bathymetry_isf_REG(iREG,jREG) = ( (jREG-jmin(8)+2*npts+1)    * Bathymetry_isf_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) ) &
        &                                      + (jmin(8)-npts-jREG) * Bathymetry_isf_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)                     
        Bathymetry_REG    (iREG,jREG) = ( (jREG-jmin(8)+2*npts+1)    * Bathymetry_ORCA025 ( NINT(FLOAT(iREG+imin_ORCA12-1-bi)/ai), NINT(FLOAT(jREG+jmin_ORCA12-1-bj)/aj) ) &
        &                                      + (jmin(8)-npts-jREG) * Bathymetry_ORCA12  ( iREG+imin_ORCA12-1, jREG+jmin_ORCA12-1 ) ) / (npts+1)                     
      enddo
    enddo

    !- fill bellingshausen:
    do kk=1,8
      isf_draft_REG      (imin(kk):imax(kk),jmin(kk):my_REG) = 0.0
      Bathymetry_isf_REG (imin(kk):imax(kk),jmin(kk):my_REG) = 0.0
      Bathymetry_REG     (imin(kk):imax(kk),jmin(kk):my_REG) = 0.0
    enddo

endif

!=================================================================================
! 6- Writing new REG bathymetry file :
!=================================================================================

write(*,*) 'Creating ', TRIM(file_bathy_out)

status = NF90_CREATE(TRIM(file_bathy_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')                     

status = NF90_DEF_DIM(fidM,"y",my_REG,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx_REG,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")

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
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imin_extraction",imin_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL3")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imax_extraction",imax_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL4")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmin_extraction",jmin_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL5")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmax_extraction",jmax_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL6")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidM,nav_lat_ID,nav_lat_REG); call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,nav_lon_REG); call erreur(status,.TRUE.,"var_nav_lon_ID")
if ( ln_isfcav) then
  status = NF90_PUT_VAR(fidM,isf_draft_ID,isf_draft_REG);           call erreur(status,.TRUE.,"var_isf_draft_ID")
  status = NF90_PUT_VAR(fidM,Bathymetry_isf_ID,Bathymetry_isf_REG); call erreur(status,.TRUE.,"var_Bathymetry_isf_ID")
endif
status = NF90_PUT_VAR(fidM,Bathymetry_ID,Bathymetry_REG); call erreur(status,.TRUE.,"var_Bathymetry_ID")

status = NF90_CLOSE(fidM)                    
call erreur(status,.TRUE.,"final")         

DEALLOCATE( nav_lat_ORCA12, nav_lon_ORCA12, isf_draft_ORCA12, Bathymetry_isf_ORCA12, Bathymetry_ORCA12)
DEALLOCATE( nav_lat_ORCA025, nav_lon_ORCA025, Bathymetry_ORCA025)
DEALLOCATE( nav_lat_REG, nav_lon_REG, isf_draft_REG, Bathymetry_isf_REG, Bathymetry_REG)

!=================================================================================
! 7- Reading coordinates on global domain
!=================================================================================

write(*,*) 'Extracting regional domain from : ', TRIM(file_in_coord_extract)

status = NF90_OPEN(TRIM(file_in_coord_extract),0,fidCOORDpar); call erreur(status,.TRUE.,"read input coordinates") 

status = NF90_INQ_DIMID(fidCOORDpar,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidCOORDpar,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
  
status = NF90_INQUIRE_DIMENSION(fidCOORDpar,dimID_y,len=my_tmp); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidCOORDpar,dimID_x,len=mx_tmp); call erreur(status,.TRUE.,"inq_dim_x")
if ( mx_tmp .ne. mx_ORCA12 .or. my_tmp .ne. my_ORCA12 ) then
  write(*,*) '~!@#$%^* ERROR : mismatch between bathymetry and coordinates dimensions >>>>>> stop !!'
  stop
endif

!ALLOCATE(  mask(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  e2f(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  e2v(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  e2u(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  e2t(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  e1f(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  e1v(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  e1u(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  e1t(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  gphif(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  gphiv(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  gphiu(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  gphit(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  glamf(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  glamv(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  glamu(mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  glamt(mx_ORCA12,my_ORCA12)  ) 

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

!ALLOCATE(  mask_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e2f_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e2v_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e2u_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e2t_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e1f_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e1v_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e1u_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e1t_REG(mx_REG,my_REG)  ) 
ALLOCATE(  gphif_REG(mx_REG,my_REG)  ) 
ALLOCATE(  gphiv_REG(mx_REG,my_REG)  ) 
ALLOCATE(  gphiu_REG(mx_REG,my_REG)  ) 
ALLOCATE(  gphit_REG(mx_REG,my_REG)  ) 
ALLOCATE(  glamf_REG(mx_REG,my_REG)  ) 
ALLOCATE(  glamv_REG(mx_REG,my_REG)  ) 
ALLOCATE(  glamu_REG(mx_REG,my_REG)  ) 
ALLOCATE(  glamt_REG(mx_REG,my_REG)  ) 

write(*,*) 'Creating ', TRIM(file_coord_out)
   
status = NF90_CREATE(TRIM(file_coord_out),NF90_NOCLOBBER,fidCOORDreg); call erreur(status,.TRUE.,'create file_coord_out') 

status = NF90_DEF_DIM(fidCOORDreg,"x",mx_REG,dimID_x); call erreur(status,.TRUE.,"def_dimID_coord_x")
status = NF90_DEF_DIM(fidCOORDreg,"y",my_REG,dimID_y); call erreur(status,.TRUE.,"def_dimID_coord_y")

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

status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"history","Created using extract_bathy_coord.f90"); call erreur(status,.TRUE.,"put_att_GLOB1")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"domain",TRIM(aaa));                                call erreur(status,.TRUE.,"put_att_GLOB2")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"imin_extraction",imin_ORCA12);                     call erreur(status,.TRUE.,"put_att_GLOB3")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"imax_extraction",imax_ORCA12);                     call erreur(status,.TRUE.,"put_att_GLOB4")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"jmin_extraction",jmin_ORCA12);                     call erreur(status,.TRUE.,"put_att_GLOB5")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"jmax_extraction",jmax_ORCA12);                     call erreur(status,.TRUE.,"put_att_GLOB6")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"ai",ai);                                           call erreur(status,.TRUE.,"put_att_GLOB7")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"bi",bi);                                           call erreur(status,.TRUE.,"put_att_GLOB8")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"aj",aj);                                           call erreur(status,.TRUE.,"put_att_GLOB9")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"bj",bj);                                           call erreur(status,.TRUE.,"put_att_GLOB0")

status = NF90_ENDDEF(fidCOORDreg); call erreur(status,.TRUE.,"end_definition_coord") 

e2f_REG=e2f(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
e2v_REG=e2v(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
e2u_REG=e2u(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
e2t_REG=e2t(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
e1f_REG=e1f(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
e1v_REG=e1v(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
e1u_REG=e1u(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
e1t_REG=e1t(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
gphif_REG=gphif(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
gphiv_REG=gphiv(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
gphiu_REG=gphiu(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
gphit_REG=gphit(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
glamf_REG=glamf(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
glamv_REG=glamv(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
glamu_REG=glamu(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
glamt_REG=glamt(imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)

!status = NF90_PUT_VAR(fidCOORDreg,mask_ID,mask_REG); call erreur(status,.TRUE.,"var_mask_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2f_ID,e2f_REG);      call erreur(status,.TRUE.,"var_e2f_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2v_ID,e2v_REG);      call erreur(status,.TRUE.,"var_e2v_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2u_ID,e2u_REG);      call erreur(status,.TRUE.,"var_e2u_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2t_ID,e2t_REG);      call erreur(status,.TRUE.,"var_e2t_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1f_ID,e1f_REG);      call erreur(status,.TRUE.,"var_e1f_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1v_ID,e1v_REG);      call erreur(status,.TRUE.,"var_e1v_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1u_ID,e1u_REG);      call erreur(status,.TRUE.,"var_e1u_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1t_ID,e1t_REG);      call erreur(status,.TRUE.,"var_e1t_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphif_ID,gphif_REG);  call erreur(status,.TRUE.,"var_gphif_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphiv_ID,gphiv_REG);  call erreur(status,.TRUE.,"var_gphiv_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphiu_ID,gphiu_REG);  call erreur(status,.TRUE.,"var_gphiu_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphit_ID,gphit_REG);  call erreur(status,.TRUE.,"var_gphit_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamf_ID,glamf_REG);  call erreur(status,.TRUE.,"var_glamf_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamv_ID,glamv_REG);  call erreur(status,.TRUE.,"var_glamv_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamu_ID,glamu_REG);  call erreur(status,.TRUE.,"var_glamu_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamt_ID,glamt_REG);  call erreur(status,.TRUE.,"var_glamt_ID")

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
