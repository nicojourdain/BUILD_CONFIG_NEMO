program modif                                         
                                                      
USE netcdf                                            
                                                      
IMPLICIT NONE                                         

INTEGER ::fidM,  iii, jjj, ki, kj, kk, rs
                                                      
INTEGER :: fidA, status, dimID_test_cast_midlocation_number, dimID_test_cast_midlocation_length, dimID_test_cast_midpressure_number, dimID_test_cast_midpressure_length, dimID_value_test_cast, dimID_Arctic_test_cast_number, dimID_Arctic_test_cast_length, dimID_test_cast_number, dimID_test_cast_length, dimID_nx_cast, dimID_ny_cast, dimID_nz_cast, dimID_nx, dimID_ny, dimID_nz, mtest_cast_midlocation_number, mtest_cast_midlocation_length, mtest_cast_midpressure_number, mtest_cast_midpressure_length, mvalue_test_cast, mArctic_test_cast_number, mArctic_test_cast_length, mtest_cast_number, mtest_cast_length, mnx_cast, mny_cast, mnz_cast, mnx, mny, mnz, ocean_ref_ID, ndepth_ref_ID, deltaSA_ref_ID, SA_ref_ID, SAAR_ref_ID

CHARACTER(LEN=150) :: file_in, file_out                     
                                                      
REAL*8,ALLOCATABLE,DIMENSION(:,:) :: ocean_ref, ndepth_ref, tmp_ocean_ref, tmp_ndepth_ref
 
REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: deltaSA_ref, SA_ref, SAAR_ref, tmp_deltaSA_ref, tmp_SA_ref, tmp_SAAR_ref
                                                      
    file_in  = 'GSW-Fortran/test/gsw_data_v3_0.nc'
    file_out = 'gsw_data_v3_0_to_be_ncks-A.nc'
                                                           
!---------------------------------------                   
! Read netcdf input file :                                 
                                                           
         status = NF90_OPEN(TRIM(file_in),0,fidA)          
         call erreur(status,.TRUE.,"read") 
                                                           
       !Lecture des ID des dimensions qui nous interessent
         status = NF90_INQ_DIMID(fidA,"test_cast_midlocation_number",dimID_test_cast_midlocation_number)
         call erreur(status,.TRUE.,"inq_dimID_test_cast_midlocation_number")
         status = NF90_INQ_DIMID(fidA,"test_cast_midlocation_length",dimID_test_cast_midlocation_length)
         call erreur(status,.TRUE.,"inq_dimID_test_cast_midlocation_length")
         status = NF90_INQ_DIMID(fidA,"test_cast_midpressure_number",dimID_test_cast_midpressure_number)
         call erreur(status,.TRUE.,"inq_dimID_test_cast_midpressure_number")
         status = NF90_INQ_DIMID(fidA,"test_cast_midpressure_length",dimID_test_cast_midpressure_length)
         call erreur(status,.TRUE.,"inq_dimID_test_cast_midpressure_length")
         status = NF90_INQ_DIMID(fidA,"value_test_cast",dimID_value_test_cast)
         call erreur(status,.TRUE.,"inq_dimID_value_test_cast")
         status = NF90_INQ_DIMID(fidA,"Arctic_test_cast_number",dimID_Arctic_test_cast_number)
         call erreur(status,.TRUE.,"inq_dimID_Arctic_test_cast_number")
         status = NF90_INQ_DIMID(fidA,"Arctic_test_cast_length",dimID_Arctic_test_cast_length)
         call erreur(status,.TRUE.,"inq_dimID_Arctic_test_cast_length")
         status = NF90_INQ_DIMID(fidA,"test_cast_number",dimID_test_cast_number)
         call erreur(status,.TRUE.,"inq_dimID_test_cast_number")
         status = NF90_INQ_DIMID(fidA,"test_cast_length",dimID_test_cast_length)
         call erreur(status,.TRUE.,"inq_dimID_test_cast_length")
         status = NF90_INQ_DIMID(fidA,"nx_cast",dimID_nx_cast)
         call erreur(status,.TRUE.,"inq_dimID_nx_cast")
         status = NF90_INQ_DIMID(fidA,"ny_cast",dimID_ny_cast)
         call erreur(status,.TRUE.,"inq_dimID_ny_cast")
         status = NF90_INQ_DIMID(fidA,"nz_cast",dimID_nz_cast)
         call erreur(status,.TRUE.,"inq_dimID_nz_cast")
         status = NF90_INQ_DIMID(fidA,"nx",dimID_nx)
         call erreur(status,.TRUE.,"inq_dimID_nx")
         status = NF90_INQ_DIMID(fidA,"ny",dimID_ny)
         call erreur(status,.TRUE.,"inq_dimID_ny")
         status = NF90_INQ_DIMID(fidA,"nz",dimID_nz)
         call erreur(status,.TRUE.,"inq_dimID_nz")
                                                               
       !Lecture des valeurs des dimensions qui nous interessent
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_test_cast_midlocation_number,len=mtest_cast_midlocation_number)
         call erreur(status,.TRUE.,"inq_dim_test_cast_midlocation_number")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_test_cast_midlocation_length,len=mtest_cast_midlocation_length)
         call erreur(status,.TRUE.,"inq_dim_test_cast_midlocation_length")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_test_cast_midpressure_number,len=mtest_cast_midpressure_number)
         call erreur(status,.TRUE.,"inq_dim_test_cast_midpressure_number")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_test_cast_midpressure_length,len=mtest_cast_midpressure_length)
         call erreur(status,.TRUE.,"inq_dim_test_cast_midpressure_length")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_value_test_cast,len=mvalue_test_cast)
         call erreur(status,.TRUE.,"inq_dim_value_test_cast")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_Arctic_test_cast_number,len=mArctic_test_cast_number)
         call erreur(status,.TRUE.,"inq_dim_Arctic_test_cast_number")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_Arctic_test_cast_length,len=mArctic_test_cast_length)
         call erreur(status,.TRUE.,"inq_dim_Arctic_test_cast_length")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_test_cast_number,len=mtest_cast_number)
         call erreur(status,.TRUE.,"inq_dim_test_cast_number")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_test_cast_length,len=mtest_cast_length)
         call erreur(status,.TRUE.,"inq_dim_test_cast_length")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_nx_cast,len=mnx_cast)
         call erreur(status,.TRUE.,"inq_dim_nx_cast")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_ny_cast,len=mny_cast)
         call erreur(status,.TRUE.,"inq_dim_ny_cast")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_nz_cast,len=mnz_cast)
         call erreur(status,.TRUE.,"inq_dim_nz_cast")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_nx,len=mnx)
         call erreur(status,.TRUE.,"inq_dim_nx")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_ny,len=mny)
         call erreur(status,.TRUE.,"inq_dim_ny")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_nz,len=mnz)
         call erreur(status,.TRUE.,"inq_dim_nz")
                              
         write(*,*) mny,mnx,mnz
 
       !Allocation of arrays : 
         !ALLOCATE( ocean_ref(mny,mnx) ) 
         !ALLOCATE( ndepth_ref(mny,mnx) )
         !ALLOCATE( SAAR_ref(mny,mnx,mnz) )
         !ALLOCATE( SA_ref(mny,mnx,mnz) )
         !ALLOCATE( deltaSA_ref(mny,mnx,mnz) )
         ALLOCATE( ocean_ref(mny,mnx) ) 
         ALLOCATE( ndepth_ref(mny,mnx) )
         ALLOCATE( SAAR_ref(mnz,mny,mnx) )
         ALLOCATE( SA_ref(mnz,mny,mnx) )
         ALLOCATE( deltaSA_ref(mnz,mny,mnx) )
                                 
       !Lecture des ID des variables qui nous interessent
         status = NF90_INQ_VARID(fidA,"ocean_ref",ocean_ref_ID)
         call erreur(status,.TRUE.,"inq_ocean_ref_ID")
         status = NF90_INQ_VARID(fidA,"ndepth_ref",ndepth_ref_ID)
         call erreur(status,.TRUE.,"inq_ndepth_ref_ID")
         status = NF90_INQ_VARID(fidA,"SAAR_ref",SAAR_ref_ID)
         call erreur(status,.TRUE.,"inq_SAAR_ref_ID")
         status = NF90_INQ_VARID(fidA,"SA_ref",SA_ref_ID)
         call erreur(status,.TRUE.,"inq_SA_ref_ID")
         status = NF90_INQ_VARID(fidA,"deltaSA_ref",deltaSA_ref_ID)
         call erreur(status,.TRUE.,"inq_deltaSA_ref_ID")
                                                              
       !Lecture des valeurs des variables qui nous interessent
         status = NF90_GET_VAR(fidA,ocean_ref_ID,ocean_ref)
         call erreur(status,.TRUE.,"getvar_ocean_ref")
         status = NF90_GET_VAR(fidA,ndepth_ref_ID,ndepth_ref)
         call erreur(status,.TRUE.,"getvar_ndepth_ref")
         status = NF90_GET_VAR(fidA,SAAR_ref_ID,SAAR_ref)
         call erreur(status,.TRUE.,"getvar_SAAR_ref")
         status = NF90_GET_VAR(fidA,SA_ref_ID,SA_ref)
         call erreur(status,.TRUE.,"getvar_SA_ref")
         status = NF90_GET_VAR(fidA,deltaSA_ref_ID,deltaSA_ref)
         call erreur(status,.TRUE.,"getvar_deltaSA_ref")
                                                      
     !Fermeture du fichier lu                         
       status = NF90_CLOSE(fidA)                      
       call erreur(status,.TRUE.,"fin_lecture")     
                                                              
!--------------------------------------------------             
! Remove NaNs by interpolation of closest neighbor

ALLOCATE( tmp_SAAR_ref(mnz,mny,mnx) )
do kk=1,mnz
  do ki=1,mnx
  do kj=1,mny
    if ( .not. SAAR_ref(kk,kj,ki) .lt. 1.e9 ) then
      do rs=1,max(mny,mnx)
            iii=MIN(ki+rs,mnx); jjj=kj               
            if ( SAAR_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=kj               
            if ( SAAR_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MIN(kj+rs,mny)
            if ( SAAR_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MAX(kj-rs,1)     
            if ( SAAR_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MIN(kj+rs,mny)
            if ( SAAR_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MAX(kj-rs,1)     
            if ( SAAR_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MIN(kj+rs,mny)
            if ( SAAR_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MAX(kj-rs,1)     
            if ( SAAR_ref(kk,jjj,iii) .lt. 1.e9 ) exit
      enddo
      tmp_SAAR_ref(kk,kj,ki)=SAAR_ref(kk,jjj,iii)
    else
      tmp_SAAR_ref(kk,kj,ki)=SAAR_ref(kk,kj,ki)
    endif
  enddo
  enddo
  SAAR_ref(:,:,:)=tmp_SAAR_ref(:,:,:)
enddo

ALLOCATE( tmp_SA_ref(mnz,mny,mnx) )
do kk=1,mnz
  do ki=1,mnx
  do kj=1,mny
    if ( .not. SA_ref(kk,kj,ki) .lt. 1.e9 ) then
      do rs=1,max(mny,mnx)
            iii=MIN(ki+rs,mnx); jjj=kj               
            if ( SA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=kj               
            if ( SA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MIN(kj+rs,mny)
            if ( SA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MAX(kj-rs,1)     
            if ( SA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MIN(kj+rs,mny)
            if ( SA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MAX(kj-rs,1)     
            if ( SA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MIN(kj+rs,mny)
            if ( SA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MAX(kj-rs,1)     
            if ( SA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
      enddo
      tmp_SA_ref(kk,kj,ki)=SA_ref(kk,jjj,iii)
    else
      tmp_SA_ref(kk,kj,ki)=SA_ref(kk,kj,ki)
    endif
  enddo
  enddo
  SA_ref(:,:,:)=tmp_SA_ref(:,:,:)
enddo

ALLOCATE( tmp_deltaSA_ref(mnz,mny,mnx) )
do kk=1,mnz
  do ki=1,mnx
  do kj=1,mny
    if ( .not. deltaSA_ref(kk,kj,ki) .lt. 1.e9 ) then
      do rs=1,max(mny,mnx)
            iii=MIN(ki+rs,mnx); jjj=kj               
            if ( deltaSA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=kj               
            if ( deltaSA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MIN(kj+rs,mny)
            if ( deltaSA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MAX(kj-rs,1)     
            if ( deltaSA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MIN(kj+rs,mny)
            if ( deltaSA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MAX(kj-rs,1)     
            if ( deltaSA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MIN(kj+rs,mny)
            if ( deltaSA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MAX(kj-rs,1)     
            if ( deltaSA_ref(kk,jjj,iii) .lt. 1.e9 ) exit
      enddo
      tmp_deltaSA_ref(kk,kj,ki)=deltaSA_ref(kk,jjj,iii)
    else
      tmp_deltaSA_ref(kk,kj,ki)=deltaSA_ref(kk,kj,ki)
    endif
  enddo
  enddo
  deltaSA_ref(:,:,:)=tmp_deltaSA_ref(:,:,:)
enddo

ALLOCATE( tmp_ocean_ref(mny,mnx) )
do kk=1,mnz
  do ki=1,mnx
  do kj=1,mny
    if ( .not. ocean_ref(kj,ki) .lt. 1.e9 ) then
      do rs=1,max(mny,mnx)
            iii=MIN(ki+rs,mnx); jjj=kj               
            if ( ocean_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=kj               
            if ( ocean_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MIN(kj+rs,mny)
            if ( ocean_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MAX(kj-rs,1)     
            if ( ocean_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MIN(kj+rs,mny)
            if ( ocean_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MAX(kj-rs,1)     
            if ( ocean_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MIN(kj+rs,mny)
            if ( ocean_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MAX(kj-rs,1)     
            if ( ocean_ref(jjj,iii) .lt. 1.e9 ) exit
      enddo
      tmp_ocean_ref(kj,ki)=ocean_ref(jjj,iii)
    else
      tmp_ocean_ref(kj,ki)=ocean_ref(kj,ki)
    endif
  enddo
  enddo
  ocean_ref(:,:)=tmp_ocean_ref(:,:)
enddo

ALLOCATE( tmp_ndepth_ref(mny,mnx) )
do kk=1,mnz
  do ki=1,mnx
  do kj=1,mny
    if ( .not. ndepth_ref(kj,ki) .lt. 1.e9 ) then
      do rs=1,max(mny,mnx)
            iii=MIN(ki+rs,mnx); jjj=kj               
            if ( ndepth_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=kj               
            if ( ndepth_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MIN(kj+rs,mny)
            if ( ndepth_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=ki           ; jjj=MAX(kj-rs,1)     
            if ( ndepth_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MIN(kj+rs,mny)
            if ( ndepth_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MIN(ki+rs,mnx); jjj=MAX(kj-rs,1)     
            if ( ndepth_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MIN(kj+rs,mny)
            if ( ndepth_ref(jjj,iii) .lt. 1.e9 ) exit
            iii=MAX(ki-rs,1) ; jjj=MAX(kj-rs,1)     
            if ( ndepth_ref(jjj,iii) .lt. 1.e9 ) exit
      enddo
      tmp_ndepth_ref(kj,ki)=ndepth_ref(jjj,iii)
    else
      tmp_ndepth_ref(kj,ki)=ndepth_ref(kj,ki)
    endif
  enddo
  enddo
  ndepth_ref(:,:)=tmp_ndepth_ref(:,:)
enddo

!---------------------------------------                      
! Writing new netcdf file :                                   
                                                              
        status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM)
        call erreur(status,.TRUE.,'create')                     
                                                                
        !Definition des dimensions du fichiers                  
         status = NF90_DEF_DIM(fidM,"test_cast_midlocation_number",mtest_cast_midlocation_number,dimID_test_cast_midlocation_number)
         call erreur(status,.TRUE.,"def_dimID_test_cast_midlocation_number")
         status = NF90_DEF_DIM(fidM,"test_cast_midlocation_length",mtest_cast_midlocation_length,dimID_test_cast_midlocation_length)
         call erreur(status,.TRUE.,"def_dimID_test_cast_midlocation_length")
         status = NF90_DEF_DIM(fidM,"test_cast_midpressure_number",mtest_cast_midpressure_number,dimID_test_cast_midpressure_number)
         call erreur(status,.TRUE.,"def_dimID_test_cast_midpressure_number")
         status = NF90_DEF_DIM(fidM,"test_cast_midpressure_length",mtest_cast_midpressure_length,dimID_test_cast_midpressure_length)
         call erreur(status,.TRUE.,"def_dimID_test_cast_midpressure_length")
         status = NF90_DEF_DIM(fidM,"value_test_cast",mvalue_test_cast,dimID_value_test_cast)
         call erreur(status,.TRUE.,"def_dimID_value_test_cast")
         status = NF90_DEF_DIM(fidM,"Arctic_test_cast_number",mArctic_test_cast_number,dimID_Arctic_test_cast_number)
         call erreur(status,.TRUE.,"def_dimID_Arctic_test_cast_number")
         status = NF90_DEF_DIM(fidM,"Arctic_test_cast_length",mArctic_test_cast_length,dimID_Arctic_test_cast_length)
         call erreur(status,.TRUE.,"def_dimID_Arctic_test_cast_length")
         status = NF90_DEF_DIM(fidM,"test_cast_number",mtest_cast_number,dimID_test_cast_number)
         call erreur(status,.TRUE.,"def_dimID_test_cast_number")
         status = NF90_DEF_DIM(fidM,"test_cast_length",mtest_cast_length,dimID_test_cast_length)
         call erreur(status,.TRUE.,"def_dimID_test_cast_length")
         status = NF90_DEF_DIM(fidM,"nx_cast",mnx_cast,dimID_nx_cast)
         call erreur(status,.TRUE.,"def_dimID_nx_cast")
         status = NF90_DEF_DIM(fidM,"ny_cast",mny_cast,dimID_ny_cast)
         call erreur(status,.TRUE.,"def_dimID_ny_cast")
         status = NF90_DEF_DIM(fidM,"nz_cast",mnz_cast,dimID_nz_cast)
         call erreur(status,.TRUE.,"def_dimID_nz_cast")
         status = NF90_DEF_DIM(fidM,"nx",mnx,dimID_nx)
         call erreur(status,.TRUE.,"def_dimID_nx")
         status = NF90_DEF_DIM(fidM,"ny",mny,dimID_ny)
         call erreur(status,.TRUE.,"def_dimID_ny")
         status = NF90_DEF_DIM(fidM,"nz",mnz,dimID_nz)
         call erreur(status,.TRUE.,"def_dimID_nz")
                                                              
        !Definition des variables                             
         status = NF90_DEF_VAR(fidM,"ocean_ref",NF90_DOUBLE,(/dimID_ny,dimID_nx/),ocean_ref_ID)
         call erreur(status,.TRUE.,"def_var_ocean_ref_ID")
         status = NF90_DEF_VAR(fidM,"ndepth_ref",NF90_DOUBLE,(/dimID_ny,dimID_nx/),ndepth_ref_ID)
         call erreur(status,.TRUE.,"def_var_ndepth_ref_ID")
         status = NF90_DEF_VAR(fidM,"SAAR_ref",NF90_DOUBLE,(/dimID_nz,dimID_ny,dimID_nx/),SAAR_ref_ID)
         call erreur(status,.TRUE.,"def_var_SAAR_ref_ID")
         status = NF90_DEF_VAR(fidM,"SA_ref",NF90_DOUBLE,(/dimID_nz,dimID_ny,dimID_nx/),SA_ref_ID)
         call erreur(status,.TRUE.,"def_var_SA_ref_ID")
         status = NF90_DEF_VAR(fidM,"deltaSA_ref",NF90_DOUBLE,(/dimID_nz,dimID_ny,dimID_nx/),deltaSA_ref_ID)
         call erreur(status,.TRUE.,"def_var_deltaSA_ref_ID")
                                   
        ! Attributs des variables :
         status = NF90_PUT_ATT(fidM,ocean_ref_ID,"standard_name","ocean_atlas_value_in_atlas_at_that_grid")
         call erreur(status,.TRUE.,"put_att_ocean_ref_ID")
         status = NF90_PUT_ATT(fidM,ndepth_ref_ID,"standard_name","number of levels_in_atlas_at_that_grid")
         call erreur(status,.TRUE.,"put_att_ndepth_ref_ID")
         status = NF90_PUT_ATT(fidM,SAAR_ref_ID,"standard_name","Absolute_Salinity_Anomaly_Ratio")
         call erreur(status,.TRUE.,"put_att_SAAR_ref_ID")
         status = NF90_PUT_ATT(fidM,SA_ref_ID,"standard_name","Reference_Salinity_atlas_value")
         call erreur(status,.TRUE.,"put_att_SA_ref_ID")
         status = NF90_PUT_ATT(fidM,deltaSA_ref_ID,"standard_name","Absolute_Salinity_Anomaly_atlas_value")
         call erreur(status,.TRUE.,"put_att_deltaSA_ref_ID")
                                                      
        !Fin des definitions                          
         status = NF90_ENDDEF(fidM)                   
         call erreur(status,.TRUE.,"fin_definition") 
                                                      
        !Valeurs prises par les variables :           
         status = NF90_PUT_VAR(fidM,ocean_ref_ID,ocean_ref)
         call erreur(status,.TRUE.,"var_ocean_ref_ID")
         status = NF90_PUT_VAR(fidM,ndepth_ref_ID,ndepth_ref)
         call erreur(status,.TRUE.,"var_ndepth_ref_ID")
         status = NF90_PUT_VAR(fidM,SAAR_ref_ID,SAAR_ref)
         call erreur(status,.TRUE.,"var_SAAR_ref_ID")
         status = NF90_PUT_VAR(fidM,SA_ref_ID,SA_ref)
         call erreur(status,.TRUE.,"var_SA_ref_ID")
         status = NF90_PUT_VAR(fidM,deltaSA_ref_ID,deltaSA_ref)
         call erreur(status,.TRUE.,"var_deltaSA_ref_ID")
                                                      
        !Fin de l'ecriture                            
         status = NF90_CLOSE(fidM)                    
         call erreur(status,.TRUE.,"final")         

end program modif



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
