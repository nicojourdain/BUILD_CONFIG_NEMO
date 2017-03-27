program test

!use netcdf

use gsw_mod_kinds
use gsw_mod_netcdf
use gsw_mod_toolbox

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

IMPLICIT NONE

REAL*8,DIMENSION(2,2) :: salout, sal, dep, lat, lon

call gsw_saar_init (.true.)

sal(:,:)=34.345367
dep(:,:)=100.0
lat(:,:)=11.0
lon(:,:)=142.0

salout = gsw_sa_from_sp(sal,dep,lon,lat)

write(*,*) salout

end program test
