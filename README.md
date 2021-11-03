## NEMO Pre-processing toolbox 

Pre-processing toolbox to build a regional NEMO configuration laterally forced by a global ocean simulation.

Designed for NEMO versions from 3.6 to 4.2.

##### Authors

Main contributor: Nicolas Jourdain (IGE-CNRS).

Thanks for useful feedbacks from: Chris Bull (U. Northumbria), Tony Payne (U. Bristol), Ute Hausmann (LOCEAN).

This package makes use of the former Gibbs Sea Water (GSW) Toolbox, which is no longer supported in fortran and therefore provided here. See the [Licence for the use of the Gibbs SeaWater (GSW) Oceanographic Toolbox](https://www.teos-10.org/pubs/gsw/html/gsw_licence.html). The other sources provided here are protected by a GNU General Public License.

##### Requirements

Before you start, make sure your system has the following:
* a fortran compiler.
* netcdf and fortran-netcdf libraries.
* nco tools (ncrcat, etc).

##### History

* NOV 21 , N. Jourdain :
	- modifications and tests (ongoing) with NEMO-4.2-RC.

* JUL 18 , N. Jourdain : 
	- new scripts to interpolate bathy and ice draft from lon/lat or stereo data.
	- can now handle interannual runoff.
	- more flexible way to provide file names in the namelist.
	- tests for the AMUXL12 configuration (Amundsen Sea).

##### User guide

The previous toolbox user guide (written for NEMO-3.6) can be found [here](https://github.com/nicojourdain/BUILD_CONFIG_NEMO/blob/master/README_OLD.txt).

A comprehensive example of use for NEMO4-XIOS2 can be found [here](https://nicojourdain.github.io/students_dir/students_nemo4_occigen) (see in particular section 5 that makes use of the present toolbox).
