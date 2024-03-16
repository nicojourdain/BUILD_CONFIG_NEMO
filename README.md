# NEMO Pre-processing toolbox 

Pre-processing toolbox to build a regional NEMO configuration laterally forced by a global ocean simulation.

Designed for NEMO versions from 3.6 to 4.2.

### Authors

Main contributor: Nicolas Jourdain (IGE-CNRS).

Thanks for useful feedbacks from: 
* Chris Bull (U. Northumbria)
* Tony Payne (U. Bristol)
* Ute Hausmann (LOCEAN)
* Christoph Kittel (IGE/U. Liège)

This package makes use of the former Gibbs Sea Water (GSW) Toolbox, which is no longer supported in fortran and therefore provided here. See the [Licence for the use of the Gibbs SeaWater (GSW) Oceanographic Toolbox](http://www.teos-10.org/pubs/gsw/html/gsw_licence.html). The other sources provided here are protected by a [GNU General Public License](http://github.com/nicojourdain/BUILD_CONFIG_NEMO/blob/master/LICENSE).

### Requirements

Before you start, make sure your system has the following:
* a fortran compiler.
* netcdf and fortran-netcdf libraries.
* nco tools (ncrcat, etc).

### History

* JAN 2022, N. Jourdain : 
        - new variable suffixes (PAR, EXT, CHLD).
        - tests for eAMUXL12.L121

* NOV 2021, N. Jourdain :
	- minor modifications and successful tests for AMUXL025.L75 with NEMO-4.2-RC.

* JUL 2018, N. Jourdain : 
	- new scripts to interpolate bathy and ice draft from lon/lat or stereo data.
	- can now handle interannual runoff.
	- more flexible way to provide file names in the namelist.
	- tests for the AMUXL12 configuration (Amundsen Sea).

* FEB 2017, N. Jourdain :
	- Initial version.

--------------------

### User guide

The previous toolbox user guide (NEMO-3.6 version) can be found [here](http://github.com/nicojourdain/BUILD_CONFIG_NEMO/blob/master/README_details.txt).

A comprehensive example of use for NEMO4-XIOS2 can be found [here](http://nicojourdain.github.io/students_dir/students_nemo4_occigen) (see in particular section 5 that makes use of the present toolbox).

The main idea of these tools is that we start from a parent grid (global or regional) and we define a child grid at a resolution that is 1, 3 or 5 times higher than the parent domain and that is embedded within the parent domain. Then we build all the input files for this child grid.

##### 1- Compilation of the GSW Toolbox

This former fortran version of the GSW toolbox is used to convert from EOS80 (potential temperature, practical salinity) to TEOS10 (conservative temperature, absolute salinity) if needed. You need to compile it as follows:

```bash
cd BUILD_CONFIG_NEMO
cd GSW-Fortran/test
vi makefile  ## check the FC and NETCDF_INCDIR variables
make
./gsw_check  ## to check (some have to work, maybe not all of them)
cd ../..
```

##### 2- Compilation of the other tools and namelist preparation 

First, define the name of your new regional configuration, and create the directory where you will store all the netcdf files created along the following steps, e.g. :
```bash
export CONFIG='AMUXL025.L75'
mkdir $WORKDIR/input/nemo_${CONFIG}
```

Then, edit your own namelist for the preprocessing, and link namelist\_pre (which is the file read by the processing tools) to it:
```bash
vi namelist_${CONFIG}  ## start from provided examples, e.g. adapt config_dir to $WORKDIR/input/nemo_${CONFIG}
ln -s -v namelist_${CONFIG} namelist_pre
```

Then compile the preprocessing tools as follows:
```bash
vi compile_ALL.sh  ## adapt fortran compiler (and maybe netcdf path), then execute it:
./compile_ALL.sh
```

The next steps make use of submit.sh to execute jobs, which you may need to adapt consistently with your architecture (the default is a simple execution):
```bash
vi submit.sh
```

##### 3- Extract coordinates and bathymetry

In the following, the new regional configuration ($CONFIG) is called the CHILD configuration (e.g. 'AMUXL12.L75'). It will be forced laterally by a PARENT configuration (e.g. ORCA025.L75) that is also used to extract the CHILD initial state and make consistent bathymetries at the boundaries. We currently use a global grid to extract the CHILD coordinates, it is referred to as the EXT grid.

_Example_: AMUXL12.L75 (CHILD) grid extracted as a limited part of ORCA12.L75 (EXT), but AMUXL12.L75 simulation embedded in ORCA025.L75 (PARENT) simulation.

*NB:* this section may be replaced by the *NESTING* tools provided with NEMO.

Fill the ```&griddata``` section of the namelist. Then do as follows to generate a coordinate and a bathy_meter files:
```bash
./submit.sh extract_bathy_coord 01
ls -al $WORKDIR/input/nemo_${CONFIG}/coordinates_${CONFIG}.nc  ## check after completion of extract_bathy_coord
ls -al $WORKDIR/input/nemo_${CONFIG}/bathy_meter_${CONFIG}.nc  ## check after completion of extract_bathy_coord
```

If your are happy with the newly created bathymetry, go directly to step 4. If you prefer to replace the bathymetry (and maybe ice shelf draft) with an interpolation from a dataset independent from the EXT grid, fill the ```&bathy_special``` section of the namelist, and do as follows:
```bash
rm -f $WORKDIR/input/nemo_${CONFIG}/bathy_meter_${CONFIG}.nc
## if the dataset is on a lon/lat grid :
./submit.sh extract_bathy_special_lonlat 03 30
## if the dataset is on a stereographic grid (use 60Gb instead of 30Gb for BedMachine):
./submit.sh extract_bathy_special_stereo 03 30
```

##### 4- Extract mesh\_mask.nc

This step makes use of NEMO directly and an exemple is provided [here](https://nicojourdain.github.io/coding_dir/coding_nemo4_occigen). Basically, you need to generate the domain_cfg.nc file before running NEMO, but the present tools make use of mesh\_mask.nc, so you need to specify ```nn_msh = 1``` in the ```&namdom```section of NEMO's namelist to also obtain a mesh\_mask.nc file. Then, for the next steps, place and rename this file to get:
```bash
ls -al $WORKDIR/input/nemo_${CONFIG}/mesh_mask_${CONFIG}.nc
```

##### 5- Create the initial state

To exctract the CHILD initial state (temperature and salinity) from the PARENT simulation, fill the ```&init``` section of the namelist:
```bash 
./submit.sh extract_istate_TS 01
ls -al $WORKDIR/input/nemo_${CONFIG}/istate_TS_${CONFIG}.nc  # check after completion of extract_istate_TS
```

If you also need an initial state for sea ice (concentration, ice thickness, snow thickness):
```bash 
./submit.sh extract_istate_sea_ice 01
ls -al $WORKDIR/input/nemo_${CONFIG}/istate_sea_ice_${CONFIG}.nc  # check after completion of extract_istate_sea_ice
```

##### 6- Create the lateral boundary conditions (u, v, T, S, SSH, sea-ice)

To define the position of the lateral boundaries, fill the ```&bdy``` section of the namelist, as well as ```&bdy_east```, ```&bdy_west```, ```&bdy_north``` and ```&bdy_south```, if relevant. Then:
```bash
./submit.sh build_coordinates_bdy 01  ## -> creates the coordinate file for lateral boundaries
ls -al $WORKDIR/input/nemo_${CONFIG}/coordinates_bdy_${CONFIG}.nc  # check after completion of build_coordinates_bdy
```

Then, fill the ```&bdy_data``` section of the namelist to select the PARENT data that you want to put at the boundaries and:
```bash
./submit.sh extract_bdy_gridT 04 15   ## -> creates T,S bdy files
./submit.sh extract_bdy_gridU 05 15   ## -> creates U   bdy files
./submit.sh extract_bdy_gridV 05 15   ## -> creates V   bdy files
./submit.sh extract_bdy_icemod 01     ## -> creates ice bdy files
./submit.sh extract_bdy_ssh 01        ## -> creates SSH bdy files
ls -lrt $WORKDIR/input/nemo_${CONFIG}/BDY ## to check progress
```

To concatenate to yearly files:
```bash
vi concatenate_yearly_BDY.sh ## adapt years and CONFIG name
./concatenate_yearly_BDY.sh
```

To generate lateral forcing for barotropic tides, fill the ```&bdy_tide``` section of the namelist, then:
```bash
./submit.sh extract_bdy_tides 01
```

##### 7- Preparation of other input files (SSS for restoring, runoff, chlorophyll, zdfiwm)

If you need to do sea surface salinity (SSS) restoring in the CHILD domain, it can be extracted from the PARENT simulation as follows:
```bash
./submit.sh extract_SSS_restoring 03 15
ls -lrt $WORKDIR/input/nemo_${CONFIG}/SSS ## to check progress
./concatenate_yearly_SSS.sh
```

To extract runoff (incl. iceberg melt) data from the PARENT simulation to the CHILD grid:
```bash
## if only liquid runoff:
./submit.sh extract_runoff 03 8
## if both liquid and solid/iceberg ruonff:
./submit.sh extract_runoff_icb 03 8
##
ls -lrt $WORKDIR/input/nemo_${CONFIG}/RNF ## to check progress
./concatenate_yearly_runoff.sh
```

To extract Chlorophyll-A data:
```bash
# to extract chlorophyll from the parent grid:
./submit.sh extract_chloro 01 8
# or, to extract from a regular lon-lat grid:
./submit.sh extract_chloro_from_lonlat 01 8
ls -lrt $WORKDIR/input/nemo_${CONFIG}/chlorophyll_${CONFIG}.nc
```

If you use the internal wave mixing parameterisation (De Lavergne et al.), you can extract it from the parent grid as follows:
```bash
vi namelist_pre # fill &zdfiwm
./submit.sh extract_zdfiwm 01 8 
ls ../nemo_${CONFIG}/zdfiwm_${CONFIG}.nc
```


##### 8- Weights for the interpolation of Surface Boundary Conditions from atmospheric reanalyses

See *WEIGHTS* tool provided with NEMO and example [here](https://nicojourdain.github.io/students_dir/students_nemo4_occigen).
