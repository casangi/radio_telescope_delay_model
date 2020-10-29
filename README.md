# Radio Telescope Delay Model
Wraps the CALC11 Libary Fortran library found [here](https://bitbucket.alma.cl/projects/ASW/repos/control/browse/thirdPartyPackages). 
Currently only the ```almacalc``` function has been wrapped.

## Installation
```sh
pip install radio_telescope_delay_model
```
## Build
```sh
git clone https://github.com/casangi/radio_telescope_delay_model
cd radio_telescope_delay_model
pip install -e ./
```
## Requirements 
gcc => 4.8.5

## Usage 
```sh
from calc11 import almacalc 
```
The [example notebook](https://colab.research.google.com/github/casangi/radio_telescope_delay_model/blob/main/example_notebooks/Radio_Telescope_Delay_Model_Example.ipynb) shows how to calculate uvw coordinates using ```almacalc``` and compares the results to CASA, Astropy, and data in a measurement set.
