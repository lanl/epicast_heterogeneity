# EpicastHeterogeneity

## Installing

***NOTE***: You **MUST** unzip the file `geo-data.zip` contained in the `EpicastGeoplot/assests` sub-directory before using this code for the first time.

## Usage

### First time
```julia
import Pkg
Pkg.activate("EpicastHeterogeneity")
Pkg.instantiate() # will download / install all deps
import EpicastHeterogeneity as EH

data_dir = "..." # full path to the directory containing the data for this paper
EH.figure3(data_dir);
```

### Subsequently

```julia
import Pkg
Pkg.activate("EpicastHeterogeneity")

import EpicastHeterogeneity as EH

data_dir = "..." # full path to the directory containing the data for this paper
EH.figure3(data_dir);
```