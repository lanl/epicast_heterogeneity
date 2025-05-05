# epicast-heterogeneity

Script ('epicast_heterogeneity.py') for interpreting EpiCast `events.bin' files and plotting racial/ethnic trends of interest. This script generates results reported in 'The value of incorporating population heterogeneity into models of respiratory pathogen spread' by Harris et al. Additional scripts included for preparing input data.

## Usage
Six main functions avaialable after script is imported into workspace:

```python
import epicast_heterogeneity as eh
data_dir = "..." # YOU MUST SET THIS to where ever you put the data

# main figure 4: 35 is the FIPS code for NM
eh.figure_4_gen(data_dir)

# main figure 5
eh.figure_5_gen(data_dir)

# main figure 6
eh.figure_6_gen(data_dir)

# supp figure S1
eh.figure_s1_gen(data_dir)
# alias for: eh.figure_4_gen(data_dir, 35, "ethnicity")

# supp figure S2
eh.figure_s2_gen(data_dir)

# supp figure S3
eh.figure_s3_gen(data_dir)
# alias for: eh.figure_4_gen(data_dir, 12, "ethnicity")
```

These can also be called from the commandline:

```bash

# YOU MUST SET THIS to where ever you put the data
datadir="..."

# figure 4 (35 -> NM)
python epicast_heterogeneity.py 4 "$datadir"

# figure 5
python epicast_heterogeneity.py 5 "$datadir"

# figure 6
python epicast_heterogeneity.py 6 "$datadir"

# figure S1
python epicast_heterogeneity.py s1 "$datadir"

# figure S2
python epicast_heterogeneity.py s2 "$datadir"

# figure S3 (12 -> FL)
python epicast_heterogeneity.py s3 "$datadir"
```

***Note***: all figures are saved as PDFs in the directory in which the main script is located

Additional scripts:

1. high_risk_industries_check.py: determines high risk industries used in analysis

2. sg_size_analysis.py: interprets schoolgroup size checkpoint file; outputs schoolgroup size data for figure 6 generation

Experimental files (TOML files for epicast initialisation; bash scripts for running EpiCast) are included. 

## Model variants
Model 1 - random NAICS code assignment:

src/initialize/urbanpop_v2/communities.h:

    static uint16_t naics_set[97] = {full urbanpop NAICS code list}

    static inline int get_naics_code(const agent_t* agent, int UNUSED(emp_status)) {
        return  naics_set[lrand48() % 97];
    }

Model 2 - static schoolgroup size target:

src/initialize/community_groups.h:

    sg->count[0] = round_max(sg->count[0], 18.44, 1); # FL state average
    sg->count[1] = round_max(sg->count[1], 18.44, 1); # FL state average
    ...
