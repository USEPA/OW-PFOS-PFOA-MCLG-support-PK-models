# Repository for PFOS/PFOA MCLG support PK modeling work
This repository contains the pharmacokinetic modeling code for calculating human equivalent doses from the Proposed Maximum Contaminant Level Goal for Perfluorooctanoic Acid (PFOA) in Drinking Water and the Proposed Maximum Contaminant Level Goal for Perfluorooctane Sulfonic Acid (PFOS) in Drinking Water (Proposed MCLGs for PFOA and PFOS in Drinking Water).

Because this analysis directly generates the numbers used in [Proposed MCLGs for PFOA and PFOS in Drinking Water](https://www.regulations.gov/docket/EPA-HQ-OW-2022-0114), we ask that you visit [CONTRIBUTING](https://github.com/USEPA/OW-PFOS-PFOA-MCLG-support-PK-models/blob/master/contributing.md) for instructions on how to contact EPA about code requests.

## Dependencies
In order to run this code, the following dependencies are needed.

### R (https://cran.r-project.org/)
deSolve: `install.packages("deSolve")`

RTools: https://cran.r-project.org/bin/windows/Rtools/. We suggest RTools v.4.0 for maximum compatibility with current code.

### RMCSim

The RMCSim pseudo package (provided in this repo) must be installed through R. This can be accomplished using the RMCSim subdirectory: 

In R:

```
setwd("/path/to/OW-PFOS-PFOA-MCLG-support-PK-models")
install.packages("RMCSim", repo=NULL, type="source")`
```

For additional examples on the use of RMCSim with PK modeling, please see [Kapraun et al., 2022](https://academic.oup.com/toxsci/article/189/2/155/6661359)

### Python (https://www.anaconda.com/products/distribution)
rpy2: `pip install rpy2==3.4.4`

Note: Must use `pip` and **not** conda distribution. Also, must use rpy2 v.3.4.4

Note: If import of rpy2 through `from rpy2.robjects import r` fails during import of `PFAS_DR.py`, add the following lines to `PFAS_DR.py` before `import rpy2`:

```
import os
os.environ['R_HOME'] = 'path/to/R/install'
```

Generally, the path to the R install is something like `C:\Program Files\R\R-<version-number>`

watermark: `pip install watermark`

### GNU MCSim (https://www.gnu.org/software/mcsim/)
The pharmacokinetic models used in this analysis are written in the GNU MCSim language (https://www.gnu.org/software/mcsim/) Bois F., 2009, GNU MCSim: Bayesian statistical inference for SBML-coded systems biology models Bioinformatics, 25:1453-1454, doi: 10.1093/bioinformatics/btp162 and as stated on the website "The basic tools needed to build and run GNU MCSim are not available to many users of Windows systems ... [but] R software when installed with its Rtools can compile and run GNU MCSim models." See [LICENSE](https://github.com/USEPA/OW-PFOS-PFOA-MCLG-support-PK-models/blob/master/license.md) for details on the GNU General Public License.

To aid in reproducibility, pertinent code from MCSim has been compiled to an executable (`mod.exe`) with `mod.exe` installed into `R/<R-version>/library/RMCSim/executables/mod.exe` where `R/<R-version>/library` can be found be running `.libPaths()` in R. For a more detailed description of using MCSim c code with RMCSim, please see install instructions from Kapraun, DF et al. A Generic Pharmacokinetic Model for Quantifying Mother-to-Offspring Transfer of Lipophilic Persistent Environmental Chemicals. 2022, Tox. Sci. 189(2), p. 155-174.

## Compiling a model
Before running code, the .model files must be compiled. This is best done in RStudio once all the R dependencies are satisfied.

### Compile model in R
Open RStudio and run the following lines

```
library(RMCSim)
Sys.setenv(PATH = paste("C:/Rtools40/usr/bin", Sys.getenv("PATH"), sep=";")) # Check RTools path. usr/bin of RTools must be added here
setwd('/path/to/model')
compile_model('model_name')
```

### Compile model in Python
Ensure the working directory is `OW-PFOS-PFOA-MCLG-support-PK-models\animal` and then run

```
from PFAS_DR import PFAS_DR
model_path = 'path/to/model_file.model'
model_compile = PFAS_DR()
model_compile.pk_run.compile_pfas_model(model_path)
```

## Pharmacokinetic Analysis

### Animal studies
Python and R are needed to run the animal analysis.

**Note**: Before running, check path for `Sys.setenv(PATH = paste("C:/Rtools40/usr/bin", Sys.getenv("PATH"), sep=";"))` in the `_define_R_code` module. Update path if needed.

The main class for animal simulations is housed in animal/PFAS_DR.py. This class is called by the various notebooks to complete animal-specific simulations.

Descriptions of each notebook and output files for the animal analysis are found in animal_README.txt

### Human analysis
Only R is needed to run the human analysis.

File descriptions for the human analysis are found human_README.txt

# Disclaimer
The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government. 
