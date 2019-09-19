# GSEA R Implementation v1.2
## Revisions to the GSEA-P R Application (Subramanian and Tamayo et al 2005)

Preliminary revisions to the original GSEA implementation in R to enable it to run on modern R distributions.
Recommended R minimum version: 3.6

Additionally, a rudimentary Collapse dataset function has been backported from the Java GSEA application, however differences in the implementation result in inconsistencies with the desktop collapse function.

GSEA-R v1.2 contains a preliminary implementation of GSEA-preranked. This implementation has not been thoroughly tested and is reliant on undocumented gene set permutation code.

The collapse dataset feature requires the <tt>dplyr</tt> package from <tt>tidyverse</tt>

You may install this requirement from CRAN with <tt>install.packages("dplyr")</tt>

The original R-GSEA documentation indicates that it supports "phenotype permutation" mode only, however, code was present in the application to perform gene set permutation testing. This code has been enabled, but not tested.

R-GSEA remains unsupported by the GSEA-MSigDB Team.

See the original R-GSEA readme here: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/R-GSEA_Readme

### Running this application
Install the GSEA package with: <tt>library("devtools")</tt> then <tt>install_github("GSEA-MSigDB/GSEA_R")</tt>

Optionally, a helper script has been provided to simplify use of this package. Initialize the helper script with the <tt>source()</tt> command by calling <tt>source(system.file('extdata', 'Run.GSEA.R', package = 'GSEA'))</tt>

The helper script will import the GSEA functions from the <tt>GSEA</tt> package, and begin a series of prompts for loading files required for GSEA and setting experiment parameters.
