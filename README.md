All code used to generate the results seen in the manuscript "High spatial overlap but diverging age-related trajectories of cortical magnetic resonance imaging markers aiming to represent intracortical myelin and microstructure"

# Dependencies
* minc-toolkit/1.9.17
* R/3.5.1
* RMINC/1.5.2.3
* Matlab/2020a

# Directory structure
AHBA -> Analysis relating gene-expression-derived cell-type densities with cortical MRI markers
- data_paper -> Vertex-wise gene-expression-derived cell-type densities

AIC_age_effect -> Analysis determining optimal shape of age trajectories with the Akaike Information Criterion (AIC) of each cortical MRI marker
- data_paper -> Vertex-wise optimal shape of age trajectories
- AIC_age_effect.R -> Analysis code

BigBrain -> Analysis relating histology-derived cell density with cortical MRI markers
- data_paper -> Vertex-wise histology-derived cell density

R1_spatial_distribution_average -> Analysis relating average R1 relaxation rates with cortical MRI markers
- data_paper -> Vertex-wise average R1 values
- R1_spatial_distribution_average.R -> Analysis code

lin_age_effect -> Linear age effect of each cortical MRI marker
- data_paper -> Vertex-wise beta coefficients of linear age effect
- lin_age_effect.R -> Analysis code

quad_age_effect -> Quadratic age effect of each cortical MRI marker
- data_paper -> Vertex-wise beta coefficients of quadratic age effect
- quad_age_effect.R -> Analysis code

spatial_distribution_average -> Computing spatial distribution average for each cortical MRI marker
- data_paper -> Vertex-wise beta average values for each cortical MRI marker
- spatial_distribution_average.R -> Analysis code

spin_test -> Spatial correlations and p-values for each analysis above using the spin test
- data -> Masks of midline vertices
- scripts -> All scripts for the spin test. Run the files starting with run_

master_anon.csv -> CSV with demographic, cognitive, and quality control (QC) data for all included subjects

# Steps to run code
1. Make sure all dependencies are loaded
2. Clone GitHub repository
```
git clone https://github.com/parent41/cortical_markers_paper
```
3. Download subject-wise vertex data for cortical MRI markers from the [Canadian Open Neuroscience Platform (CONP)](https://portal.conp.ca/index). Put all files in a directory named "vertex_files_20mm_anon" located in the main directory
4. Run R scripts to generate marker maps (for spatial distribution average, AIC of age effect shape, linear age effect, quadratic age effect)
```
# First set working directory to the specific analysis folder (e.g., spatial_distribution_average)
cd spatial_distribution_average
# Then run R script
Rscript spatial_distribution_average.R
cd ..

cd AIC_age_effect
Rscript AIC_age_effect.R
cd ..

cd lin_age_effect
Rscript lin_age_effect.R
cd ..

cd quad_age_effect
Rscript quad_age_effect.R
cd ..

cd R1_spatial_distribution_average
Rscript R1_spatial_distribution_average.R
cd ..
```
5. Run spin tests in Matlab (scripts to run begin with run_). Set working directory to ./spin_test/scripts/
```
matlab run_spatial_distribution_average.m
matlab run_R1_spatial_distribution_average.m
matlab run_lin_age_effect.m
matlab run_BigBrain.m
matlab run_AHBA.m
```

