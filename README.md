All code used to generate the results seen in the manuscript ____

# Dependencies
* minc-toolkit/1.9.17
* R/3.5.1
* RMINC/1.5.2.3
* Matlab/2020a

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

