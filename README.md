# voxelwise_PET-BEH

## DESCRIPTION
Code for doing a voxelwise correlation analysis of some score or scores
with dopamine synthesis capacity.

The user indicates which variable(s) to use as covariate(s), along with
other settings. Depending on the option selected a multiple regression
analysis is performed including all covariates selected, or individual
single regressions are run for each selected covariate separately.

Test type: one sample t-tests

INPUT
Comma separated values (csv) file containing the behavioral covariate
data in columns. Include one column with subject numbers.

OUTPUT
Estimated GLM, beta and contrast images and Tmaps per covariate. Binary
images of the significant voxels in each contrast at specified thresholds
are written too if requested.

-------------------------------------------------------------------------
Ruben van den Bosch
Donders Institute for Brain, Cognition and Behaviour
Radboud University
Nijmegen, The Netherlands
April 2019
