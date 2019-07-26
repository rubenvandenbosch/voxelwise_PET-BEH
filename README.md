# voxelwise_PET-BEH

## DESCRIPTION
Code for doing a voxelwise correlation analysis of some score or scores
with dopamine synthesis capacity.

The user indicates which variable to use as covariate. A separate
analysis design is created and run for each covariate. 

Test type: one sample t-tests

INPUT  
Comma separated values (csv) file containing the behavioral covariate
data in columns.
 
OUTPUT  
Estimated GLM, beta and contrast images and Tmaps per covariate. Binary
images of the significant voxels in each contrast at specified thresholds
are written too. 
 
-------------------------------------------------------------------------
Ruben van den Bosch
Donders Institute for Brain, Cognition and Behaviour
Radboud University
Nijmegen, The Netherlands
April 2019

