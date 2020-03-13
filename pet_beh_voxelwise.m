% pet_beh_voxelwise
%
% DESCRIPTION
% Code for doing a voxelwise correlation analysis of some score or scores
% with dopamine synthesis capacity.
% 
% The user indicates which variable(s) to use as covariate(s), along with
% other settings. Depending on the option selected a multiple regression
% analysis is performed including all covariates selected, or individual
% single regressions are run for each selected covariate separately.
% 
% Test type: one sample t-tests
%
% INPUT
% A separated values (csv or tsv) file containing the behavioral covariate
% data in columns. Include one column with subject numbers.
% 
% OUTPUT
% Estimated GLM, beta and contrast images and Tmaps per covariate. Binary
% images of the significant voxels in each contrast at specified thresholds
% are written too if requested.
% In the specified behavioral derivatives direcotry (derivBEHdir in i/o 
% structure below in user input) a new directory "pet/<covariateName>" is 
% created, in which the output is saved.
% 
% -------------------------------------------------------------------------
% Ruben van den Bosch
% Donders Institute for Brain, Cognition and Behaviour
% Radboud University
% Nijmegen, The Netherlands
% April 2019
%

% USER INPUT
% =========================================================================
% Directory structure
% -------------------------------------------------------------------------
if isunix
    projectDir = '';
    addpath /home/common/matlab/spm12/
elseif ispc
    projectDir = '';
    addpath H:\common\matlab\spm12
else
    error('Unknown OS');
end

iostruct = struct('projectDir', projectDir, ...
                  'derivBEHdir', fullfile(projectDir, 'bids','derivatives','beh', '<task_name>'), ...
                  'derivPETdir', fullfile(projectDir, 'bids','derivatives','pet'), ...
                  'codeDir', fullfile(projectDir,'code','analysis','beh','<task_name>','pet'), ...
                  'functionDir', fileparts(mfilename('fullpath'))); % dir containing required function create_significant_voxels_binary (included in repo; now set to this mfile's dir)

% Subject to run
% -------------------------------------------------------------------------
subjectIx = 1:100;

% Steps to do
% -------------------------------------------------------------------------
todo.DesignSpecification            = true;
todo.ModelEstimation                = true;
todo.Contrasts_and_ResultsExport    = true;

% Single or multiple regression?
% -------------------------------------------------------------------------
% Do individual regressions for each predictor separately or do a multiple
% regression with all predictors in one model?
% 
% todo.glmType : char; 'single' OR 'multi'
todo.glmType = 'single';

% Define csv or tsv file with behavioral data to be used as covariates
% -------------------------------------------------------------------------
% Make sure the column for subject numbers is named "subject"
behDataFile = fullfile(iostruct.derivBEHdir, '<file_name>.csv');

% Which scores to use as covariate
% -------------------------------------------------------------------------
% Make sure the score names match the column headers from the csv file with
% covariate data (the code loops over the fieldnames of the covariate
% structure)
%
% Syntax:
% covariate.<score1>    = true;
% covariate.<score2>    = false

covariate.EXAMPLE        = false;

% Settings for results to export
% -------------------------------------------------------------------------
% exportBinary    : export a binary nifti file for each contrast with the 
%                   significant voxels in that contrast?
%                   This binary image contains the significant voxels for
%                   effects in both the positive and negative direction. In
%                   order to create this image, a negative (inverse)
%                   contrast is created for each requested regular
%                   contrast.
% thresholdType   : 'uncorrected' OR 'fwe'
% threshold       : cell array of p-value thresholds for multiple outputs
% extent          : double; minumum cluster size
results.exportBinary                = true;
results.significance.thresholdType  = 'uncorrected';
results.significance.threshold      = {0.001};
results.significance.extent         = 0;

% Open SPM GUI? 
gui = false;

% =========================================================================
% END USER INPUT
% =========================================================================

% Combine user input in settings structure
% -------------------------------------------------------------------------
settings = struct('io',iostruct, ...
                  'subjectIx',subjectIx, ...
                  'todo',todo, ...
                  'covariate',covariate, ...
                  'results',results, ...
                  'gui',gui);

% Add function dir to matlab path
addpath(settings.io.functionDir);

% Load SPM PET defaults and open GUI if set to true
spm('defaults','PET');
if settings.gui
    spm pet
else
    spm_jobman('initcfg');
end

% Define jobdir
dirs.jobs = fullfile(settings.io.codeDir, 'jobs');

% Make dirs if they don't exist
dirNames = fieldnames(dirs);
for iDir = 1:numel(dirNames)
    if ~exist(dirs.(dirNames{iDir}),'dir')
        mkdir(dirs.(dirNames{iDir}));
    end
end

% Get covariate names
covs = fieldnames(settings.covariate);

% If design specification is set to true, pre-collect available ki images
% to prevent redoing this on every loop in case of individual regressions.
% Also preload the behavior data file.
% -------------------------------------------------------------------------
if settings.todo.DesignSpecification
    
    % Collect smoothed normalized Ki images
    % ---------------------------------------------------------------------
    % Log for which subjects ki imgs are available
    ki_imgs.subs = [];
    ki_imgs.imgs = '';

    for ix = 1:numel(settings.subjectIx)
        iSubject = settings.subjectIx(ix);

        % If exists, add smoothed normalized Ki img and log subject num
        % .................................................................
        img = fullfile(settings.io.derivPETdir, sprintf('sub-%.3d',iSubject), 'Ki', 'MRIspace', 'swKi_map_brain.nii');

        if exist(img, 'file')
            ki_imgs.subs = [ki_imgs.subs; iSubject];
            ki_imgs.imgs = strvcat(ki_imgs.imgs, img);
        end
    end
    ki_imgs.imgs = cellstr(ki_imgs.imgs);
    
    % Load behavior data file
    % ---------------------------------------------------------------------
    [~,~,ext] = fileparts(behDataFile);
    if strcmp(ext,'.csv')
        delim = ',';
    elseif strcmp(ext,'.tsv')
        delim = '\t';
    end
    behdata  = readtable(behDataFile,'Delimiter',delim,'TreatAsEmpty','NA');
end

% Switch between individual regressions for each predictor or one multiple
% regression with all predictors
% -------------------------------------------------------------------------
% If only one covariate is set to true, set glmType to single
% .........................................................................
TrueFalse = cell2mat(struct2cell(settings.covariate));
if numel(TrueFalse(TrueFalse==1)) == 1 && strcmpi(settings.todo.glmType,'multi')
    warning('glmType is set to multi, but only one covariate is set to true. Changing glmType to single')
    settings.todo.glmType = 'single';
elseif numel(TrueFalse(TrueFalse==1)) == 0
    error('No covariate has been set to true for processing')
end

switch lower(settings.todo.glmType)
    case 'single'

        % In case single regressions for each individual predictor, loop
        % over predictors and do design specification, estimation and
        % contrast and results export in each loop.
        % =================================================================
        for icov = 1:numel(covs)
            covname = covs{icov};

            % Skip this covariate if it's not set to true in settings
            if ~settings.covariate.(covname)
                continue
            end
            
            % Output Directory
            % -------------------------------------------------------------
            outDir = fullfile(settings.io.derivBEHdir, 'pet', covname);

            if settings.todo.DesignSpecification

                % Make output directory if it does not exist
                % ---------------------------------------------------------
                if ~exist(outDir,'dir')
                    mkdir(outDir)
                end

                % Get covariate vector
                % ---------------------------------------------------------
                cov.subs = behdata.subject;
                cov.vec  = behdata.(covname);

                % Keep only cov data for subjects that have a ki img
                % .........................................................
                ix_tokeep = find(ismember(cov.subs, ki_imgs.subs));
                cov.subs  = cov.subs(ix_tokeep);
                cov.vec   = cov.vec(ix_tokeep);

                % Remove subjects with missing data in covariate vector
                % .........................................................
                ix_nans = find(isnan(cov.vec));
                cov.subs(ix_nans)       = [];
                cov.vec(ix_nans)        = [];

                % Keep only Ki data for those subjects that have covariate 
                % data
                % .........................................................
                ix_tokeep = find(ismember(ki_imgs.subs,cov.subs));
                ki_imgs_to_include = ki_imgs.imgs(ix_tokeep);

                % Fill Batch and run
                % ---------------------------------------------------------
                clear jobs
                jobs{1}.spm.stats.factorial_design.dir                       = cellstr(outDir);
                jobs{1}.spm.stats.factorial_design.des.t1.scans              = ki_imgs_to_include;
                jobs{1}.spm.stats.factorial_design.cov.c                     = cov.vec;
                jobs{1}.spm.stats.factorial_design.cov.cname                 = covname;
                jobs{1}.spm.stats.factorial_design.cov.iCFI                  = 1;
                jobs{1}.spm.stats.factorial_design.cov.iCC                   = 1;
                jobs{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
                jobs{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
                jobs{1}.spm.stats.factorial_design.masking.im                = 1;
                jobs{1}.spm.stats.factorial_design.masking.em                = {''};
                jobs{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
                jobs{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
                jobs{1}.spm.stats.factorial_design.globalm.glonorm           = 1;

                % Save and run job
                % ---------------------------------------------------------
                jobName = 'design_oneSttest_singleCovariate';
                run_spm_jobs(jobName,dirs.jobs,jobs);
            end

            if settings.todo.ModelEstimation

                % Model Estimation
                % ---------------------------------------------------------
                SPMmat = fullfile(outDir,'SPM.mat');
                model_estimation(SPMmat,dirs.jobs);

            end
            if settings.todo.Contrasts_and_ResultsExport

                % Contrast and results export
                % ---------------------------------------------------------
                SPMmat = fullfile(outDir,'SPM.mat');
                contrasts_and_ResultsExport(settings,SPMmat,cellstr(covname),dirs.jobs);
            end
        end
    case 'multi'
        
        % In case of a multiple regression with all predictors in one
        % model, loop over predictors in design specification stage to
        % collect them in one model. Then run model estimation and
        % contrast and result export once.
        % =================================================================
        
        % Output directory
        % -----------------------------------------------------------------
        outDir = fullfile(settings.io.derivBEHdir, 'pet', 'multiple');
        
        if settings.todo.DesignSpecification

            % If the output directory already exists, ask what to do
            % -------------------------------------------------------------
            if exist(outDir,'dir')
                answer = questdlg('The output directory called multiple already exists. Do you want to create a new directory with a tracking number (e.g. multi_3) or overwrite existing directory?', ...
                                   'Output directory exists', ...
                                   'New directory', 'Cancel', 'Overwrite', 'New directory');
                switch answer
                    case 'New directory'
                        num_existing = numel(dir([outDir '*']));
                        outDir = sprintf('%s_%d',outDir,num_existing+1);
                        mkdir(outDir);
                    case 'Overwrite'
                        warning('Output directory will be overwritten: %s',outDir)
                        if exist(fullfile(outDir,'SPM.mat'),'file')
                            delete(fullfile(outDir,'SPM.mat'))
                        end
                    case 'Cancel'
                        return
                end
            else
                % Else, make the output directory
                mkdir(outDir);
            end
            
            % Get covariate vectors
            % -------------------------------------------------------------
            % Save which covariates are included in a text file
            fid = fopen(fullfile(outDir,'includedRegressors.txt'),'w');
            
            % Loop over covariates
            ix = 1;
            for icov = 1:numel(covs)
                covname = covs{icov};
                
                % Skip covariate if not set to true in settings
                if ~settings.covariate.(covname)
                    continue
                end
                
                % Save covariate name in includedRegressors file
                fprintf(fid,[covname '\n']);

                % Add to covariates
                % .....................................................
                cov(ix).subs = behdata.subject;
                cov(ix).vec  = behdata.(covname);
                cov(ix).name = covname;

                ix = ix+1;
            end
            
            % Close includedRegressors file
            fclose(fid);
            
            % Remove missing values and corresponding subjects from 
            % covariate data
            % .............................................................
            for icov = 1:numel(cov)
                cov(icov).subs(isnan(cov(icov).vec)) = [];
                cov(icov).vec(isnan(cov(icov).vec))  = [];
            end
            
            % Keep only subjects that have no missing data in any of the
            % covariates
            % .............................................................
            if numel(cov) == 2
                includeSubs = intersect(cov(1).subs,cov(2).subs);
            elseif numel(cov) > 2
                cmd = 'intersect(cov(1).subs,cov(2).subs)';
                for icov = 3:numel(cov)
                    cmd = sprintf('intersect(%s,cov(%d).subs)',cmd,icov);
                end
                includeSubs = eval(cmd);
            end
            
            % Keep only subjects for whom ki data is available
            % .............................................................
            includeSubs = includeSubs(ismember(includeSubs, ki_imgs.subs));
            
            % Remove not selected subjects and cov data
            % .............................................................
            for icov = 1:numel(cov)
                ix_tokeep = find(ismember(cov(icov).subs,includeSubs));
                cov(icov).subs = cov(icov).subs(ix_tokeep);
                cov(icov).vec  = cov(icov).vec(ix_tokeep);
            end
            
            % Keep only Ki data for those subjects that have covariate data
            % .............................................................
            ix_tokeep = find(ismember(ki_imgs.subs,includeSubs));
            ki_imgs_to_include = ki_imgs.imgs(ix_tokeep);
            
            % Fill batch and run
            % -------------------------------------------------------------
            clear jobs
            jobs{1}.spm.stats.factorial_design.dir                       = cellstr(outDir);
            jobs{1}.spm.stats.factorial_design.des.t1.scans              = ki_imgs_to_include;
            for icov = 1:numel(cov)
                jobs{1}.spm.stats.factorial_design.cov(icov).c           = cov(icov).vec;
                jobs{1}.spm.stats.factorial_design.cov(icov).cname       = cov(icov).name;
                jobs{1}.spm.stats.factorial_design.cov(icov).iCFI        = 1;
                jobs{1}.spm.stats.factorial_design.cov(icov).iCC         = 1;
            end
            jobs{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
            jobs{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
            jobs{1}.spm.stats.factorial_design.masking.im                = 1;
            jobs{1}.spm.stats.factorial_design.masking.em                = {''};
            jobs{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
            jobs{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
            jobs{1}.spm.stats.factorial_design.globalm.glonorm           = 1;

            % Save and run job
            % -------------------------------------------------------------
            jobName = 'design_oneSttest_multiCovariates';
            run_spm_jobs(jobName,dirs.jobs,jobs);
        end
        
        if settings.todo.ModelEstimation
            
            % Model Estimation
            % -------------------------------------------------------------
            SPMmat = fullfile(outDir, 'SPM.mat');
            model_estimation(SPMmat,dirs.jobs);
                    
        end
        
        if settings.todo.Contrasts_and_ResultsExport
            
            % Contrast and results export
            % -------------------------------------------------------------
            SPMmat = fullfile(outDir, 'SPM.mat');
            
            % Covariates to include
            covnames = [];
            for icov = 1:numel(covs)
                if settings.covariate.(covs{icov})
                    covnames = strvcat(covnames,covs{icov});
                end
            end
            covnames = cellstr(covnames);
            
            % Create contrasts and export results
            contrasts_and_ResultsExport(settings,SPMmat,covnames,dirs.jobs);
        end
end

function model_estimation(SPMmat,jobDir)
% model_estimation(SPMmat)
% 
% Run SPM model estimation
% 
% INPUT
% SPMmat : char; path to SPM.mat file of model to be estimated
% jobDir : char; directory where spm job file will be saved
% 
% =========================================================================
% 

% Fill in batch
clear jobs
jobs{1}.spm.stats.fmri_est.spmmat           = cellstr(SPMmat);
jobs{1}.spm.stats.fmri_est.write_residuals  = 0;
jobs{1}.spm.stats.fmri_est.method.Classical = 1;

% Save and run job
jobName = 'estimate';
run_spm_jobs(jobName,jobDir,jobs);
end

function contrasts_and_ResultsExport(settings, SPMmat, covnames, jobDir)

% SUMMARY
% Function to create the contrasts for the covariates set to true in
% settings, and if requested export binary images of the significant voxels
% at the specified significance thresholds (export uses external function).
% 
% INPUTS
% settings  : struct with user specified settings
% SPMmat    : char; path to SPM.mat file
% covnames  : cellstr; names of covariates to make contrasts for
% jobDir    : char; path to directory where spm job file is saved
% 
% =========================================================================
% 

% Create contrasts
% -------------------------------------------------------------------------
% Prevent erroneous accumulation of jobs
clear jobs

% Base contrast weights vector of zeroes for all regressors (1 constant +
% all covariates)
baseCon = zeros(1,numel(covnames)+1);

% Fill job
% .........................................................................
% SPM.mat file
jobs{1}.spm.stats.con.spmmat = cellstr(SPMmat);

% Loop over requested covariates
for icov = 1:numel(covnames)

    % Covariate contrast name
    jobs{1}.spm.stats.con.consess{icov}.tcon.name     = covnames{icov};

    % Contrast weight of 1
    % Index of cov regressor is cov number + 1 for the constant column
    weights         = baseCon;
    weights(icov+1) = 1;
    jobs{1}.spm.stats.con.consess{icov}.tcon.weights  = weights;
    
    % Don't repeat over sessions
    jobs{1}.spm.stats.con.consess{icov}.tcon.sessrep  = 'none';
end

% Delete existing contrasts
jobs{1}.spm.stats.con.delete = 1;

% Save and run contrast job
% .........................................................................
jobName = 'create_contrasts';
run_spm_jobs(jobName,jobDir,jobs);

% Export significant voxels as binary images
% -------------------------------------------------------------------------
if settings.results.exportBinary
    % Call external function
    modality = 'pet';
    create_significant_voxels_binary(SPMmat,covnames,modality,settings.results.significance)
end
end

function run_spm_jobs(jobName,jobDir,jobs)
% RUN_SPM_JOBS(jobName,jobDir,jobs)
% 
% Subfunction to save and run spm jobs.
% 

% Filename for save job
jobFile = fullfile(jobDir,[jobName,'_',datestr(now,'yyyymmddTHHMMSS'),'.m']);

% Initialise job
jobId = cfg_util('initjob', jobs);

% If successful save and run job
sts = cfg_util('isjob_id', jobId);
if sts
    cfg_util('savejob', jobId, jobFile);
    cfg_util('run', jobId);
else
    error('Error in initialising %s job.',jobName)
end
end
