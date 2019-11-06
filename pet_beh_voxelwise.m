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
% Comma separated values (csv) file containing the behavioral covariate
% data in columns. Include one column with subject numbers.
% 
% OUTPUT
% Estimated GLM, beta and contrast images and Tmaps per covariate. Binary
% images of the significant voxels in each contrast at specified thresholds
% are written too if requested.
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
spm pet

iostruct = struct('projectDir', projectDir, ...
                  'derivBEHdir', fullfile(projectDir, 'bids','derivatives','beh', '<task_name>'), ...
                  'derivPETdir', fullfile(projectDir, 'bids','derivatives','pet'), ...
                  'codeDir', fullfile(projectDir,'code','analysis','beh','<task_name>','pet'), ...
                  'batchDir', fullfile(projectDir,'code','analysis','beh','<task_name>','pet'));

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
todo.glmType = 'multi';

% Define csv file with behavioral data to be used as covariates
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
% Type            : 'uncorrected', 'fwe'
% Threshold       : cell array of thresholds for multiple outputs
% includeNegative : include the inverse (negative) contrast?
% exportBinary    : export a binary nifti file for each contrast with the 
%                   significant voxels in that contrast?
%                   Indicate whether the inverse (negative) contrast and 
%                   the combination of positive and negative should also be
%                   exported.
%                   includeNegativeContrast must be true in order to export
%                   the negative and combined binary images.

results.thresholdType             = 'uncorrected';
results.threshold                 = {0.001};
results.includeNegativeContrast   = false;
results.exportBinary              = true;
results.exportBinaryNegative      = false;
results.exportBinaryCombined      = false;

% Use inclusive mask?
results.mask.use      = false;
results.mask.image    = fullfile(iostruct.codeDir,'single_subj_T1_brain_mask.nii');

% Open SPM GUI? 
gui = false;

% Combine everything in settings structure
% -------------------------------------------------------------------------
settings = struct('io',iostruct, ...
                  'subjectIx',subjectIx, ...
                  'todo',todo, ...
                  'covariate',covariate, ...
                  'results',results, ...
                  'gui',gui);

% =========================================================================
% END USER INPUT
% =========================================================================

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
    behdata  = dataset('file', behDataFile, 'delim', ',','TreatAsEmpty','NA');
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

            if settings.covariate.(covname)
                
                % Output Directory
                % ---------------------------------------------------------
                outDir = fullfile(settings.io.derivBEHdir, 'pet', covname);

                if settings.todo.DesignSpecification
                    
                    % Make output directory if it does not exist
                    % -----------------------------------------------------
                    if ~exist(outDir,'dir')
                        mkdir(outDir)
                    end
            
                    % Get covariate vector
                    % -----------------------------------------------------
                    cov.subs = behdata.subject;
                    cov.vec  = behdata.(covname);

                    % Keep only cov data for subjects that have a ki img
                    % .....................................................
                    ix_tokeep = find(ismember(cov.subs, ki_imgs.subs));
                    cov.subs  = cov.subs(ix_tokeep);
                    cov.vec   = cov.vec(ix_tokeep);

                    % Remove subjects with missing data in covariate vector
                    % .....................................................
                    ix_nans = find(isnan(cov.vec));
                    cov.subs(ix_nans)       = [];
                    cov.vec(ix_nans)        = [];
                    
                    % Keep only Ki data for those subjects that have 
                    % covariate data
                    % .....................................................
                    ix_tokeep = find(ismember(ki_imgs.subs,cov.subs));
                    ki_imgs_to_include = ki_imgs.imgs(ix_tokeep);
                    
                    % Fill Batch and run
                    % -----------------------------------------------------
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
                    % -----------------------------------------------------
                    jobName = 'design_oneSttest_singleCovariate';
                    run_spm_jobs(jobName,dirs.jobs,jobs);
                end
                
                if settings.todo.ModelEstimation
                    
                    % Model Estimation
                    % -----------------------------------------------------
                    SPMmat = cellstr(fullfile(outDir,'SPM.mat'));
                    model_estimation(SPMmat,dirs.jobs);
                    
                end
                if settings.todo.Contrasts_and_ResultsExport
                    
                    % Contrast and results export
                    % -----------------------------------------------------
                    SPMmat = cellstr(fullfile(outDir,'SPM.mat'));
                    contrasts_and_ResultsExport(settings,SPMmat,cellstr(covname),dirs.jobs);
                end 
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
                
                % Add covariate only if set to true in settings
                if settings.covariate.(covname)
                    
                    % Save covariate name in includedRegressors file
                    fprintf(fid,[covname '\n']);
                    
                    % Add to covariates
                    % .....................................................
                    cov(ix).subs = behdata.subject;
                    cov(ix).vec  = behdata.(covname);
                    cov(ix).name = covname;
                    
                    ix = ix+1;
                end
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
            SPMmat = cellstr(fullfile(outDir, 'SPM.mat'));
            model_estimation(SPMmat,dirs.jobs);
                    
        end
        
        if settings.todo.Contrasts_and_ResultsExport
            
            % Contrast and results export
            % -------------------------------------------------------------
            SPMmat = cellstr(fullfile(outDir, 'SPM.mat'));
            
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
% SPMmat : cellstr; path to SPM.mat file of model to be estimated
% jobDir : char; directory where spm job file will be saved
% 
% =========================================================================
% 

% Fill in batch
clear jobs
jobs{1}.spm.stats.fmri_est.spmmat = SPMmat;
jobs{1}.spm.stats.fmri_est.write_residuals = 0;
jobs{1}.spm.stats.fmri_est.method.Classical = 1;

% Save and run job
jobName = 'estimate';
run_spm_jobs(jobName,jobDir,jobs);
end

function contrasts_and_ResultsExport(settings, SPMmat, covnames, jobDir)

% INPUTS
% settings  : struct with user specified settings
% SPMmat    : cellstr; path to SPM.mat file
% covnames  : cellstr; names of covariates to make contrasts for
% jobDir    : char; path to directory where spm job file is saved
% 
% =========================================================================
% 

% Loop in case multiple results exports at different p-thresholds are 
% requested
% -------------------------------------------------------------------------
for ip = 1:numel(settings.results.threshold)
    
    % Base contrast weights vector of zero for all regressors
    baseCon = zeros(1,numel(covnames)+1);
    
    % Prevent erroneous accumulation of jobs
    clear jobs

    % Create Contrasts
    % -----------------------------------------------------------------
    jobs{1}.spm.stats.con.spmmat = SPMmat;

    % Loop over requested covariates
    for icov = 1:numel(covnames)
        
        % Define covariate contrast
        jobs{1}.spm.stats.con.consess{icov}.tcon.name     = covnames{icov};
        
        % Index of cov regressor is cov number + 1 for the constant column
        weights         = baseCon;
        weights(icov+1) = 1;
        
        jobs{1}.spm.stats.con.consess{icov}.tcon.weights  = weights;
        jobs{1}.spm.stats.con.consess{icov}.tcon.sessrep  = 'none';
    end
    
    if settings.results.includeNegativeContrast
        
        % Add inverse contrasts at the end
        for icov = 1:numel(covnames)
            
            % Define negative covariate contrast
            % .............................................................
            % Append to list of contrasts
            conIx = numel(jobs{1}.spm.stats.con.consess) + 1;
            jobs{1}.spm.stats.con.consess{conIx}.tcon.name     = ['negative_' covnames{icov}];
            
            % Ix of cov regressor is cov number + 1 for the constant column
            weights         = baseCon;
            weights(icov+1) = -1;
            
            jobs{1}.spm.stats.con.consess{conIx}.tcon.weights  = weights;
            jobs{1}.spm.stats.con.consess{conIx}.tcon.sessrep  = 'none';
        end
    end

    % Delete existing contrasts
    jobs{1}.spm.stats.con.delete = 1;
       
    % ---------------------------------------------------------------------
    % Export Result
    % ---------------------------------------------------------------------
    if settings.results.exportBinary
      
        % Loop over contrasts and add export jobs to joblist
        for icon = 1:numel(jobs{1}.spm.stats.con.consess)
            
            % Path to SPM, and empty titlestring
            jobs{icon+1}.spm.stats.results.spmmat = SPMmat;
            jobs{icon+1}.spm.stats.results.conspec.titlestr = '';

            % Contrast index
            jobs{icon+1}.spm.stats.results.conspec.contrasts = icon;

            % Threshold type, value and min cluster size
            if strcmpi(settings.results.thresholdType,'uncorrected')
                jobs{icon+1}.spm.stats.results.conspec.threshdesc = 'none';
            elseif strcmpi(settings.results.thresholdType,'fwe')
                jobs{icon+1}.spm.stats.results.conspec.threshdesc = 'FWE';
            end
            jobs{icon+1}.spm.stats.results.conspec.thresh = settings.results.threshold{ip};
            jobs{icon+1}.spm.stats.results.conspec.extent = 0;

            % Export as binary. Define basename. Use inclusive mask if 
            % applicable
            % .............................................................
            jobs{icon+1}.spm.stats.results.units = 1;

            % Get p as string for use in basename
            p = regexp(num2str(settings.results.threshold{ip}), '\.', 'split');
            p = p{2};

            % Inclusive mask or not
            if settings.results.mask.use
                jobs{icon+1}.spm.stats.results.conspec.mask.image.name  = cellstr(settings.results.mask.image);
                jobs{icon+1}.spm.stats.results.conspec.mask.image.mtype = 0;
            else
                jobs{icon+1}.spm.stats.results.conspec.mask.none = 1;
            end

            % Output base name (SPM prepends spmT_xxx to the name)
            % If icon > numel(covnames) it's a negative contrast
            if ~(icon > numel(covnames))
                jobs{icon+1}.spm.stats.results.export{1}.binary.basename = sprintf('significant_voxels_%s_%s_p%s',covnames{icon},settings.results.thresholdType,p);
            elseif icon > numel(covnames)
                jobs{icon+1}.spm.stats.results.export{1}.binary.basename = sprintf('negativeCon_significant_voxels_%s_%s_p%s',covnames{icon-numel(covnames)},settings.results.thresholdType,p);
            end
        end
        
        % If requested, create combined binary mask of significant clusters
        % of boththe positive and negative covariate correlation contrast
        % -----------------------------------------------------------------
        if settings.results.exportBinaryCombined
    
            % Output dir
            [outDir,~,~] = fileparts(SPMmat{1});
            
            jobNr = numel(jobs) + 1;
            
            for icov = 1:numel(covnames)
                
                % Find images. Skip if one of the two does not exist
                pos = spm_select('FPList',outDir,sprintf('^spmT_.*_significant_voxels_%s_%s_p%s.nii$',covnames{icov},settings.results.thresholdType,p));
                neg = spm_select('FPList',outDir,sprintf('^spmT_.*_negativeCon_significant_voxels_%s_%s_p%s.nii$',covnames{icov},settings.results.thresholdType,p));
            
                if isempty(pos) || isempty(neg)
                    continue
                end
                
                % Output file name
                outputname = sprintf('combPosNeg_significant_voxels_%s_%s_p%s',covnames{icov},settings.results.thresholdType,p);
                
                % Fill imcalc job
                jobs{jobNr}.spm.util.imcalc.input = cellstr(strvcat(pos,neg));
                jobs{jobNr}.spm.util.imcalc.output = outputname;
                jobs{jobNr}.spm.util.imcalc.outdir = cellstr(outDir);
                jobs{jobNr}.spm.util.imcalc.expression = 'i1 | i2';
                jobs{jobNr}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                jobs{jobNr}.spm.util.imcalc.options.dmtx = 0;
                jobs{jobNr}.spm.util.imcalc.options.mask = 0;
                jobs{jobNr}.spm.util.imcalc.options.interp = 0;
                jobs{jobNr}.spm.util.imcalc.options.dtype = 2;
                
                jobNr = jobNr + 1;
            end
        end
    end
    
    % Save and run job
    % =============================================================
    jobName = 'Contrast_and_ResultsExport';
    run_spm_jobs(jobName,jobDir,jobs);
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
