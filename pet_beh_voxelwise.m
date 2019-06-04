% pet_beh_voxelwise
%
% DESCRIPTION
% Code for doing a voxelwise correlation analysis of some score or scores
% with dopamine synthesis capacity.
% 
% The user indicates which variable to use as covariate. A separate
% analysis design is created and run for each covariate. 
% 
% Test type: one sample t-tests
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
                  'derivBEHdir', fullfile(projectDir, 'data','derivatives','beh', 'RL'), ...
                  'derivPETdir', fullfile(projectDir, 'bids','derivatives','pet'), ...
                  'codeDir', fullfile(projectDir,'code','analysis','beh','RL','pet'), ...
                  'batchDir', fullfile(projectDir,'code','analysis','beh','RL','pet'));

% Subject to run
% -------------------------------------------------------------------------
subjectIx = 1:100;

% Steps to do
% -------------------------------------------------------------------------
todo.DesignSpecification            = true;
todo.ModelEstimation                = true;
todo.Contrasts_and_ResultsExport    = true;

% Define csv file with behavioral data to be used as covariates
% -------------------------------------------------------------------------
% Make sure the column for subject numbers is named "subject"
behDataFile = fullfile(iostruct.derivBEHdir, 'allbeh_exOutliers_nativeKi.csv');

% Which scores to use as covariate
% -------------------------------------------------------------------------
% Make sure the score names match the column headers from the csv file with
% covariate data (the code loops over the fieldnames of the covariate
% structure)
%
% Syntax:
% covariate.<score1>    = true;
% covariate.<score2>    = false

covariate.subjective_value        = false;

% Settings for results to export
% -------------------------------------------------------------------------
% Type options : 'uncorrected', 'fwe'
% Threshold    : cell array of thresholds for multiple outputs
results.thresholdType = 'uncorrected';
results.threshold     = {0.001, 0.01};

% Use inclusive mask?
results.mask.use      = false;
results.mask.image    = '';

% Combine everything in settings structure
% -------------------------------------------------------------------------
settings = struct('io',iostruct, ...
                  'subjectIx',subjectIx, ...
                  'todo',todo, ...
                  'covariate',covariate, ...
                  'results',results);

% END USER INPUT
% =========================================================================

% Define jobdir
dirs.jobs = fullfile(settings.io.codeDir, 'jobs');

% Make dirs if they don't exist
dirNames = fieldnames(dirs);
for iDir = 1:numel(dirNames)
    if ~exist(getfield(dirs,dirNames{iDir}),'dir')
        mkdir(getfield(dirs,dirNames{iDir}));
    end
end

% Loop over covariates
covs = fieldnames(settings.covariate);

for icov = 1:numel(covs)
    
    covname = covs{icov};
    
    if settings.covariate.(covname)
        
        if settings.todo.DesignSpecification
            
        % =================================================================
        % Design Specification
        % =================================================================
        % Load job
        % -----------------------------------------------------------------
        jobTemplate = fullfile(settings.io.batchDir,'design_oneSttest_covariate.m');
        jobId       = cfg_util('initjob', jobTemplate);
        inputs      = cell(4,1);

        % Define inputs
        % =================================================================
        % 1. Output directory for SPM.mat
        % 2. Smoothed normalized Ki images
        % 3. Covariate vector
        % 4. Name of covariate

        % 1. Define output directory
        % -----------------------------------------------------------------
        outDir = fullfile(settings.io.derivBEHdir, 'pet', covname);
        if ~exist(outDir,'dir')
            mkdir(outDir)
        end

        % 2. Select images and log for which subjects ki imgs are available
        % -----------------------------------------------------------------
        ki_imgs.subs = [];
        ki_imgs.imgs = '';

        for ix = 1:numel(settings.subjectIx)
            iSubject = settings.subjectIx(ix);
            
            % If exists add smoothed normalized Ki img and log subject num
            % .............................................................
            img = fullfile(settings.io.derivPETdir, sprintf('sub-%.3d',iSubject), 'Ki', 'MRIspace', 'swKi_map.nii');
            
            if exist(img, 'file')
                ki_imgs.subs = [ki_imgs.subs; iSubject];
                ki_imgs.imgs = strvcat(ki_imgs.imgs, img);
            end
        end
        ki_imgs.imgs = cellstr(ki_imgs.imgs);
        
        % 3. Get covariate vector
        % -----------------------------------------------------------------
        behdata = dataset('file', RLdatafile, 'delim', ',','TreatAsEmpty','NA');
        cov.subs = behdata.subject;
        cov.vec = behdata.(covs{icov});

        % Keep only cov data for subjects that have a ki img
        % .................................................................
        ix_tokeep = find(ismember(cov.subs, ki_imgs.subs));
        cov.subs = cov.subs(ix_tokeep);
        cov.vec = cov.vec(ix_tokeep);
        
        % Remove subjects with missing data in covariate vector
        % .................................................................
        ix_nans = find(isnan(cov.vec));
        cov.subs(ix_nans)       = [];
        cov.vec(ix_nans)        = [];
        ki_imgs.subs(ix_nans)   = [];
        ki_imgs.imgs(ix_nans)   = [];
        
        % Fill in batch
        % -----------------------------------------------------------------
        inputs{1} = cellstr(outDir);
        inputs{2} = ki_imgs.imgs;
        inputs{3} = cov.vec;
        inputs{4} = covname;
     
        % Save and run job
        % =================================================================
        [~,nm,ext] = fileparts(jobTemplate);
        jobFile = fullfile(dirs.jobs,[nm,'_',covname,'_',datestr(now,'yyyymmddTHHMMSS'),ext]);

        sts = cfg_util('filljob', jobId, inputs{:});
        if sts
            cfg_util('savejob', jobId, jobFile);
            spm_jobman('run', jobFile, inputs{:});
        else 
            error('Error in design specification for covariate %s: job could not be filled',covname)
        end
        end

        if settings.todo.ModelEstimation
            
        % =================================================================
        % Model Estimation
        % =================================================================
        % Load job
        % -----------------------------------------------------------------
        jobTemplate = fullfile(settings.io.batchDir,'estimate.m');
        jobId       = cfg_util('initjob', jobTemplate);
        inputs      = cell(1,1);

        % Define inputs
        % =================================================================
        % 1. SPM.mat file containing the design matrix
        outDir = fullfile(settings.io.derivBEHdir, 'pet', covname);
        SPMmat = fullfile(outDir,'SPM.mat');

        % Fill in batch
        % -----------------------------------------------------------------
        inputs{1} = cellstr(SPMmat);

        % Save and run job
        % =================================================================
        [~,nm,ext] = fileparts(jobTemplate);
        jobFile = fullfile(dirs.jobs,[nm,'_',covname,'_',datestr(now,'yyyymmddTHHMMSS'),ext]);

        sts = cfg_util('filljob', jobId, inputs{:});
        if sts
            cfg_util('savejob', jobId, jobFile);
            spm_jobman('run', jobFile, inputs{:});
        else 
            error('Error in model estimation for contrast %s: job could not be filled',covname)
        end
        end
        
        if todo.Contrasts_and_ResultsExport
        
        % =================================================================
        % Create Contrasts And Export Result
        % =================================================================
        
        % Loop in case multiple results exports at different p-thresholds
        % are requested
        for ip = 1:numel(settings.results.threshold)
        
            % Create Contrasts
            % -------------------------------------------------------------
            % Define SPM
            outDir = fullfile(settings.io.derivBEHdir, 'pet', covname);
            SPMmat = fullfile(outDir,'SPM.mat');
            jobs{1}.spm.stats.con.spmmat = cellstr(SPMmat);
            
            % Define covariate contrast
            jobs{1}.spm.stats.con.consess{1}.tcon.name     = covname;
            jobs{1}.spm.stats.con.consess{1}.tcon.weights  = [0 1];
            jobs{1}.spm.stats.con.consess{1}.tcon.sessrep  = 'none';
            
            % Define negative covariate contrast
            jobs{1}.spm.stats.con.consess{2}.tcon.name     = ['negative_' covname];
            jobs{1}.spm.stats.con.consess{2}.tcon.weights  = [0 -1];
            jobs{1}.spm.stats.con.consess{2}.tcon.sessrep  = 'none';
                        
            % Delete existing contrasts
            jobs{1}.spm.stats.con.delete = 1;

            % Export Result
            % -------------------------------------------------------------
            % Loop twice, to export both positive and negative covariate
            % correlation contrast
            
            for icon = 1:2
                % Path to SPM, and empty titlestring
                jobs{icon+1}.spm.stats.results.spmmat = cellstr(SPMmat);
                jobs{icon+1}.spm.stats.results.conspec.titlestr = '';

                % Contrast index
                jobs{icon+1}.spm.stats.results.conspec.contrasts = icon;

                % Threshold type, value and min cluster size
                if strcmpi(settings.results.thresholdType,'uncorrected')
                    jobs{icon+1}.spm.stats.results.conspec.threshdesc = 'none';
                elseif strcmpi(settings.results.thresholdType,'fwe')
                    jobs{icon+1}.spm.stats.results.conspec.threshdesc = 'fwe';
                end
                jobs{icon+1}.spm.stats.results.conspec.thresh = settings.results.threshold{ip};
                jobs{icon+1}.spm.stats.results.conspec.extent = 0;

                % Export as binary. Define basename. Use inclusive mask if 
                % applicable
                % .........................................................
                jobs{icon+1}.spm.stats.results.units = 1;

                % Get p as string for use in basename
                p = regexp(num2str(settings.results.threshold{ip}), '\.', 'split');
                p = p{2};

                if ~settings.results.mask.use
                    jobs{icon+1}.spm.stats.results.conspec.mask.none = 1;
                    if icon == 1
                        jobs{icon+1}.spm.stats.results.export{1}.binary.basename = sprintf('sigClust_p%s_%s_binary',p,settings.results.thresholdType);
                    elseif icon == 2
                        jobs{icon+1}.spm.stats.results.export{1}.binary.basename = sprintf('negativeCon_sigClust_p%s_%s_binary',p,settings.results.thresholdType);
                    end
                else
                    jobs{icon+1}.spm.stats.results.conspec.mask.image.name  = cellstr(settings.results.mask.image);
                    jobs{icon+1}.spm.stats.results.conspec.mask.image.mtype = 0;
                    if icon == 1
                        jobs{icon+1}.spm.stats.results.export{1}.binary.basename = sprintf('sigClust_p%s_%s_binary_inclMask',p,settings.results.thresholdType);
                    elseif icon == 2
                        jobs{icon+1}.spm.stats.results.export{1}.binary.basename = sprintf('negativeCon_sigClust_p%s_%s_binary_inclMask',p,settings.results.thresholdType);
                    end
                end
            end

            % Create combined binary mask of significant clusters of both
            % the positive and negative covariate correlation contrast
            % -------------------------------------------------------------
            jobNr = numel(jobs) + 1;
            jobs{jobNr}.spm.util.imcalc.input = cellstr(strvcat(fullfile(outDir,['spmT_0001_' jobs{2}.spm.stats.results.export{1}.binary.basename '.nii']), ...
                                                                fullfile(outDir,['spmT_0002_' jobs{3}.spm.stats.results.export{1}.binary.basename '.nii'])));
            jobs{jobNr}.spm.util.imcalc.output = ['combPosNeg_' jobs{2}.spm.stats.results.export{1}.binary.basename];
            jobs{jobNr}.spm.util.imcalc.outdir = cellstr(outDir);
            jobs{jobNr}.spm.util.imcalc.expression = 'i1 | i2';
            jobs{jobNr}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            jobs{jobNr}.spm.util.imcalc.options.dmtx = 0;
            jobs{jobNr}.spm.util.imcalc.options.mask = 0;
            jobs{jobNr}.spm.util.imcalc.options.interp = 0;
            jobs{jobNr}.spm.util.imcalc.options.dtype = 2;

            % Save and run job
            % =============================================================
            jobName = 'Contrast_and_ResultsExport';
            jobFile = fullfile(dirs.jobs,[jobName,'_',covname,'_',datestr(now,'yyyymmddTHHMMSS'),'.m']);

            % Initialise job
            jobId = cfg_util('initjob', jobs);

            % If successful save and run job
            sts = cfg_util('isjob_id', jobId);
            if sts
                cfg_util('savejob', jobId, jobFile);
                cfg_util('run', jobId);
            else
                error('Error in initialising job for contrast specification and result export for covariate %s.',covname)
            end
        end
        end
    end
end
