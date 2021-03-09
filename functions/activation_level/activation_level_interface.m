function [subject,properties] = activation_level_interface(subject,properties)

band                                                        = properties.sensor_level_out.band;
%%
%% Defining path
%%
disp('=================================================================');
if(isfield(band,'f_bin'))    
    disp(strcat( 'BC-V-->> Activation level for frequency band: (' , band.name , ') bin ->>>' , string(band.f_bin), 'Hz') );    
    properties.str_band                                     = strcat( band.name,'_',string(band.f_bin),'Hz');
else
    disp(strcat( 'BC-V-->> Activation level for frequency band: (' , band.name , ') ' , string(band.f_start), 'Hz-->' , string(band.f_end) , 'Hz') );
    properties.str_band                                     = strcat( band.name,'_',string(band.f_start),'Hz_',string(band.f_end),'Hz');
end
text_level                                                  = 'Activation_level'; 

%%
%% Band Analysis, activation level
%%
for m=1:length(properties.activation_params.methods)
    analysis_method                                         = properties.activation_params.methods{m};
    fields                                                  = fieldnames(analysis_method);
    method_name                                             = fields{1};
    if(analysis_method.(method_name).run)        
        disp('-----------------------------------------------------------------');
        disp(strcat("-->> Start time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
        switch method_name
            case 'sssblpp'                
                if(properties.BC_V_info.properties.general_params.run_by_trial.value)
                    trial_name                              = properties.trial_name;
                    properties.pathname                     = fullfile(subject.subject_path,trial_name,text_level,'sSSBLpp',band.name);
                else
                    properties.pathname                     = fullfile(subject.subject_path,text_level,'sSSBLpp',band.name);
                end
                if(~isfolder(properties.pathname))
                    mkdir(properties.pathname);
                end
                properties.activation_params.sssblpp_th     = analysis_method.(method_name).sssblpp_th;
                [stat,J,T,indms,properties]                 = activation_level_sssblpp(subject,properties);
                properties.activation_params                = rmfield(properties.activation_params,'sssblpp_th');
            case 'eloreta'
                if(properties.general_params.run_by_trial.value)
                    trial_name                              = properties.trial_name;
                    properties.pathname                     = fullfile(subject.subject_path,trial_name,text_level,'eLORETA',band.name);
                else
                    properties.pathname                     = fullfile(subject.subject_path,text_level,'eLORETA',band.name);
                end
                if(~isfolder(properties.pathname))
                    mkdir(properties.pathname);
                end
                properties.activation_params.gamma1         = analysis_method.(method_name).gamma1;
                properties.activation_params.gamma2         = analysis_method.(method_name).gamma2;
                properties.activation_params.delta_gamma    = analysis_method.(method_name).delta_gamma;
                properties.activation_params.eloreta_th     = analysis_method.(method_name).eloreta_th;
                [stat,J,T,indms,properties]                 = activation_level_eloreta(subject,properties);
                properties.activation_params                = rmfield(properties.activation_params,'gamma1');
                properties.activation_params                = rmfield(properties.activation_params,'gamma2');
                properties.activation_params                = rmfield(properties.activation_params,'delta_gamma');
                properties.activation_params                = rmfield(properties.activation_params,'eloreta_th');
            case 'lcmv'                
                if(properties.general_params.run_by_trial.value)
                    trial_name                              = properties.trial_name;
                    properties.pathname                     = fullfile(subject.subject_path,trial_name,text_level,'LCMV',band.name);
                else
                    properties.pathname                     = fullfile(subject.subject_path,text_level,'LCMV',band.name);
                end
                if(~isfolder(properties.pathname))
                    mkdir(properties.pathname);
                end
                properties.activation_params.gamma1         = analysis_method.(method_name).gamma1;
                properties.activation_params.gamma2         = analysis_method.(method_name).gamma2;
                properties.activation_params.delta_gamma    = analysis_method.(method_name).delta_gamma;
                properties.activation_params.lcmv_th        = analysis_method.(method_name).lcmv_th;
                [stat,J,T,indms,properties]                 = activation_level_lcmv(subject,properties);
                properties.activation_params                = rmfield(properties.activation_params,'gamma1');
                properties.activation_params                = rmfield(properties.activation_params,'gamma2');
                properties.activation_params                = rmfield(properties.activation_params,'delta_gamma');
                properties.activation_params                = rmfield(properties.activation_params,'lcmv_th');
        end
        
        disp(strcat("-->> End time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
        
        reference_path                                                                                          = strsplit(properties.pathname,subject.name);
        if(properties.BC_V_info.properties.general_params.run_by_trial.value)
            if(properties.BC_V_info.properties.general_params.run_frequency_bin.value)
                f_bin                                                                                           = replace(num2str(band.f_bin),'.','_');
                f_bin                                                                                           = strcat(band.name,'_',f_bin);
                properties.BC_V_info.(trial_name).activation_level.(method_name).(band.name).(f_bin).name       = properties.file_name;
                properties.BC_V_info.(trial_name).activation_level.(method_name).(band.name).(f_bin).ref_path   = reference_path{2};
            else
                properties.BC_V_info.(trial_name).activation_level.(method_name).(band.name).name               = properties.file_name;
                properties.BC_V_info.(trial_name).activation_level.(method_name).(band.name).ref_path           = reference_path{2};
            end
        else
            if(properties.BC_V_info.properties.general_params.run_frequency_bin.value)
                f_bin                                                                                           = replace(num2str(band.f_bin),'.','_');
                f_bin                                                                                           = strcat(band.name,'_',f_bin);
                properties.BC_V_info.activation_level.(method_name).(band.name).(f_bin).name                    = properties.file_name;
                properties.BC_V_info.activation_level.(method_name).(band.name).(f_bin).ref_path                = reference_path{2};
            else
                properties.BC_V_info.activation_level.(method_name).(band.name).name                            = properties.file_name;
                properties.BC_V_info.activation_level.(method_name).(band.name).ref_path                        = reference_path{2};
            end
        end
    end
end

end

                                                                                                                                                                           lysis
try
    correlationMats = ROInets.run_individual_network_analysis(dataFile, Settings, resultsName);
catch ME
    fprintf('%s: Test 2 failed. \n', mfilename);
    rethrow(ME);
end%try

%% TEST 3
% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.timeRange                = {[0 120]};                             % range of times to use for analysis
Settings.Regularize.do            = false;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.leakageCorrectionMethod  = 'pairwise';                      % choose from 'closest', 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
Settings.EnvelopeParams.windowLength = 2; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.overlap   = 0.5;
Settings.EnvelopeParams.useFilter    = false;                        % use a more sophisticated filter than a sliding window average
Settings.EnvelopeParams.takeLogs  = false;
Settings.frequencyBands           = {[]};                       % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'mean';                          % 'PCA', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'fixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = sessionName; 
Settings.SaveCorrected.timeCourse    = false;
Settings.SaveCorrected.envelopes      = false;
Settings.SaveCorrected.variances     = false;
Settings.SaveCorrected.ROIweightings = false;

% run the ROI network analysis
try
    correlationMats = ROInets.run_individual_network_analysis(dataFile, Settings, resultsName);
catch ME
    fprintf('%s: Test 3 failed. \n', mfilename);
    rethrow(ME);
end%try



%% Test 4 trial data, one subj
dataFile = '/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession12_spm_meeg.mat';
% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.timeRange                = {[0.01 3.99]};                             % range of times to use for analysis
Settings.Regularize.do            = false;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.leakageCorrectionMethod  = 'closest';                      % choose from 'closest', 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 1;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.EnvelopeParams.windowLength = 1/40; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.useFilter = true;                        % use a more sophisticated filter than a sliding window average
Settings.EnvelopeParams.takeLogs  = true;
Settings.frequencyBands           = {[4 30]};                       % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'PCA';                          % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'mixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = sessionName; 
Settings.SaveCorrected.timeCourse    = false;
Settings.SaveCorrected.envelopes      = false;
Settings.SaveCorrected.variances     = false;
Settings.SaveCorrected.ROIweightings = false;
Settings.SubjectLevel.conditionLabel   = {'longvalidR', 'longvalidL', 'ShortValidRight', 'ShortValidLeft'};
Settings.SubjectLevel.designSummary    = {[1 0 0 0]', [0 1 0 0], [0 0 1 0]', [0 0 0 1]};
Settings.SubjectLevel.contrasts        = {[1 0 1 0]; [1 -1 1 -1]};

Settings = ROInets.check_inputs(Settings);
try
    CorrMats = ROInets.run_individual_network_analysis_task(dataFile, ...
                                                           Settings, ...
                                                            resultsName);
catch ME
    fprintf('%s: Test 4 failed. \n', mfilename);
    rethrow(ME);
end%try

%% Test 4 trial data, several subj
dataFiles = {'/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession12_spm_meeg.mat'; ...
             '/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession120_spm_meeg.mat'; ...
             '/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession121_spm_meeg.mat'; ...
             '/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession122_spm_meeg.mat'};
% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.timeRange                = [0.01 3.99];                             % range of times to use for analysis
Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.Regularize.path          = 0.001;              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
Settings.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
Settings.Regularize.adaptivePath  = false;                           % adapth the regularization path if necessary
Settings.leakageCorrectionMethod  = 'closest';                      % choose from 'closest', 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 1;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.EnvelopeParams.windowLength = 1/40; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.useFilter = true;                        % use a more sophisticated filter than a sliding window average
Settings.EnvelopeParams.takeLogs  = true;
Settings.frequencyBands           = {[]};                       % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'PCA';                          % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'mixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = {'sess1', 'sess2', 'sess3', 'sess4'}; 
Settings.SaveCorrected.timeCourse    = false;
Settings.SaveCorrected.envelopes      = false;
Settings.SaveCorrected.variances     = false;
Settings.SaveCorrected.ROIweightings = false;
Settings.SubjectLevel.conditionLabel   = {'longvalidR', 'longvalidL', 'ShortValidRight', 'ShortValidLeft'};
Settings.SubjectLevel.designSummary    = {[1 0 0 0]', [0 1 0 0], [0 0 1 0]', [0 0 0 1]};
Settings.SubjectLevel.contrasts        = {[1 0 1 0]; [1 -1 1 -1]};
Settings.GroupLevel.designMatrix       = [1 0
                                          1 0
                                          0 1
                                          0 1];
Settings.GroupLevel.contrasts          = [1  1;  % contrast 1
                                          1 -1]; % contrast 2

try
    CorrMats = run_network_analysis(dataFiles, Settings);
catch ME
    fprintf('%s: Test 4 failed. \n', mfilename);
    rethrow(ME);
end%try

%% We've got to the end!
FAIL = 0;
ROInets.call_fsl_wrapper(['rm -rf ' dataDir]);

end%test_me
% [EOF]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          on = Settings.nSessions:-1:1,
    fprintf('\n\n%s: Individual correlation analysis for file %d out of %d\n', ...
            mfilename, Settings.nSessions - iSession + 1, Settings.nSessions);

    D                          = Dlist{iSession};
    sessionName                = Settings.sessionName{iSession};
    matsSaveFileName{iSession} = fullfile(outputDirectory,                                      ...
                                          sprintf('%s_single_session_correlation_mats_tmp.mat', ...
                                                  sessionName));

    if strcmpi(Settings.paradigm, 'task'),
        mats{iSession} = ROInets.run_individual_network_analysis_task(D,                          ...
                                                                 Settings,                   ...
                                                                 matsSaveFileName{iSession}, ...
                                                                 iSession);
    elseif strcmpi(Settings.paradigm, 'rest'),
        mats{iSession} = ROInets.run_individual_network_analysis(D,                          ...
                                                                 Settings,                   ...
                                                                 matsSaveFileName{iSession}, ...
                                                                 iSession);
    else
        error([mfilename ':BadParadigm'], ...
              'Unrecognised paradigm %s. \n', Settings.paradigm);
    end%if
end%for

% reformat results - correlationMats is a cell array of frequency bands
correlationMats = ROInets.reformat_results(mats, Settings);

% save current results: will write over later
% just in case of crash at group level
saveFileName = fullfile(outputDirectory, 'ROInetworks_correlation_mats.mat');
save(saveFileName, 'correlationMats');
clear mats

%% Subject-level analysis to average over sessions in a fixed-effects manner
correlationMats = ROInets.do_subject_level_glm(correlationMats, Settings);

%% Group-level analysis
% Find whole group means
if strcmpi(Settings.paradigm, 'rest'),
    correlationMats = ROInets.do_group_level_statistics(correlationMats, Settings);
end%if

% Perform group-level GLM
if ~isempty(Settings.GroupLevel),
    correlationMats = ROInets.do_group_level_glm(correlationMats, Settings);
end%if

%% save matrices
fprintf('\n%s: Saving Results. \n', mfilename);

% save collected results
save(saveFileName, 'correlationMats');

% we stored individual results as we went along, in case of crash. Delete
% them if we've safely got to this stage. 
for iSession = 1:length(matsSaveFileName),
    delete(matsSaveFileName{iSession});
end%for
 
% tidy output of funciton
Settings.correlationMatsFile = saveFileName;
save(fullfile(outputDirectory, 'ROInetworks_settings.mat'), 'Settings');

fprintf('%s: Analysis complete. \n\n\n', mfilename);
end%run_network_analysis
% [EOF]
                                                                                                                  time,            ...
                                                           Settings.EnvelopeParams); %#ok<ASGLU>
    end%for loop over trials
    
    % save power envelopes
    if Settings.SaveCorrected.envelopes,
        saveDir = fullfile(Settings.outputDirectory, ...
                           'corrected-ROI-timecourses', filesep);
        ROInets.make_directory(saveDir);
        saveFile = fullfile(saveDir,                                      ...
                            sprintf('%s_%s_ROI_envelope_timecourses.mat', ...
                                    sessionName, bandName));
        save(saveFile, 'nodeEnv', 'time_ds');
    end%if
        
    %% Run correlation analysis 
    % calculate correlation matrices. 
    CorrMats{iFreq} = ROInets.run_correlation_analysis([],           ...
                                                       nodeEnv,      ...
                                                       Settings.Regularize);
    
    
    % Use an empirical null to enable conversion to z-stats
    transformSurrogates = ~Settings.EnvelopeParams.takeLogs;
    RegParams           = struct('do', Settings.Regularize.do, ...
                                 'rho', CorrMats{iFreq}.Regularization.mean);
    sigma = ROInets.find_permutation_H0_distribution_width(nodeEnv,                    ...
                                                           Settings.nEmpiricalSamples, ...
                                                           RegParams,                  ...
                                                           transformSurrogates);
          
    CorrMats{iFreq}.H0Sigma = sigma;
    
    % Store session name
    CorrMats{iFreq}.sessionName = sessionName;
    CorrMats{iFreq}.timeWindow  = timeRange;
    
    % free up some memory
    clear nodeData nodeEnv
    
    %% conversion of correlations to z-stats
    fprintf(' Converting correlations to normal z-stats\n');
    CorrMats{iFreq} = ROInets.convert_correlations_to_normal_variables(CorrMats{iFreq}, ...
                                                                       sigma,      ...
                                                                       doRegularize);
    
    %% Run first-level GLM
    CorrMats{iFreq}.firstLevel = run_first_level_glm(CorrMats{iFreq},    ...
                                                     designMat,          ...
                                                     Settings.SubjectLevel.contrasts, ...
													 D.fname);
    % clean up filtered object
    delete(cleanD);
end%loop over freq bands

%% save results to disc to provide backup of completed session analyses
save(resultsSaveName, 'CorrMats');




%%% END OF FUNCTION PROPER %%%

end%run_individual_correlation_analysis
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
function FirstLevel = run_first_level_glm(CorrMats, designMat, contrasts, fileName)
%RUN_FIRST_LEVEL_GLM
%


% input checking
[nTrials, nRegressors] = size(designMat);
nContrasts             = length(contrasts);
[~, nModes, checkMe]   = size(CorrMats.envCorrelation_z);
assert(checkMe == nTrials,         ...
      [mfilename ':LostTrials'],   ...
      'Number of trials must match columns of design matrix. \n');

assert(iscell(contrasts),               ...
       [mfilename ':NonCellContrasts'], ...
       'Contrasts must be a cell array. \n');
assert(all(cellfun(@length, contrasts) == nRegressors), ...
       [mfilename ':BadContrastFormat'],                ...
       'All contrasts must have the same length as number of regressors. \n');
   
% make sure contrasts are formatted as a cell array of column vectors
useContrasts = cell(1,nContrasts);
for iContrast = 1:nContrasts,
    useContrasts{iContrast} = contrasts{iContrast}(:);
end%for
   
% Precompute some helpful things
XtX       = designMat' * designMat;
[RXtX, p] = chol(XtX);
if ~p,
	invXtX = RXtX \ (RXtX' \ eye(nRegressors));
	pinvX  = RXtX \ (RXtX' \ designMat');
	hasBadEVs    = false;
	badContrasts = false(nContrasts, 1);
else
	% design matrix was rank deficient
	% is that because we have missing information for certain trial types?
	badEVs    = all(0 == designMat);
	hasBadEVs = any(badEVs);
	if hasBadEVs,
		warning([mfilename ':MissingTrials'],                   ...
			    '%s: file %s is missing trials for %d EVs. \n', ...
				mfilename, fileName, sum(badEVs));
		badContrasts = logical(cellfun(@(C) any(C(badEVs)), useContrasts));
		invXtX = pinv(XtX);
		pinvX  = invXtX * designMat';
	else
		error([mfilename ':RankDeficientDesign'],                     ...
			  ['%s: the design matrix is rank deficient. ',           ...
			   'Check that you''ve specified your EVs sensibly. \n'], ...
			  mfilename);
	end%if
end%if

% declare memory
[rho, prho, prhoReg] = deal(zeros(nModes, nModes, nContrasts));

% run GLM on each edge
for i = 1:nModes,
    for j = i+1:nModes,
        rho(i,j,:) = glm_fast_for_meg(squeeze(CorrMats.envCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
        prho(i,j,:) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
								  
	    % fill in uninformative values with NaN.
		if hasBadEVs,
			rho(i,j,badContrasts)  = NaN;
			prho(i,j,badContrasts) = NaN;
		end%if
        if isfield(CorrMats, 'envPartialCorrelationRegularized_z'),
        prhoReg(i,j,1:nContrasts) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelationRegularized_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0); 
        prhoReg(i,j,badContrasts) = NaN;
        else
            prhoReg(i,j) = 0;
        end%if
    end%for
end%for

% symmetrise and reformat
for iContrast = nContrasts:-1:1,
    FirstLevel(iContrast).cope.correlation                   = rho(:,:,iContrast) + rho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelation            = prho(:,:,iContrast) + prho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelationRegularized = prhoReg(:,:,iContrast) + prhoReg(:,:,iContrast)';
end%for

end%run_first_level_glm






%--------------------------------------------------------------------------
function [designMat, goodTrials, trialID, nConditions] = set_up_first_level(D, Settings, excludeTrials)
%SET_UP_FIRST_LEVEL creates the design matrix and trial identifiers

nConditions = length(Settings.SubjectLevel.conditionLabel);

% hold the relevant indices for trials in each condition
for iCondition = nConditions:-1:1,
    tI = D.indtrial(Settings.SubjectLevel.conditionLabel{iCondition}, ...
                    'GOOD');
    trialInds{iCondition} = ROInets.setdiff_pos_int(tI, excludeTrials);
    % check we've found something                               
    if isempty(trialInds{iCondition}), 
        warning([mfilename ':EmptyCondition'], ...
                'No good trials found in %s for condition %s. \n', ...
                D.fname, Settings.SubjectLevel.conditionLabel{iCondition});
    end%if
end%for

% extract a list of all good, relevant trials
goodTrials = sort([trialInds{:}]);

% generate a set of IDs linking trials to condition number
trialID = zeros(length(goodTrials), 1);
for iCondition = nConditions:-1:1,
    trialID(ismember(goodTrials,trialInds{iCondition})) = iCondition;
end%for

% check design matrix size
assert(all(cellfun(@length, Settings.SubjectLevel.designSummary) == nConditions), ...
       [mfilename ':DesignMatrixSizeFault'],                                      ...
       ['The design matrix summary must, in each cell, contain a vector of the ', ...
        'same length as the number of conditions. \n']);
% use an OSL function to generate the subject-specific design matrix
designMat = oat_setup_designmatrix(struct('Xsummary', {Settings.SubjectLevel.designSummary}, ...
                                          'trialtypes', trialID));
end%set_up_first_level
    
%--------------------------------------------------------------------------
function [t, tI, tR] = time_range(time, timeRange, iSession)
%TIME_RANGE selects time range for analysis of each session
% TIME is a vector of times
% TIMERANGE is either a cell array of two-component vectors, a single
% two-component vector, or a null vector

if isempty(timeRange),
    % use the whole time range
    t = time;
    tI = true(size(time));
    tR = [];
else
    % subselect time range
    if iscell(timeRange),
        tR = timeRange{iSession};
    else
        tR = timeRange;
    end%if
    validateattributes(tR, {'numeric'}, ...
                       {'vector', 'numel', 2, 'nondecreasing'}, ... % can have negative times in task data. 
                       'time_range', 'timeRange', 2);
                   
    tI = (time <= tR(2)) & (time >= tR(1));
    t  = time(tI);
end%if
end%time_range



%--------------------------------------------------------------------------
% [EOF]
                                                                                                                                                                                                                                                                                                                                                                                                        (std(voxelData, [], 2), eps);
        voxelDataScaled = ROInets.demean(voxelData, 2);
        clear voxelData
        
        % pre-allocate PCA weightings for each parcel
        voxelWeightings = zeros(size(spatialBasis));
        
        % perform PCA on each parcel and select 1st PC scores to represent
        % parcel
        for iParcel = nParcels:-1:1,
%             progress = nParcels - iParcel + 1;
%             ft_progress(progress / nParcels, ...
%                         [mfilename ...
%                          ':    Finding PCA time course for ROI %d out of %d'], ...
%                         iParcel, nParcels);
                
            thisMask = spatialBasis(:, iParcel);
            if any(thisMask), % non-zero
                parcelData = voxelDataScaled(thisMask, :);
                
                [U, S, V]  = ROInets.fast_svds(parcelData, 1);
                PCAscores  = S * V';
                
                % restore sign and scaling of parcel time-series
                % U indicates the weight with which each voxel in the
                % parcel contributes to the 1st PC
                TSsign          = sign(mean(U));
                relVoxelWeights = abs(U) ./ sum(abs(U)); % normalise the linear combination
                % weight the temporal STDs from the ROI by the proportion used in 1st PC
                TSscale         = dot(relVoxelWeights, temporalSTD(thisMask)); 
                nodeTS          = TSsign .*                               ...
                                  (TSscale / max(std(PCAscores), eps)) .* ... 
                                  PCAscores;
                      
                % return the linear operator which is applied to the data
                % to retrieve the nodeTS
                voxelWeightings(thisMask, iParcel) = TSsign .* ...
                                                     (TSscale / max(std(PCAscores), eps)) ...
                                                     .* U';
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeTS = zeros(1, ROInets.cols(voxelDataScaled));
            end%if
            
            nodeData(iParcel,:) = nodeTS;
        end%for
        
        clear parcelData voxelDataScaled

    case 'peakvoxel'
        if any(spatialBasis(:)~=0 & spatialBasis(:)~=1),
            warning([mfilename ':NonBinaryParcelMask'],    ...
                    ['Input parcellation is not binary. ', ...
                     'It will be binarised. \n']);
        end%if
        spatialBasis = logical(spatialBasis);
        
        % find rms power in each voxel
        voxelPower = sqrt(ROInets.row_sum(voxelData.^2) ./ ...
                          ROInets.cols(voxelData));
                      
        % pre-allocate weightings for each parcel
        voxelWeightings = zeros(size(spatialBasis));
                      
        % take peak voxel in each parcel
        for iParcel = nParcels:-1:1,
%             progress = nParcels - iParcel + 1;
%             ft_progress(progress / nParcels, ...
%                         [mfilename ...
%                          ':    Finding peak voxel time course for ROI %d out of %d'], ...
%                         iParcel, nParcels);
            
            thisMask = spatialBasis(:, iParcel);
            
            if any(thisMask), % non-zero
                % find index of voxel with max power
                thisParcPower            = voxelPower;
                thisParcPower(~thisMask) = 0;
                [~, maxPowerInd]         = max(thisParcPower);
                
                % select voxel timecourse
                nodeData(iParcel,:) = voxelData(maxPowerInd,:);
                
                % save which voxel was used
                voxelWeightings(maxPowerInd, iParcel) = 1;
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeData(iParcel,:) = zeros(1, ROInets.cols(voxelData));
            end%if
        end%loop over parcels
        
        clear voxelData parcelData
        
    case 'spatialbasis'
        % scale group maps so all have a positive peak of height 1
        % in case there is a very noisy outlier, choose the sign from the
        % top 5% of magnitudes
        top5pcInd = abs(spatialBasis) >=                        ...
                         repmat(prctile(abs(spatialBasis), 95), ...
                                [ROInets.rows(spatialBasis), 1]);
        for iParcel = nParcels:-1:1,
            mapSign(iParcel) = sign(mean(...
                              spatialBasis(top5pcInd(:,iParcel), iParcel)));
        end%for
        scaledSpatialMaps = ROInets.scale_cols(spatialBasis, ...
                                   mapSign ./                ...
                                   max(max(abs(spatialBasis), [], 1), eps));
        
        % find time-course for each spatial basis map
        for iParcel = nParcels:-1:1, % allocate memory on the fly
%             progress = nParcels - iParcel + 1;
%             ft_progress(progress / nParcels, ...
%                         [' ' mfilename ...
%                          ':    Finding spatial basis time course for ROI %d out of %d'], ...
%                         iParcel, nParcels);
            
            % extract the spatial map of interest
            thisMap     = scaledSpatialMaps(:, iParcel);
            parcelMask  = logical(thisMap);
            
            % estimate temporal-STD for normalisation
            temporalSTD = max(std(voxelData, [], 2), eps);
            
            % variance-normalise all voxels to remove influence of
            % outliers. - remove this step 20 May 2014 for MEG as data are
            % smooth and little risk of high-power outliers. Also, power is
            % a good indicator of sensible signal. 
            % Weight all voxels by the spatial map in question
            weightedTS  = ROInets.scale_rows(voxelData, thisMap);
            
            % perform svd and take scores of 1st PC as the node time-series
            % U is nVoxels by nComponents - the basis transformation
            % S*V holds nComponents by time sets of PCA scores - the 
            % timeseries data in the new basis
            [U, S, V]   = ROInets.fast_svds(weightedTS(parcelMask,:), 1);
            clear weightedTS
            
            PCAscores   = S * V';
            maskThresh  = 0.5; % 0.5 is a decent arbitrary threshold chosen by Steve Smith and MJ after playing with various maps.
            thisMask    = thisMap(parcelMask) > maskThresh;   
            
            if any(thisMask), % the mask is non-zero
                % U is the basis by which voxels in the mask are weighted
                % to form the scores of the 1st PC
                relativeWeighting = abs(U(thisMask)) ./ ...
                                    sum(abs(U(thisMask)));
                
                TSsign  = sign(mean(U(thisMask)));
                TSscale = dot(relativeWeighting, temporalSTD(thisMask));       
                nodeTS  = TSsign .*                               ...
                          (TSscale / max(std(PCAscores), eps)) .* ...      
                          PCAscores;
                      
                % for Mark: this is the linear operator which is applied to
                % the voxel data to get nodeTS.
                voxelWeightings(parcelMask,iParcel) = TSsign .* ...
                                             (TSscale / max(std(PCAscores), eps)) ...
                                             .* (U' .* thisMap(parcelMask)');
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeTS = zeros(1, ROInets.cols(weightedTS));
                voxelWeightings(~thisMask, iParcel) = zeros(length(thisMask), 1);
            end%if
            
            nodeData(iParcel, :) = nodeTS;
            
        end%loop over parcels
        
        clear voxelData 
        
        
    otherwise
        error([mfilename ':UnrecognisedTimeCourseMethod'],            ...
              ['Unrecognised method for finding ROI time-course. \n', ...
               'Expected ''PCA'', ''spatialBasis'', or ''mean''. \n']);
end%switch

% ft_progress('close');

end%get_node_tcs
%--------------------------------------------------------------------------
end%run_individual_correlation_analysis
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
function balanced = balance_correlations(x)
%BALANCE_CORRELATIONS makes matrices symmetric

if ismatrix(x) && ~isvector(x),
    balanced = (x + x') ./ 2.0;
else
    balanced = x;
end%if
end%balance_correlations






%--------------------------------------------------------------------------
function [t, tI, tR] = time_range(time, timeRange, iSession)
%TIME_RANGE selects time range for analysis of each session
% TIME is a vector of times
% TIMERANGE is either a cell array of two-component vectors, a single
% two-component vector, or a null vector

if isempty(timeRange),
    % use the whole time range
    t = time;
    tI = true(size(time));
    tR = [];
else
    % subselect time range
    if iscell(timeRange),
        tR = timeRange{iSession};
    else
        tR = timeRange;
    end%if
    validateattributes(tR, {'numeric'}, ...
                       {'vector', 'numel', 2, 'nonnegative', 'nondecreasing'}, ...
                       'time_range', 'timeRange', 2);
                   
    tI = (time <= tR(2)) & (time >= tR(1));
    t  = time(tI);
end%if
end%time_range





%--------------------------------------------------------------------------
function [] = save_corrected_timecourse_results(nodeData, ...
    allROImask, voxelWeightings, Settings, sessionName, protocol, bandName)
% Save the weightings over voxels used to calculate the time-course
% for each ROI
if Settings.SaveCorrected.ROIweightings,
    allVoxelWeightings                = zeros(ROInets.rows(allROImask), ...
                                              ROInets.cols(voxelWeightings));
    allVoxelWeightings(allROImask, :) = voxelWeightings;

    saveDir             = fullfile(Settings.outputDirectory, ...
                                   'spatialBasis-ROI-weightings', filesep);
    ROItcWeightSaveFile = fullfile(saveDir,                      ...
                                   [sessionName '_' protocol '_' ...
                                    bandName '_ROI_timecourse_weightings']);
    ROInets.make_directory(saveDir);
    try
        nii.quicksave(allVoxelWeightings, ROItcWeightSaveFile, Settings.gridStep);
    catch % perhaps we have a weird number of voxels
        save(ROItcWeightSaveFile, 'allVoxelWeightings');
    end%try
end%if

% save node data
if Settings.SaveCorrected.timeCourses,
    saveDir = fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses', filesep);
    ROInets.make_directory(saveDir);
    saveFile = fullfile(saveDir, [sessionName '_correction-' protocol '_' bandName '_ROI_timecourses.mat']);
    save(saveFile, 'nodeData');
end%if

% save variance in each ROI
if Settings.SaveCorrected.variances,
    saveDir = fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses', filesep);
    ROInets.make_directory(saveDir);
    varSaveFile = fullfile(saveDir, [sessionName '_correction-' protocol '_' bandName '_ROI_variances.mat']);
    ROIvariances = var(nodeData, [], 2);                                   %#ok<NASGU>
    save(varSaveFile, 'ROIvariances');
end%if
end%save_corrected_timecourse_results
%--------------------------------------------------------------------------
% [EOF]
                                                                    martshare_history: historycount = %d, peercount = %d
 smartmem: host->memavail = %llu
 smartmem: MemTotal       = %llu (%f GB)
 smartmem: MemFree        = %llu (%f GB)
 smartmem: Buffers        = %llu (%f GB)
 smartmem: Cached         = %llu (%f GB)
 smartmem: NumPeers       = %u
 smartmem: MemReserved    = %llu (%f GB)
 smartmem: MemSuggested   = %llu (%f GB)
 smartcpu_update: switching to zombie
 smartcpu_update: ProcessorCount = %d
 smartcpu_update: NumPeers       = %d
 smartcpu_update: BogoMips       = %.2f
 smartcpu_update: AvgLoad        = %.2f
 smartcpu_update: CpuLoad        = %.2f %%
 smartcpu_update: host->status   = %u
 smartcpu_update: switching to idle
 open_uds_connection socket error: open_uds_connection socket
 open_uds_connection connect error: open_uds_connection connect
 open_uds_connection: connected to %s on socket %d
 open_tcp_connection: using direct memory copy
 open_tcp_connection: server = %s, port = %u
 open_tcp_connection: nslookup1 failed on '%s'
 open_tcp_connection: nslookup2 failed on '%s'
 open_tcp_connection: socket = %d
 open_tcp_connection error: open_tcp_connection
 open_tcp_connection: connectioncount = %d
 open_tcp_connection: connected to %s:%u on socket %d
 close_connection: socket = %d
 close_connection error: close_connection
 close_connection: connectioncount = %d
                          †  4   4   _ö      4                                   zR xê  $      Zˇˇˇˇˇˇ=        AÜC       $   D   0Zˇˇˇˇˇˇ∏        AÜC       $   l   »ZˇˇˇˇˇˇÓ       AÜC       $   î   ê\ˇˇˇˇˇˇÇ        AÜC       $   º   ¯\ˇˇˇˇˇˇ2	       AÜC              zR xê  $      ¯eˇˇˇˇˇˇi        AÜC       $   D   @fˇˇˇˇˇˇK       AÜC       $   l   hgˇˇˇˇˇˇÎ       AÜC       $   î   0jˇˇˇˇˇˇñ       AÜC              zR xê  $      êmˇˇˇˇˇˇ       AÜC       $   D   àoˇˇˇˇˇˇ(
       AÜC              zR xê  $      xyˇˇˇˇˇˇ÷        AÜC       $   D   0zˇˇˇˇˇˇR       AÜC       $   l   hÄˇˇˇˇˇˇ¨       AÜC              zR xê  $      ÿÇˇˇˇˇˇˇV        AÜC       $   D   ÉˇˇˇˇˇˇU       AÜC       $   l   Hâˇˇˇˇˇˇ2       AÜC              zR xê  $      Hãˇˇˇˇˇˇm        AÜC       $   D   êãˇˇˇˇˇˇ:       AÜC       $   l   ®åˇˇˇˇˇˇg       AÜC       $   î   çˇˇˇˇˇˇ"       AÜC       $   º   ¯èˇˇˇˇˇˇp        AÜC       $   ‰   @êˇˇˇˇˇˇp        AÜC       $     àêˇˇˇˇˇˇ_        AÜC       $   4  ¿êˇˇˇˇˇˇ¬        AÜC       $   \  hëˇˇˇˇˇˇF       AÜC       $   Ñ  êíˇˇˇˇˇˇó        AÜC       $   ¨  ìˇˇˇˇˇˇø        AÜC       $   ‘  †ìˇˇˇˇˇˇø        AÜC       $   ¸  8îˇˇˇˇˇˇø        AÜC       $   $  –îˇˇˇˇˇˇø        AÜC       $   L  hïˇˇˇˇˇˇø        AÜC       $   t   ñˇˇˇˇˇˇø        AÜC       $   ú  òñˇˇˇˇˇˇ        AÜC       $   ƒ  ÄñˇˇˇˇˇˇÄ        AÜC              zR xê  $      ¿ñˇˇˇˇˇˇ>       AÜC       $   D   ÿóˇˇˇˇˇˇv       AÜC              zR xê  $      ûˇˇˇˇˇˇ"       AÜC       $   D    üˇˇˇˇˇˇø       AÜC       $   l   ∏•ˇˇˇˇˇˇ#        AÜC              zR xê  $      ®•ˇˇˇˇˇˇ       AÜC       $   D   †¶ˇˇˇˇˇˇœ       AÜC              zR xê  $      0∂ˇˇˇˇˇˇ         AÜC       $   D   ÿ∂ˇˇˇˇˇˇP       AÜC       $   l    ∏ˇˇˇˇˇˇP       AÜC       $   î   (πˇˇˇˇˇˇP       AÜC              zR xê  $      8∫ˇˇˇˇˇˇÖ       AÜC              zR xê  $      àºˇˇˇˇˇˇb        AÜC       $   D   –ºˇˇˇˇˇˇJ       AÜC       $   l   ¯¿ˇˇˇˇˇˇR       AÜC              zR xê  $      √ˇˇˇˇˇˇì        AÜC       $   D   ê√ˇˇˇˇˇˇ∆
       AÜC              zR xê  $      0Œˇˇˇˇˇˇû        AÜC       $   D   ®Œˇˇˇˇˇˇ       AÜC              zR xê  $      à‘ˇˇˇˇˇˇÁ       AÜC       $   D   P÷ˇˇˇˇˇˇ)       AÜC       $   l   XŸˇˇˇˇˇˇ#        AÜC       $   î   `Ÿˇˇˇˇˇˇ       AÜC        ¿_ˇ  ¿_ˇ          ¯õ      ú      ú      ú      (ú      4ú      @ú      Lú      Xú      dú      pú      |ú      àú      îú      †ú      ¨ú      ∏ú      ƒú      –ú      ‹ú      Ëú      Ùú       ù      ù      ù      $ù      0ù      <ù      Hù      Tù      `ù      lù      xù      Ñù      êù      úù      ®ù      ¥ù      ¿ù      Ãù      ÿù      ‰ù      ù      ¸ù      û      û       û      ,û      8û      Dû      Pû      \û      hû      tû      Äû      åû      òû      §û      ∞û      ºû      »û      ‘û      ‡û      Ïû      ¯û      ü      ü      ü              ®ü      ≠ü              ß´™2                                                            ß´™2                                                                   ª±∞<                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                            ß´™2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (     0     8     @     H     P     X     `     h     p     x     Ä     à     ê     ò     †     ®     ∞     ∏     ¿     »     –     ÿ     ‡     Ë          ¯                           (    0    8    @    H    P    X    `    h    p    x    Ä    à    ê    ò    †    ®    ∞    ∏    ¿    »    –    ÿ    ‡    Ë        ¯                          (    0    H    @    †@¿ê¿p–††∞‡‡∞`‡¿p¿∞pp`––†¿¿¿¿¿¿Ä¿Ä∞¿0†–––––êp–‡†‡††∞0        d           +   d           8   f ¡∑U       .  ¿      x   $  ¿         $   @          N  @          .         å   $            $   ¿          N  ¿          .  ¿      ü   $  ¿         $            N           .  ∞      ´   $  ∞         $   ê          N  ê          .  @      ¥   $  @         $   2	         N  2	      ¡               œ               ‡               È               Ú                                d             d             d           !  f Ä∑U       .  Ä      T  $  Ä         $   p          N  p          .        [  $           $   P         N  P         .  @      m  $  @         $            N           .  0!      w  $  0!         $   ñ         N  ñ         d             d           Ü  d           ë  f Ä∑U       .  –$      ƒ  $  –$         $             N            .  &      ÷  $  &         $   (
         N  (
         d             d           ‡  d           È  f Ä∑U       .   1        $   1         $   ‡          N  ‡          .   2      *  $   2         $   `         N  `         .  `8      2  $  `8         $   ¨         N  ¨         d             d           B  d           K  f Å∑U    |              ä              ñ              £              µ              «              Ÿ              Ô              ˙              	                            +              @              U              k                            î              £              ≤              ¡              œ              ‡              Ò                                          "              0              =              J              W              h              n              x              Å              ê              †              ∞              ¡              –              ‡                            ˙                                                           d             d           #  d           .  f Å∑U       .  ;      a  $  ;         $   `          N  `          .  p;      g  $  p;         $   `         N  `         .  –A      q  $  –A         $   2         N  2         d             d           {  d           Ç  f Å∑U       .  D      ±  $  D         $   p          N  p          .  ÄD      æ  $  ÄD         $   @         N  @         .  ¿E      «  $  ¿E         $   p         N  p         .  0G      —  $  0G         $   0         N  0         .  `I      Ÿ  $  `I         $   p          N  p          .  –I      „  $  –I         $   p          N  p          .  @J      Ó  $  @J         $   `          N  `          .  †J      ˙  $  †J         $   –          N  –          .  pK      
  $  pK         $   P         N  P         .  ¿L        $  ¿L         $   †          N  †          .  `M      /  $  `M         $   ¿          N  ¿          .   N      D  $   N         $   ¿          N  ¿          .  ‡N      Z  $  ‡N         $   ¿          N  ¿          .  †O      o  $  †O         $   ¿          N  ¿          .  `P      Ö  $  `P         $   ¿          N  ¿          .   Q      ú  $   Q         $   ¿          N  ¿          .  ‡Q      ≤  $  ‡Q         $             N            .  Q      √  $  Q         $   Ä          N  Ä          d             d           À  d           ◊  f Å∑U       .  pR        $  pR         $   @         N  @         .  ∞S        $  ∞S         $   v         N  v         d             d           )  d           5  f Å∑U       .  0Z      i  $  0Z         $   0         N  0         .  `[      |  $  `[         $   ¿         N  ¿         .   b      á  $   b      ï  Ñ              $   #          N  #          d             d              d           ,  f Å∑U       .  Pb      `  $  Pb         $             N            .  pc      s  $  pc         $   œ         N  œ         d             d           ~  d           â  f Ç∑U       .  @s      º  $  @s         $   –          N  –          .  t      Ã  $  t         $   P         N  P         .  `u      ﬂ  $  `u         $   P         N  P         .  ∞v      Û  $  ∞v         $   P         N  P         d             d           	  d           	  f Ç∑U       .   x      F	  $   x         $   Ö         N  Ö         d             d           W	  d           d	  f Ç∑U       .  êz      ô	  $  êz         $   p          N  p          .   {      ´	  $   {         $   P         N  P         .  P      Ω	  $  P         $   R         N  R         d             d           —	  d           ‹	  f Ç∑U       .  ∞Å      
  $  ∞Å         $   †          N  †          .  PÇ      
  $  PÇ         $   ÿ
         N  ÿ
         d             d           /
  d           :
  f Ç∑U       .  0ç      m
  $  0ç         $   †          N  †          .  –ç      |
  $  –ç         $            N           d             d           ç
  d           ó
  f Ç∑U       .  ì      …
  $  ì         $            N           .  ‡ï      ﬁ
  $  ‡ï         $   0         N  0         .  ô      Û
  $  ô         $   0          N  0          .  @ô        $  @ô         $            N           d              †      ,    ¥      ?    ¿      S           f    ¿      r    ∞      {    Ä      Ç          î    @      û    0!      ≠    –$      ø    &      …     1      Ÿ     2      ·    `8      Ò    ;      ˜    p;          –A          D          ÄD      !    ¿E      +    0G      3    `I      =    –I      H    @J      T    †J      d    pK      s    ¿L      â    `M      û     N      ¥    ‡N      …    †O      ﬂ    `P      ˆ     Q          ‡Q          Q      %    pR      8    ∞S      C    0Z      V    `[      a     b      o    Pb      Ç    pc      ç    @s      ù    t      ∞    `u      ƒ    ∞v      ◊     x      Ë    êz      ˙     {          P           ∞Å      /    PÇ      @    0ç      O    –ç      `    ì      u    ‡ï      ä    ô      ò    @ô      ™    X¬      ∏    ò¬      …    ÿ¬      ◊    ‡¬      „    √          P√          ê√          –√      &    ƒ      <    Pƒ      G    êƒ      V    –ƒ      d    ≈      x    P≈      ç    ê≈      ¢    –≈      ∏    ∆      Ã    P∆      ·    ê∆          –∆      ˇ    «          P«          ê«      -    –«      6    ÿ«      ?    ‡«      Q    Ë«      c    «      t    Ù«      Ö    ¯«      ï    ¸«      •     »      ≥    »      ¿    »      Õ    »      ⁄    »      Î    »      Ò     »      ˚    (»          0»          8»      #    @»      3    H»      D    P»      S    X»      c    `»      s    h»      }    à»      á    ò»      ë    ®»      ö    »»      ¶    @      ≥            º                         Ÿ            Î            ˛                                    #            )            0            7            @            G            Q            W            ^            d            q            z            â            ñ            ¢            Ø            ∑            ¡            Ã            ◊            ﬂ            Ô            ˜            ˇ            
                        #            4            J            d            n            Ü            ë            ô            °            ±            ∆            ÷            Ê            Ù                                    ,            @            F            O            U            _            e            m            y            Å            à            î            ú            §            ±            ∫            ≈            À            ”            Ì  Ó  Ô    Ò  Í  Î  Ï  À  Ã  Õ  Œ  –  —  “  ”  ‘  ’  ÷  ◊  ÿ  Ÿ  ⁄  €  ‹  ›  ﬁ  ﬂ  ‡  ·  ‚  „  ‰  Â  Ê  Á  Ë  È  Ú  Û  Ù  ı  ˆ  ˜  ¯  ˘  ˙  ˚  ¸  ˝  ˛  ˇ                     	  
            œ  Ì  Ó  Ô    Ò  Í  Î  Ï  À  Ã  Õ  Œ  –  —  “  ”  ‘  ’  ÷  ◊  ÿ  Ÿ  ⁄  €  ‹  ›  ﬁ  ﬂ  ‡  ·  ‚  „  ‰  Â  Ê  Á  Ë  È  Ú  Û  Ù  ı  ˆ  ˜  ¯  ˘  ˙  ˚  ¸  ˝  ˛  ˇ                     	  
              /Users/roboos/matlab/fieldtrip/peer/src/ memprofile.c /Users/roboos/matlab/fieldtrip/peer/src/../private/memprofile.o _memprofile_cleanup _memprofile_sample _memprofile _exitFun _mexFunction _mutexmemlist _mutexmemprofile _reftime _memlist _memprofileStatus _memprofileThread announce.c /Users/roboos/matlab/fieldtrip/peer/src/announce.o _frand _cleanup_announce _announce _announce_once discover.c /Users/roboos/matlab/fieldtrip/peer/src/discover.o _cleanup_discover _discover expire.c /Users/roboos/matlab/fieldtrip/peer/src/expire.o _cleanup_expire _expire _check_watchdog extern.c /Users/roboos/matlab/fieldtrip/peer/src/extern.o _syslog_level _condstatus _mutexstatus _mutexappendcount _mutexsocketcount _mutexthreadcount _mutexconnectioncount _mutexhost _mutexpeerlist _mutexjoblist _mutexallowuserlist _mutexrefuseuserlist _mutexallowgrouplist _mutexrefusegrouplist _mutexallowhostlist _mutexrefusehostlist _mutexwatchdog _mutexsmartmem _mutexsmartcpu _mutexprevcpu _mutexsmartshare _udsserverStatus _tcpserverStatus _announceStatus _discoverStatus _expireStatus _appendcount _socketcount _threadcount _connectioncount _host _peerlist _joblist _allowuserlist _refuseuserlist _allowgrouplist _refusegrouplist _allowhostlist _refusehostlist _smartsharelist _watchdog _smartmem _smartcpu _prevcpu _smartshare peerinit.c /Users/roboos/matlab/fieldtrip/peer/src/peerinit.o _hash _peerinit _peerexit util.c /Users/roboos/matlab/fieldtrip/peer/src/util.o _threadsleep _bufread _bufwrite _append _jobcount _peercount _hoststatus _clear_peerlist _clear_joblist _clear_smartsharelist _clear_allowuserlist _clear_allowgrouplist _clear_allowhostlist _clear_refuseuserlist _clear_refusegrouplist _clear_refusehostlist _check_datatypes _getmem udsserver.c /Users/roboos/matlab/fieldtrip/peer/src/udsserver.o _cleanup_udsserver _udsserver tcpserver.c /Users/roboos/matlab/fieldtrip/peer/src/tcpserver.o _cleanup_tcpserver _tcpserver __OSSwapInt16 /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/include/libkern/i386/_OSByteOrder.h tcpsocket.c /Users/roboos/matlab/fieldtrip/peer/src/tcpsocket.o _cleanup_tcpsocket _tcpsocket security.c /Users/roboos/matlab/fieldtrip/peer/src/security.o _security_check _ismember_userlist _ismember_grouplist _ismember_hostlist localhost.c /Users/roboos/matlab/fieldtrip/peer/src/localhost.o _check_localhost smartshare.c /Users/roboos/matlab/fieldtrip/peer/src/smartshare.o _smartshare_reset _smartshare_check _smartshare_history smartmem.c /Users/roboos/matlab/fieldtrip/peer/src/smartmem.o _smartmem_info _smartmem_update smartcpu.c /U