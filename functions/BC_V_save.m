function subject = BC_V_save(properties,subject,process,varargin)
%% Preparing parameters
for k = 4 : nargin
    eval([inputname(k) '=  varargin{k-3};']);
end


disp('-->> Saving file')
if(isequal(process,'common'))
    %%
    %% Saving general variables for analysis
    %%
    pathname_common                                     = fullfile(subject.subject_path,'Common');
    if(~isfolder(pathname_common))
        mkdir(pathname_common);
    end
    Sscalp                                              = subject.Shead;
    file_name                                           = strcat('Sscalp.mat');
    reference_path                                      = strsplit(pathname_common,subject.name);
    subject.BC_V_info.common(1).Comment                 = 'Surfaces Scalp';
    subject.BC_V_info.common(1).Ref_path                = strrep(reference_path{2},'\','/');
    subject.BC_V_info.common(1).Name                    = file_name;
    disp(strcat("File: ", file_name));
    save(fullfile(pathname_common ,file_name ),'-struct','Sscalp');

    Souter                                              = subject.Sout;
    file_name                                           = strcat('Souterskull.mat');
    reference_path                                      = strsplit(pathname_common,subject.name);
    subject.BC_V_info.common(2).Comment                 = 'Surfaces Outerskull';
    subject.BC_V_info.common(2).Ref_path                = strrep(reference_path{2},'\','/');
    subject.BC_V_info.common(2).Name                    = file_name;
    disp(strcat("File: ", file_name));
    save(fullfile(pathname_common ,file_name ),'-struct','Souter');

    Sinner                                              = subject.Sinn;
    file_name                                           = strcat('Sinnerskull.mat');
    reference_path                                      = strsplit(pathname_common,subject.name);
    subject.BC_V_info.common(3).Comment                 = 'Surfaces Innerskull';
    subject.BC_V_info.common(3).Ref_path                = strrep(reference_path{2},'\','/');
    subject.BC_V_info.common(3).Name                    = file_name;
    disp(strcat("File: ", file_name));
    save(fullfile(pathname_common ,file_name ),'-struct','Sinner');

    cortex                                              = subject.Scortex;
    file_name                                           = strcat('Cortex.mat');
    reference_path                                      = strsplit(pathname_common,subject.name);
    subject.BC_V_info.common(4).Comment                 = 'Surfaces Cortex';
    subject.BC_V_info.common(4).Ref_path                = strrep(reference_path{2},'\','/');
    subject.BC_V_info.common(4).Name                    = file_name;
    disp(strcat("File: ", file_name));
    save(fullfile(pathname_common ,file_name ),'-struct','cortex');

    Channels                                            = subject.Cdata;
    file_name                                           = strcat('Channels.mat');
    reference_path                                      = strsplit(pathname_common,subject.name);
    subject.BC_V_info.common(5).Comment                 = 'Channels data';
    subject.BC_V_info.common(5).Ref_path                = strrep(reference_path{2},'\','/');
    subject.BC_V_info.common(5).Name                    = file_name;
    disp(strcat("File: ", file_name));
    save(fullfile(pathname_common ,file_name ),'-struct','Channels');

    Leadfiled                                           = subject.Headmodel;
    file_name                                           = strcat('Leadfield.mat');
    reference_path                                      = strsplit(pathname_common,subject.name);
    subject.BC_V_info.common(6).Comment                 = 'Leadfiled';
    subject.BC_V_info.common(6).Ref_path                = strrep(reference_path{2},'\','/');
    subject.BC_V_info.common(6).Name                    = file_name;
    disp(strcat("File: ", file_name));
    save(fullfile(pathname_common ,file_name ),'-struct','Leadfiled');
end
if(isequal(process,'fuctional'))
    if(properties.general_params.run_by_trial.value)
        trial = trial_info.id;
        trial_name = trial_info.name;
        pathname                                            = fullfile(subject.subject_path,trial_name);
    else
        pathname                                            = subject.subject_path;
    end
    %%
    %% Saving general variables for sensor level
    %%
    pathname                                            = fullfile(pathname,'Generals','Functional');
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    if(properties.general_params.run_by_trial.value)
        meeg                                                = subject.MEEG;
        file_name                                           = strcat('MEEG.mat');
        reference_path                                      = strsplit(pathname,subject.name);
        subject.BC_V_info.trials(trial).generals(1).Comment    = 'MEEG data';
        subject.BC_V_info.trials(trial).generals(1).Ref_path   = strrep(reference_path{2},'\','/');
        subject.BC_V_info.trials(trial).generals(1).Name       = file_name;
        disp(strcat("File: ", file_name));
        save(fullfile(pathname ,file_name ),'-struct','meeg');

        file_name                                           = strcat('Data_spectrum.mat');
        reference_path                                      = strsplit(pathname,subject.name);
        subject.BC_V_info.trials(trial).generals(2).Comment  = 'PSD';
        subject.BC_V_info.trials(trial).generals(2).Ref_path = strrep(reference_path{2},'\','/');
        subject.BC_V_info.trials(trial).generals(2).Name     = file_name;
        disp(strcat("File: ", file_name));
        parsave(fullfile(pathname ,file_name ),Svv_channel,PSD);
    else
        meeg                                                = subject.MEEG;
        file_name                                           = strcat('MEEG.mat');
        reference_path                                      = strsplit(pathname,subject.name);
        subject.BC_V_info.generals(1).Comment                 = 'MEEG data';
        subject.BC_V_info.generals(1).Ref_path                = strrep(reference_path{2},'\','/');
        subject.BC_V_info.generals(1).Name                    = file_name;
        disp(strcat("File: ", file_name));
        save(fullfile(pathname ,file_name ),'-struct','meeg');

        file_name                                           = strcat('Data_spectrum.mat');
        reference_path                                      = strsplit(pathname,subject.name);
        subject.BC_V_info.generals(2).Comment               = 'PSD';
        subject.BC_V_info.generals(2).Ref_path              = strrep(reference_path{2},'\','/');
        subject.BC_V_info.generals(2).Name                  = file_name;
        disp(strcat("File: ", file_name));
        parsave(fullfile(pathname ,file_name ),Svv_channel,PSD);
    end
end

%%
%% Sensor level outputs 
%%
if(isequal(process,'sensor'))
    if(properties.general_params.run_by_trial.value)
        trial = trial_info.id;
        trial_name = trial_info.name;
        pathname                                            = fullfile(subject.subject_path,trial_name);
    else
        pathname                                            = subject.subject_path;
    end
    level                  = 'Sensor_level';
    pathname = fullfile(pathname,level,band.name);
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    file_name                                           = strcat('Sensor_level_',band.str_band,'.mat');
    reference_path                                      = strsplit(pathname,subject.name);

    [~,band_name,~]                                     = fileparts(reference_path{2});
    if(properties.general_params.run_by_trial.value)
        subject.BC_V_info.trials(trial).sensor_level(pos).Comment         = 'Sensor_level';
        subject.BC_V_info.trials(trial).sensor_level(pos).Band            = band_name;
        subject.BC_V_info.trials(trial).sensor_level(pos).Freq            = char(band.str_band);
        subject.BC_V_info.trials(trial).sensor_level(pos).Ref_path        = strrep(reference_path{2},'\','/');
        subject.BC_V_info.trials(trial).sensor_level(pos).Name            = file_name;
        disp(strcat("File: ", file_name));
        parsave(fullfile(pathname ,file_name ),Svv,peak_pos,Nseg,band);
    else
        subject.BC_V_info.sensor_level(pos).Comment         = 'Sensor_level';
        subject.BC_V_info.sensor_level(pos).Band            = band_name;
        subject.BC_V_info.sensor_level(pos).Freq            = char(band.str_band);
        subject.BC_V_info.sensor_level(pos).Ref_path        = strrep(reference_path{2},'\','/');
        subject.BC_V_info.sensor_level(pos).Name            = file_name;
        disp(strcat("File: ", file_name));
        parsave(fullfile(pathname ,file_name ),Svv,peak_pos,Nseg,band);
    end
end

%%
%% Sensor mean 
%%
if(isequal(process,'sensor_mean'))
    pathname                                            = subject.subject_path;

    level                  = 'Sensor_level';
    pathname = fullfile(pathname,level,band.name);
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    file_name                                           = strcat('Sensor_level_',band.str_band,'.mat');
    reference_path                                      = strsplit(pathname,subject.name);
    [~,band_name,~]                                     = fileparts(reference_path{2});
    subject.BC_V_info.sensor_level(pos).Comment         = 'Sensor_level';
    subject.BC_V_info.sensor_level(pos).Band            = band_name;
    subject.BC_V_info.sensor_level(pos).Freq            = char(band.str_band);
    subject.BC_V_info.sensor_level(pos).Ref_path        = strrep(reference_path{2},'\','/');
    subject.BC_V_info.sensor_level(pos).Name            = file_name;
    disp(strcat("File: ", file_name));
    parsave(fullfile(pathname ,file_name ),Svv,peak_pos,Nseg,band);    
end

if(isequal(process,'level1'))
    disp('-->> Saving BC-VARETA Information file.')
    subject.BC_V_info.Processes(1).name         = 'Sensor_level';
    subject.BC_V_info.Processes(1).completed    = true;
    BC_V_info                                   = subject.BC_V_info;
    saveJSON(BC_V_info,fullfile(subject.subject_path ,strcat(subject.name,'.json')));
end

%%
%% Activation priors
%%
if(isequal(process,'a_priors'))
    file_name                                           = strcat('W.mat');
    W                                                   = subject.W;
    pathname                                            = fullfile(subject.subject_path,'Generals','Structural','sSSBL');
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    reference_path                                      = strsplit(pathname,subject.name);
    subject.BC_V_info.generals(3).Comment               = 'Generals';
    subject.BC_V_info.generals(3).Ref_path              = strrep(reference_path{2},'\','/');
    subject.BC_V_info.generals(3).Name                  = file_name;
    disp(strcat("File: ", file_name));
    if(getGlobalGuimode)
        dlg = msgbox('Save operation in progress...');
    end
    parsave(fullfile(pathname ,file_name ),W);
    if exist('dlg','var')
        delete(dlg);
    end
end
if(isequal(process,'activation'))
    level                                               = 'Activation_level';    
    if(properties.general_params.run_by_trial.value && ~isequal(properties.general_params.run_by_trial.level,'sensor'))
        trial = trial_info.id;
        trial_name = trial_info.name;
        pathname                                            = fullfile(subject.subject_path,trial_name);
    else
        pathname                                            = subject.subject_path;
    end
    pathname                                            = fullfile(pathname,level);
    file_name                                           = strcat('MEEG_source_',band.str_band,'.mat');
    switch method
        case 'sssblpp'
            pathname                                    = fullfile(pathname,'sSSBLpp',band.name);
            if(~isfolder(pathname))
                mkdir(pathname);
            end
        case 'eloreta'
            pathname                                    = fullfile(pathname,'eLORETA',band.name);
            if(~isfolder(pathname))
                mkdir(pathname);
            end
        case 'lcmv'
            pathname                                    = fullfile(pathname,'LCMV',band.name);
            if(~isfolder(pathname))
                mkdir(pathname);
            end
    end
    disp(strcat("File: ", file_name));
    save(fullfile(pathname ,file_name ),"-struct",'outputs');

    %%
    reference_path = strsplit(pathname,subject.name);
    if(properties.general_params.run_by_trial.value && ~isequal(properties.general_params.run_by_trial.level,'sensor'))
        subject.BC_V_info.trials(trial).activation_level(pos).Comment     = level;
        subject.BC_V_info.trials(trial).activation_level(pos).Band        = band.name;
        subject.BC_V_info.trials(trial).activation_level(pos).Method      = lower(method);
        subject.BC_V_info.trials(trial).activation_level(pos).Freq        = char(band.str_band);
        subject.BC_V_info.trials(trial).activation_level(pos).Ref_path    = strrep(reference_path{2},'\','/');
        subject.BC_V_info.trials(trial).activation_level(pos).Name        = file_name;
    else
        subject.BC_V_info.activation_level(pos).Comment     = level;
        subject.BC_V_info.activation_level(pos).Band        = band.name;
        subject.BC_V_info.activation_level(pos).Method      = lower(method);
        subject.BC_V_info.activation_level(pos).Freq        = char(band.str_band);
        subject.BC_V_info.activation_level(pos).Ref_path    = strrep(reference_path{2},'\','/');
        subject.BC_V_info.activation_level(pos).Name        = file_name;
    end

end

if(isequal(process,'level2'))
    disp('-->> Saving BC-VARETA Information file.')
    subject.BC_V_info.Processes(2).name      = 'Activation_level';
    subject.BC_V_info.Processes(2).completed = true;
    BC_V_info                                = subject.BC_V_info;
    saveJSON(BC_V_info,fullfile(subject.subject_path ,strcat(subject.name,'.json')));
end

%%
%% Connectvity priors
%%
if(isequal(process,'c_priors'))
    file_name                                           = strcat('W.mat');
    W                                                   = subject.W;
    pathname                                            = fullfile(subject.subject_path,'Generals','Structural','HiGSS');
    if(~isfolder(pathname))
        mkdir(pathname);
    end

    reference_path                                      = strsplit(pathname,subject.name);
    subject.BC_V_info.generals(4).Comment               = 'Generals';
    subject.BC_V_info.generals(4).Ref_path              = strrep(reference_path{2},'\','/');
    subject.BC_V_info.generals(4).Name                  = file_name;

    disp(strcat("File: ", file_name));
    if(getGlobalGuimode)
        dlg = msgbox('Save operation in progress...');
    end
    parsave(fullfile(pathname ,file_name ),W);
    if exist('dlg','var')
        delete(dlg);
    end
end
%%
%% Connectvity level Outputs
%%
if(isequal(process,'connectivity'))
    level                                               = 'Connectivity_level';
    if(properties.general_params.run_by_trial.value)
        pathname                                            = fullfile(subject.subject_path,properties.trial_name);
    else
        pathname                                            = subject.subject_path;
    end
    pathname                                            = fullfile(pathname,level);
    file_name                                           = strcat('MEEG_source_',band.str_band,'.mat');
    switch method
        case 'higgs'
            pathname                                    = fullfile(pathname,'HiGGS',band.name);
            if(~isfolder(pathname))
                mkdir(pathname);
            end
        case 'hg_lasso'
            pathname                                    = fullfile(pathname,'hg_LASSO',band.name);
            if(~isfolder(pathname))
                mkdir(pathname);
            end
    end
    disp(strcat("File: ", file_name));
    save(fullfile(pathname ,file_name ),"-struct",'outputs');

    %%
    reference_path = strsplit(pathname,subject.name);
    if(properties.general_params.run_by_trial.value)
        subject.BC_V_info.trials(trial).connectivity_level(pos).Comment     = level;
        subject.BC_V_info.trials(trial).connectivity_level(pos).Band        = band.name;
        subject.BC_V_info.trials(trial).connectivity_level(pos).Method      = lower(method);
        subject.BC_V_info.trials(trial).connectivity_level(pos).Freq        = char(band.str_band);
        subject.BC_V_info.trials(trial).connectivity_level(pos).Ref_path    = strrep(reference_path{2},'\','/');
        subject.BC_V_info.trials(trial).connectivity_level(pos).Name        = file_name;
    else
        subject.BC_V_info.connectivity_level(pos).Comment     = level;
        subject.BC_V_info.connectivity_level(pos).Band        = band.name;
        subject.BC_V_info.connectivity_level(pos).Method      = lower(method);
        subject.BC_V_info.connectivity_level(pos).Freq        = char(band.str_band);
        subject.BC_V_info.connectivity_level(pos).Ref_path    = strrep(reference_path{2},'\','/');
        subject.BC_V_info.connectivity_level(pos).Name        = file_name;
    end
end

if(isequal(process,'level3'))
    disp('-->> Saving BC-VARETA Information file.')
    subject.BC_V_info.Processes(3).name         = 'Connectivity_level';
    subject.BC_V_info.Processes(3).completed    = true;
    BC_V_info                                   = subject.BC_V_info;
    saveJSON(BC_V_info,fullfile(subject.subject_path ,strcat(subject.name,'.json')));
end
