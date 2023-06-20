function subject = BC_V_save(properties,subject,process,varargin)
%% Preparing parameters
for k = 4 : nargin
    eval([inputname(k) '=  varargin{k-3};']);
end
if(properties.general_params.run_by_trial.value)
    pathname                                            = fullfile(subject.subject_path,properties.trial_name);
else
    pathname                                            = subject.subject_path;
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
    file_name                                           = strcat('Souterskull.mat');
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

    meeg                                                = subject.MEEG;
    file_name                                           = strcat('MEEG.mat');
    reference_path                                      = strsplit(pathname_common,subject.name);
    subject.BC_V_info.common(5).Comment                 = 'MEEG data';
    subject.BC_V_info.common(5).Ref_path                = strrep(reference_path{2},'\','/');
    subject.BC_V_info.common(5).Name                    = file_name;
    disp(strcat("File: ", file_name));
    save(fullfile(pathname_common ,file_name ),'-struct','meeg');

    Channels                                            = subject.Cdata;
    file_name                                           = strcat('Channels.mat');
    reference_path                                      = strsplit(pathname_common,subject.name);
    subject.BC_V_info.common(6).Comment                 = 'Channels data';
    subject.BC_V_info.common(6).Ref_path                = strrep(reference_path{2},'\','/');
    subject.BC_V_info.common(6).Name                    = file_name;
    disp(strcat("File: ", file_name));
    save(fullfile(pathname_common ,file_name ),'-struct','Channels');
end
if(isequal(process,'fuctional'))
    %%
    %% Saving general variables for sensor level
    %%    
    pathname                                            = fullfile(pathname,'Generals','Functional');   
    if(~isfolder(pathname))
        mkdir(pathname);
    end    
    file_name                                           = strcat('Data_spectrum.mat');
    reference_path                                      = strsplit(pathname,subject.name);
    subject.BC_V_info.generals(1).Comment               = 'Generals';
    subject.BC_V_info.generals(1).Ref_path              = strrep(reference_path{2},'\','/');
    subject.BC_V_info.generals(1).Name                  = file_name;
    disp(strcat("File: ", file_name));
    parsave(fullfile(pathname ,file_name ),Svv_channel,PSD);
end
if(isequal(process,'sensor'))
    level                  = 'Sensor_level';
    pathname = fullfile(pathname,level,band.name);     
    if(~isfolder(pathname))
        mkdir(pathname);
    end  
    file_name                                           = strcat('Sensor_level_',band.str_band,'.mat');
    reference_path                                      = strsplit(pathname,subject.name);
    subject.BC_V_info.sensor_level(pos).Comment         = 'Sensor_level';
    [~,band_name,~]                                     = fileparts(reference_path{2});
    subject.BC_V_info.sensor_level(pos).Band            = band_name;
    subject.BC_V_info.sensor_level(pos).Freq            = char(band.str_band);
    subject.BC_V_info.sensor_level(pos).Ref_path        = strrep(reference_path{2},'\','/');
    subject.BC_V_info.sensor_level(pos).Name            = file_name;
    disp(strcat("File: ", file_name));
    parsave(fullfile(pathname ,file_name ),Svv,peak_pos,Nseg,band);
end
if(isequal(process,'priors'))
    file_name                                           = strcat('W.mat');
    W                                                   = subject.W;
    if(properties.general_params.run_by_trial.value)
        pathname                                        = fullfile(subject.subject_path,properties.trial_name,'Generals','Structural','sSSBL');
        reference_path                                  = strsplit(pathname,subject.name);
        if(~isfolder(pathname))
            mkdir(pathname);
        end
    else
        pathname                                        = fullfile(subject.subject_path,'Generals','Structural','sSSBL');
        reference_path                                  = strsplit(pathname,subject.name);
        if(~isfolder(pathname))
            mkdir(pathname);
        end
    end
    subject.BC_V_info.generals(2).Comment               = 'Generals';
    subject.BC_V_info.generals(2).Ref_path              = strrep(reference_path{2},'\','/');
    subject.BC_V_info.generals(2).Name                  = file_name;
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
    pathname                                            = fullfile(pathname,level);    
    file_name                                           = strcat('MEEG_source_',band.str_band,'.mat');
    switch method
        case 'sssblpp'
            pathname                                    = fullfile(pathname,'sSSBLpp',band.name);
            if(~isfolder(pathname))
                mkdir(pathname);
            end            
            s2j                                         = outputs.s2j;
            sigma2j                                     = outputs.sigma2j;
            T                                           = outputs.T;
            scaleSvv                                    = outputs.scaleSvv;
            scaleKe                                     = outputs.scaleKe;
            stat                                        = outputs.stat;
            J                                           = outputs.J;
            Jsp                                         = outputs.Jsp;
            indms                                       = outputs.indms;
            J_FSAve                                     = outputs.J_FSAve;
            Jsp_FSAve                                   = outputs.Jsp_FSAve;
            disp(strcat("File: ", file_name));
            parsave(fullfile(pathname ,file_name ),s2j,sigma2j,T,scaleSvv,scaleKe,stat,J,Jsp,indms,J_FSAve,Jsp_FSAve);
        case 'eloreta'
            pathname                                    = fullfile(pathname,'eLORETA',band.name);            
            if(~isfolder(pathname))
                mkdir(pathname);
            end
            s2j                                         = outputs.s2j;
            T                                           = outputs.T;
            scaleSvv                                    = outputs.scaleSvv;
            scaleKe                                     = outputs.scaleKe;
            stat                                        = outputs.stat;
            J                                           = outputs.J;
            Jsp                                         = outputs.Jsp;
            indms                                       = outputs.indms;
            J_FSAve                                     = outputs.J_FSAve;
            Jsp_FSAve                                   = outputs.Jsp_FSAve;
            disp(strcat("File: ", file_name));
            parsave(fullfile(pathname ,file_name ),s2j,sigma2j_post,T,scaleSvv,scaleKe,stat,J,Jsp,indms,J_FSAve,Jsp_FSAve);
        case 'lcmv'
            pathname                                    = fullfile(pathname,'LCMV',band.name);            
            if(~isfolder(subject.pathname))
                mkdir(subject.pathname);
            end
            s2j                                         = outputs.s2j;
            T                                           = outputs.T;
            scaleSvv                                    = outputs.scaleSvv;
            scaleKe                                     = outputs.scaleKe;
            stat                                        = outputs.stat;
            J                                           = outputs.J;
            Jsp                                         = outputs.Jsp;
            indms                                       = outputs.indms;
            J_FSAve                                     = outputs.J_FSAve;
            Jsp_FSAve                                   = outputs.Jsp_FSAve;
            disp(strcat("File: ", file_name));
            parsave(fullfile(pathname ,file_name ),s2j,sigma_post,T,scaleSvv,scaleKe,stat,J,Jsp,indms,J_FSAve,Jsp_FSAve);
    end
    %% 
    reference_path = strsplit(pathname,subject.name);    
    subject.BC_V_info.activation_level(pos).Comment     = 'Activation_level';
    subject.BC_V_info.activation_level(pos).Band        = band.name;
    subject.BC_V_info.activation_level(pos).Method      = lower(method);
    subject.BC_V_info.activation_level(pos).Freq        = char(band.str_band);
    subject.BC_V_info.activation_level(pos).Ref_path    = strrep(reference_path{2},'\','/');
    subject.BC_V_info.activation_level(pos).Name        = file_name;

end

if(isequal(process,'level1'))
    disp('-->> Saving BC-VARETA Information file.')
    subject.BC_V_info.Processes(1).name         = 'Sensor_level';
    subject.BC_V_info.Processes(1).completed    = true;
    BC_V_info                                   = subject.BC_V_info;
    save(fullfile(subject.subject_path ,'BC_V_info.mat'),'-struct','BC_V_info');
end
if(isequal(process,'level2'))
    disp('-->> Saving BC-VARETA Information file.')
    subject.BC_V_info.Processes(2).name      = 'Activation_level';
    subject.BC_V_info.Processes(2).completed = true;
    BC_V_info                                = subject.BC_V_info;
    save(fullfile(subject.subject_path ,'BC_V_info.mat'),'-struct','BC_V_info');
end
if(isequal(process,'level3'))
    disp('-->> Saving BC-VARETA Information file.')
    subject.BC_V_info.Processes(3).name         = 'Connectivity_level';
    subject.BC_V_info.Processes(3).completed    = true;
    BC_V_info                                   = subject.BC_V_info;
    save(fullfile(subject.subject_path ,'BC_V_info.mat'),'-struct','BC_V_info');
end
