function [subject,checked,error_msg_array] = checked_subject_data(subject_file_info,properties)
checked         = true;
subject         = struct;
subject_info    = load(fullfile(subject_file_info.folder,subject_file_info.name));
subject.name    = subject_info.name;
root_dir        = subject_file_info.folder;
error_msg_array = [];

if(~subject_info.completed)
    error_msg       = strcat("The subject information is not completed");
    error_msg_array = [error_msg_array; error_msg];
    checked         = false;
    return;
end

if(~isfield(subject_info,'name')...
        || ~isfield(subject_info,'modality')...
        || ~isfield(subject_info,'leadfield_dir')...
        || ~isfield(subject_info,'surf_dir')...
        || ~isfield(subject_info,'scalp_dir')...
        || ~isfield(subject_info,'innerskull_dir')...
        || ~isfield(subject_info,'outerskull_dir')...
        || ~isfield(subject_info,'channel_dir')...
        || ~isfield(subject_info,'meeg_dir'))
    checked = false;
    return;
else 
    subject_info.meeg_dir       = replace(subject_info.meeg_dir,'\','/');
    subject_info.channel_dir    = replace(subject_info.channel_dir,'\','/');
    subject_info.leadfield_dir  = replace(subject_info.leadfield_dir,'\','/');
    subject_info.surf_dir       = replace(subject_info.surf_dir,'\','/');
    subject_info.scalp_dir      = replace(subject_info.scalp_dir,'\','/');
    subject_info.innerskull_dir = replace(subject_info.innerskull_dir,'\','/');
    subject_info.outerskull_dir = replace(subject_info.outerskull_dir,'\','/');
end
disp("=====================================================================");
disp(strcat('BC-V-->>Processing subject:',subject_info.name));
disp("=====================================================================");

if(~isfile(fullfile(root_dir,subject_info.meeg_dir)))
    error_msg       = strcat("The meeg file do not exist");
    error_msg_array = [error_msg_array; error_msg];
end
if(~isfile(fullfile(root_dir,subject_info.leadfield_dir)))
    error_msg       = strcat("The leadfield file do not exist");
    error_msg_array = [error_msg_array; error_msg];
end
if(~isfile(fullfile(root_dir,subject_info.channel_dir)))
    error_msg       = strcat("The channel file do not exist");
    error_msg_array = [error_msg_array; error_msg];
end
if(~isfile(fullfile(root_dir,subject_info.surf_dir)))
    error_msg       = strcat("The surface file do not exist");
    error_msg_array = [error_msg_array; error_msg];
end
if(~isfile(fullfile(root_dir,subject_info.scalp_dir)))
    error_msg       = strcat("The scalp file do not exist");
    error_msg_array = [error_msg_array; error_msg];
end
if(~isfile(fullfile(root_dir,subject_info.innerskull_dir)))
    error_msg       = strcat("The innerskull file do not exist");
    error_msg_array = [error_msg_array; error_msg];
end
if(~isfile(fullfile(root_dir,subject_info.outerskull_dir)))
    error_msg       = strcat("The outerskull file do not exist");
    error_msg_array = [error_msg_array; error_msg];
end

if(isempty(error_msg_array))    
    Scortex     = load(fullfile(root_dir,subject_info.surf_dir));
    Cdata       = load(fullfile(root_dir,subject_info.channel_dir));
    Shead       = load(fullfile(root_dir,subject_info.scalp_dir));
    Sinn        = load(fullfile(root_dir,subject_info.innerskull_dir));
    Sout        = load(fullfile(root_dir,subject_info.outerskull_dir));
    Headmodels  = load(fullfile(root_dir,subject_info.leadfield_dir));  
    MEEG        = load(fullfile(root_dir,subject_info.meeg_dir));

    if(properties.general_params.run_by_trial.value && ~iscell(MEEG.data))
        error_msg       = strcat("Process run by trial and the data format is not a cellarray");
        error_msg_array = [error_msg_array; error_msg];
    end
    
    subject.name            = subject_info.name;
    subject.modality        = subject_info.modality;  
    subject.Headmodel       = Headmodels.HeadModel(Headmodels.iHeadModel);   
    subject.Scortex         = Scortex.Sc(Scortex.iCortex); 
    if(isfield(Scortex,'sub_to_FSAve'))
        subject.sub_to_FSAve    = Scortex.sub_to_FSAve;
    else
        subject.sub_to_FSAve    = [];
    end
    subject.Shead           = Shead;
    subject.Cdata           = Cdata;    
    subject.Sinn            = Sinn;    
    subject.Sout            = Sout;
    subject.MEEG            = MEEG;      
end
end