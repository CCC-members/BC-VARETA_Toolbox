function [subject,checked,error_msg_array] = ischecked_subject_data(subject_file_info,properties)
checked         = true;
subject         = struct;
subject_info    = load(fullfile(subject_file_info.folder,subject_file_info.name));
subject.name    = subject_info.name;
root_dir        = subject_file_info.folder;
error_msg_array = [];

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
    subject_info.surf_dir       = replace(subject_info.surf_dir,'\','/');
    subject_info.scalp_dir      = replace(subject_info.scalp_dir,'\','/');
    subject_info.innerskull_dir = replace(subject_info.innerskull_dir,'\','/');
    subject_info.outerskull_dir = replace(subject_info.outerskull_dir,'\','/');
end
disp('=================================================================');
disp(strcat('BC-V-->>Processing subject:',subject_info.name));
disp('-----------------------------------------------------------------');

if(~isfile(fullfile(root_dir,subject_info.leadfield_dir)))
    error_msg = strcat("The leadfield file do not exist");
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
    subject.name        = subject_info.name;
    subject.modality    = subject_info.modality;
    
    surf        = load(fullfile(root_dir,subject_info.surf_dir));
    channel     = load(fullfile(root_dir,subject_info.channel_dir));
    scalp       = load(fullfile(root_dir,subject_info.scalp_dir));
    innerskull  = load(fullfile(root_dir,subject_info.innerskull_dir));
    outerskull  = load(fullfile(root_dir,subject_info.outerskull_dir));
    
    if(~isstruct(subject_info.leadfield_dir))
        checked         = false;
        error_msg       = strcat("The leadfield file do not have a correct format. It have to be a structure.");
        error_msg_array = [error_msg_array; error_msg];
        return;
    end    
    
    % Matching leadfield and cortex
    sources = properties.general_params.data_resolution.sources;    
    method  = properties.general_params.data_resolution.method;
    if(isequal(sources,'default'))
        subject.Scortex     = surf.Sc(surf.iCortex);        
        leadfield = load(fullfile(root_dir,subject_info.leadfield_dir.path));            
    else
        if(ischar(sources))
            sources = str2double(sources);
        end
        for i=1:length(surf.Sc)
            if(isequal(size(surf.Sc(i).Vertices,1),sources))
                subject.Scortex = surf.Sc(i);
            end           
        end
        if(isfield(subject,'Scortex'))
            for i=1:length(subject_info.leadfield_dir)                
                leadfield = load(fullfile(root_dir,subject_info.leadfield_dir));
                if(isequal(size(leadfield.Ke,2)/3,sources) && contains(leadfield.Comment,method))
                    break;
                else
                    leadfield = [];
                end
            end
            if(isempty(leadfield))
                checked         = false;
                error_msg       = strcat("There is not any leadfiel with the specified sources or method defined in <<bcv_general_params>>.");
                error_msg_array = [error_msg_array; error_msg];
                return;
            end
        else
        checked         = false;
        error_msg       = strcat("There is not any surface on data with the specified number of vertices defined in <<bcv_general_params>>.");
        error_msg_array = [error_msg_array; error_msg];
        return; 
        end
    end
    subject.Ke              = leadfield.Ke;
    subject.GridOrient      = leadfield.GridOrient;
    subject.GridAtlas       = leadfield.GridAtlas;    
    subject.sub_to_FSAve    = surf.sub_to_FSAve;  
    
    if(~isequal(size(subject.Ke,2)/3,size(subject.Scortex.Vertices,1)))
        checked = false;
        error_msg = strcat("The selected cortex and leadfield do not have the same number of vertices.");
        error_msg_array = [error_msg_array; error_msg];
        return;
    end
    if(isempty(subject.Scortex.Atlas(subject.Scortex.iAtlas).Scouts) && properties.activation_params.IsParcel)        
        checked = false;
        error_msg = strcat("The selected cortex do not have Atlas or the iAtlas is wrong and IsParcel activation param is true.");
        error_msg_array = [error_msg_array; error_msg];
        return;
    end
    
    subject.Shead   = scalp;
    subject.Cdata   = channel;    
    subject.Sinn    = innerskull;    
    subject.Sout    = outerskull;
    
    MEEG = load(fullfile(root_dir,subject_info.meeg_dir));
    if(isequal(subject.modality,'EEG'))
        subject.data = MEEG.data;
    else
        if(properties.general_params.run_by_trial.value)           
            subject.data = MEEG.data;          
        else
            subject.data = MEEG.data;
            if(iscell(MEEG.data))
                subject.data = cell2mat(MEEG.data(1,1:end));
            end
        end
    end
end
end