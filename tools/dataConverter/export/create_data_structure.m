function [subject_path] = create_data_structure(varargin)
%CREATE_DATA_STRUCTURE Summary of this function goes here
%   Detailed explanation goes here

for i=1:length(varargin)
   eval([inputname(i) '= varargin{i};']); 
end

subject_path = fullfile(base_path,subID);
meeg_path = fullfile(subject_path,'meeg');
leadfield_path = fullfile(subject_path,'leadfield');
scalp_path = fullfile(subject_path,'scalp');
channel_path = fullfile(subject_path,'channel');
surf_path = fullfile(subject_path,'surf');
if(isequal(action,'all'))    
    if(~isfolder(subject_path))
        mkdir(subject_path);
    end    
    if(~isfolder(meeg_path))
        mkdir(meeg_path);
    end
    if(~isfolder(leadfield_path))
        mkdir(leadfield_path);
    end
    if(~isfolder(scalp_path))
        mkdir(scalp_path);
    end
    if(~isfolder(channel_path))
        mkdir(channel_path);
    end
    if(~isfolder(surf_path))
        mkdir(surf_path);
    end    
end
if(isequal(action,'meeg'))    
    if(~isfolder(subject_path))
        mkdir(subject_path);
    end    
    if(~isfolder(meeg_path))
        mkdir(meeg_path);
    end      
end
if(isequal(action,'anat'))
    if(~isfolder(subject_path))
        mkdir(subject_path);
    end        
    if(~isfolder(leadfield_path))
        mkdir(leadfield_path);
    end
    if(~isfolder(scalp_path))
        mkdir(scalp_path);
    end
    if(~isfolder(channel_path))
        mkdir(channel_path);
    end
    if(~isfolder(surf_path))
        mkdir(surf_path);
    end  
end

end

