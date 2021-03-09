function [output_subject] = create_data_structure(root_path,subject_name)
%CREATE_DATA_STRUCTURE Summary of this function goes here
%   Detailed explanation goes here


if(nargin == 1)
    guiHandle = total_subjects_guide;
    
    disp('------Waitintg for frequency_bands------');
    uiwait(guiHandle.UIFigure);
    
    if(guiHandle.canceled)
        delete(guiHandle);
        output_subject = null;
        return;
    else
        if(~isfolder(strcat(root_path,filesep,'Data')))
            mkdir(root_path,'Data');            
        end
        for i = 1:guiHandle.total_subjects
            subject_name = strcat( 'Subject#',string(i));
            mkdir(strcat(root_path,filesep,'Data'),subject_name);
            mkdir(strcat(root_path,filesep,'Data',filesep,subject_name),'eeg');
            mkdir(strcat(root_path,filesep,'Data',filesep,subject_name),'leadfield');
            mkdir(strcat(root_path,filesep,'Data',filesep,subject_name),'scalp');
            mkdir(strcat(root_path,filesep,'Data',filesep,subject_name),'surf');
        end
        delete(guiHandle);
    end
else
   
    output_subject  = strcat(root_path,filesep,subject_name);
    if(~isfolder(output_subject))
        mkdir(output_subject);
        mkdir(strcat(output_subject),'eeg');
        mkdir(strcat(output_subject),'leadfield');
        mkdir(strcat(output_subject),'scalp');
        mkdir(strcat(output_subject),'surf');
    end
   
end

