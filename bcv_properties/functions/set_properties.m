function status = set_properties(properties)
%SET_PROPERTIES Summary of this function goes here
%   Detailed explanation goes here
status = true;
disp("-->> Checking properties");
if(~isfolder(properties.general_params.bcv_workspace.BCV_input_dir))
    fprintf(2,strcat('\nBC-V-->> Error: The param BCV_input_dir defined on bcv_properties/general_params.json file: \n'));
    disp(properties.general_params.bcv_workspace.BCV_input_dir);
    fprintf(2,strcat('It is not a correct adreess directory. \n'));
    disp('Please verify the location path.');
    status = false;
    return;
end
if(~isfolder(properties.general_params.bcv_workspace.BCV_work_dir))
    fprintf(2,strcat('\nBC-V-->> Error: The param BCV_work_dir defined on bcv_properties/general_params.json file: \n'));
    disp(properties.general_params.bcv_workspace.BCV_work_dir);
    fprintf(2,strcat('It is not a correct adreess directory. \n'));
    disp('Please verify the location path.');
    status = false;
    return;
else
    [status,values] = fileattrib(properties.general_params.bcv_workspace.BCV_work_dir);
    if(~values.UserWrite)
        fprintf(2,strcat('\nBC-V-->> Error: The current user do not have write permissions on: \n'));
        disp(properties.general_params.bcv_workspace.BCV_work_dir);
        disp('Please check the folder permission.');
        status = false;
        return;
    end
end
frequencies = properties.sensor_params.frequencies;
for i=1:length(frequencies)
    freq = frequencies(i);
    if(freq.run && freq.f_start > freq.f_end)
        fprintf(2,strcat('\nBC-V-->> Error: The current frequency is not well configured: \n'));
        disp(freq.name);
        disp('Please check the <<f_start>> and <<f_end>> params.');
        status = false;
        return;
    end
end

disp("-->> Setting User properties");

% saving property files
saveJSON(properties.general_params,'bcv_properties/general_params.json');
h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/general_params.json'));
h.smartIndentContents
h.save
h.close
saveJSON(properties.sensor_params,'bcv_properties/sensor_params.json');
h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/sensor_params.json'));
h.smartIndentContents
h.save
h.close
saveJSON(properties.activation_params,'bcv_properties/activation_params.json');
h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/activation_params.json'));
h.smartIndentContents
h.save
h.close
saveJSON(properties.connectivity_params,'bcv_properties/connectivity_params.json');
h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/connectivity_params.json'));
h.smartIndentContents
h.save
h.close
end

