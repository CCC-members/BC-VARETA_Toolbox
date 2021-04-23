function status = set_properties(properties)
%SET_PROPERTIES Summary of this function goes here
%   Detailed explanation goes here
status = true;
disp("-->> Checking properties");
if(~isfolder(properties.general_params.params.bcv_workspace.BCV_input_dir))
    fprintf(2,strcat('\nBC-V-->> Error: The param BCV_input_dir defined on bcv_properties/general_params.json file: \n'));
    disp(properties.general_params.params.bcv_workspace.BCV_input_dir);
    fprintf(2,strcat('It is not a correct adreess directory. \n'));
    disp('Please verify the location path.');
    status = false;
    return;
end
if(~isfolder(properties.general_params.params.bcv_workspace.BCV_work_dir))
    fprintf(2,strcat('\nBC-V-->> Error: The param BCV_work_dir defined on bcv_properties/general_params.json file: \n'));
    disp(properties.general_params.params.bcv_workspace.BCV_work_dir);
    fprintf(2,strcat('It is not a correct adreess directory. \n'));
    disp('Please verify the location path.');
    status = false;
    return;
else
    [status,values] = fileattrib(properties.general_params.params.bcv_workspace.BCV_work_dir);
    if(~values.UserWrite)
        fprintf(2,strcat('\nBC-V-->> Error: The current user do not have write permissions on: \n'));
        disp(properties.general_params.params.bcv_workspace.BCV_work_dir);
        disp('Please check the folder permission.');
        status = false;
        return;
    end
end
frequencies = properties.spectral_params.params.frequencies;
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
pred_folder = strcat('bcv_predefinition/',properties.run_bash_mode.predefinition_params);
if(~isfolder(pred_folder))
   mkdir(pred_folder); 
end
properties.general_params_file.file_path    = strcat(pred_folder,'/general_params.json');
properties.module_param_files(1).file_path  = strcat(pred_folder,'/sensor_params.json');
properties.module_param_files(2).file_path  = strcat(pred_folder,'/activation_params.json');
properties.module_param_files(3).file_path  = strcat(pred_folder,'/connectivity_params.json');
properties.module_param_files(4).file_path  = strcat(pred_folder,'/spectral_params.json');
pred_options                                = jsondecode(fileread(strcat('bcv_predefinition/pred_properties.json')));
pred_options.params.predefinition.option    = properties.run_bash_mode.predefinition_params;

% saving property files
saveJSON(pred_options,strcat('bcv_predefinition/pred_properties.json'));
saveJSON(properties.general_params,strcat(pred_folder,'/general_params.json'));
saveJSON(properties.sensor_params,strcat(pred_folder,'/sensor_params.json'));
saveJSON(properties.activation_params,strcat(pred_folder,'/activation_params.json'));
saveJSON(properties.connectivity_params,strcat(pred_folder,'/connectivity_params.json'));
saveJSON(properties.spectral_params,strcat(pred_folder,'/spectral_params.json'));

properties = rmfield(properties,{'general_params','sensor_params','activation_params','connectivity_params','spectral_params'});
saveJSON(properties,strcat(pred_folder,'/properties.json'));

end

