function [properties] = get_properties()
try
    pred_options = jsondecode(fileread(strcat('bcv_predefinition/pred_properties.json')));
    if(~isequal(pred_options.params.predefinition.option,'default'))
        properties = jsondecode(fileread(strcat('bcv_predefinition/',pred_options.params.predefinition.option,'/properties.json')));
    else
        properties = jsondecode(fileread(strcat('app/properties.json')));
    end
catch ME
    fprintf(2,strcat('\nBC-V-->> Error: Loading the property files: \n'));
    fprintf(2,strcat(ME.message,'\n'));
    fprintf(2,strcat('Cause in file app\properties.json \n'));
    disp('Please verify the json format in the file.');
    properties = 'canceled';
    return;
end
try
    general_params              = jsondecode(fileread(properties.general_params_file.file_path));
    properties.general_params   = general_params;
catch ME
    fprintf(2,strcat('\nBC-V-->> Error: Loading the property files: \n'));
    fprintf(2,strcat(ME.message,'\n'));
    fprintf(2,strcat('Cause in file', properties.general_params_file.file_path , '\n'));
    disp('Please verify the json format in the file.');
    properties = 'canceled';
    return;
end
param_files = properties.module_param_files;
for i=1:length(param_files)
    try        
        module_params                                   = jsondecode(fileread(param_files(i).file_path));
        properties.(param_files(i).module_id)           = module_params;        
    catch ME
        fprintf(2,strcat('\nBC-V-->> Error: Loading the property files: \n'));
        fprintf(2,strcat(ME.message,'\n'));
        fprintf(2,strcat('Cause in file', param_files(i).file_path , '\n'));
        disp('Please verify the json format in the file.');
        properties = 'canceled';
        return;
    end
end
end

