function [properties] = get_properties()
try    
    properties = jsondecode(fileread(strcat('app/properties.json')));    
catch ME
    fprintf(2,strcat('\nBC-V-->> Error: Loading the property files: \n'));
    fprintf(2,strcat(ME.message,'\n'));
    fprintf(2,strcat('Cause in file app\properties.json \n'));
    disp('Please verify the json format in the file.');
    properties = 'canceled';
    return;
end
defaults_param_files = properties.default_param_files;
for i=1:length(defaults_param_files)
    try        
        module_params                                               = jsondecode(fileread(defaults_param_files(i).file_path));
        properties.defaults.(defaults_param_files(i).module_id)     = module_params;        
    catch ME
        fprintf(2,strcat('\nBC-V-->> Error: Loading the property files: \n'));
        fprintf(2,strcat(ME.message,'\n'));
        fprintf(2,strcat('Cause in file', defaults_param_files(i).file_path , '\n'));
        disp('Please verify the json format in the file.');
        properties = 'canceled';
        return;
    end
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
properties.general_params       = properties.general_params.params;
properties.sensor_params        = properties.sensor_params.params;
properties.activation_params    = properties.activation_params.params;
properties.connectivity_params  = properties.connectivity_params.params;
color_map                       = load(properties.general_params.colormap_path);
properties.cmap                 = color_map.cmap;
properties.cmap_a               = color_map.cmap_a;
properties.cmap_c               = color_map.cmap_c;
end

