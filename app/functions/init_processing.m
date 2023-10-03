function app_properties = init_processing()
%% Printing data information
app_properties = jsondecode(fileread(strcat('app/properties.json')));
disp(strcat("-->> Name:",app_properties.generals.name));
disp(strcat("-->> Version:",app_properties.generals.version));
disp(strcat("-->> Version date:",app_properties.generals.version_date));
disp("=====================================================================");

%% ------------ Checking MatLab compatibility ----------------
if(app_properties.check_matlab_version)
    disp('-->> Checking installed matlab version');
    if(~check_matlab_version())
        return;
    end
end

%% ------------  Checking updates --------------------------
if(app_properties.check_app_update)
    disp('-->> Checking latest project version');
    if(isequal(check_version,true))
        return;
    end
end

% Remove fieldtrip path for override functions 
warning off;
rmpath(genpath(fullfile('external/fieldtrip')));
warning on;

end

