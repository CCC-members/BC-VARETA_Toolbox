function update_version()
try
    if(isnetav())
        disp('-->> Comparing local and master version');
        if(local.generals.version_number < online.generals.version_number)
            %% Download lasted version
            filename = strcat('BCV_last_version.zip');
            disp(strcat("-->> Downloading latest version......."));
            url = local.generals.base_url;
            matlab.net.http.HTTPOptions.VerifyServerName = false;
            options = weboptions('Timeout',Inf,'RequestMethod','auto');
            websave(filename,url,options);

            %% Unzip lasted version
            pause(1);
            disp(strcat("-->> Unpacking files..."));

            unzip(filename,pwd);
            pause(1);
            delete(filename);
            % Backuping configuration files
            mkdir('tmp');
            copyfile('bcv_properties/general_params.json','tmp/');
            copyfile('bcv_properties/sensor_params.json','tmp/');
            copyfile('bcv_properties/activation_params.json','tmp/');
            copyfile('bcv_properties/connectivity_params.json','tmp/');

            movefile( strcat('BC-VARETA_Toolbox-master',filesep,'*'), pwd);
            rmdir BC-VARETA_Toolbox-master ;

            % Restoring configuration files
            movefile('tmp/general_params.json','bcv_properties/');
            movefile('tmp/sensor_params.json','bcv_properties/');
            movefile('tmp/activation_params.json','bcv_properties/');
            movefile('tmp/connectivity_params.json','bcv_properties/');
            rmdir('tmp');

            disp('-->> The project is already update with the latest version.');
            disp('-->> The process was stoped to refresh all file');
            disp('-->> Please configure the app properties file, before restart the process.');
            pause(5);
            exit;
        else
            disp('-->> Nothing to update');
            disp('-->> You have already the BC-VARETA latest version')
            disp('=================================================================');
        end
    end
catch
    return;
end

end

