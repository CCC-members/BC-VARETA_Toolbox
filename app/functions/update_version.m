function status = update_version()
try
    if(isnetav())
        % loading local data        
        local = jsondecode(fileread(strcat('app/properties.json')));
        % finding online data
        url = local.generals.url_check;
        matlab.net.http.HTTPOptions.VerifyServerName = false;
        options = weboptions('ContentType','json','Timeout',Inf,'RequestMethod','auto');
        online = webread(url,options);
        disp('-->> Comparing local and master version');
        if(local.generals.version_number < online.generals.version_number)
            if(getGlobalGuimode())
                fig = uifigure;
                d = uiprogressdlg(fig,'Title','Please Wait',...
                    'Message','Downloading the latest version');
            end
            pause(.5)
            %% Download lasted version
            filename = strcat('BCV_last_version.zip');
            disp(strcat("-->> Downloading latest version......."));
            url = local.generals.base_url;
            matlab.net.http.HTTPOptions.VerifyServerName = false;
            options = weboptions('Timeout',Inf,'RequestMethod','auto');
            websave(filename,url,options);

            if(getGlobalGuimode())
                d.Value = .33;
                d.Message = 'Unpacking files...';
            end
            pause(1)
            %% Unzip lasted version            
            disp(strcat("-->> Unpacking files..."));
            unzip(filename,pwd);
            pause(1);
            delete(filename);
            % Backuping configuration files
            if(getGlobalGuimode())
                d.Value = .67;
                d.Message = 'Deploying packages';
            end
            pause(1)
            mkdir('tmp');
            copyfile('bcv_properties/general_params.json','tmp/');
            copyfile('bcv_properties/sensor_params.json','tmp/');
            copyfile('bcv_properties/activation_params.json','tmp/');
            copyfile('bcv_properties/connectivity_params.json','tmp/');
            movefile( strcat('BC-VARETA_Toolbox-master',filesep,'*'), pwd,'f');
            rmdir BC-VARETA_Toolbox-master ;

            % Restoring configuration files
            movefile('tmp/general_params.json','bcv_properties/','f');
            movefile('tmp/sensor_params.json','bcv_properties/','f');
            movefile('tmp/activation_params.json','bcv_properties/','f');
            movefile('tmp/connectivity_params.json','bcv_properties/','f');
            rmdir('tmp');
            if(getGlobalGuimode())
                d.Value = 1;
                d.Message = 'Finishing';
            end
            disp('-->> The project is already update with the latest version.');
            disp('-->> The process was stoped to refresh all file');
            disp('-->> Please configure the app properties file, before restart the process.');
            pause(5);
            if(getGlobalGuimode())
                % Close dialog box
                close(d);
            end
            exit;
        else
            disp('-->> Nothing to update');
            disp('-->> You have already the BC-VARETA latest version')
            disp('=================================================================');            
        end
        status = true;
    else
        status = false;
    end
catch
    status = false;
    return;
end

end

