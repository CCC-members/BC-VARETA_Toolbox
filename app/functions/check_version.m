function check_version()

%  BC-VARETA check version
%
%
% Authors:
%   -   Ariosky Areces Gonzalez
%   -   Deirel Paz Linares
%   -   Eduardo Gonzalez Moreaira
%   -   Pedro Valdes Sosa

% - Date: November 15, 2019
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
                answer = questdlg({'There a new version available of BC-VARETA Toolbox.', ...
                    strcat('Version:',online.generals.version), ...
                    strcat('Date:',online.generals.version_date), ...
                    ' Do you want to update the laster version?'}, ...
                    'Update', ...
                    'Yes','No','Close');
                % Handle response
                switch answer
                    case 'Yes'
                        f = dialog('Position',[300 300 250 80]);

                        iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
                        iconsSizeEnums = javaMethod('values',iconsClassName);
                        SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
                        jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Starting update.');  % icon, label

                        jObj.setPaintsWhenStopped(true);  % default = false
                        jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
                        javacomponent(jObj.getComponent, [50,10,150,80], f);
                        jObj.start;
                        pause(1);


                        %% Download lasted version
                        filename = strcat('BCV_last_version.zip');
                        disp(strcat("-->> Downloading latest version......."));
                        jObj.setBusyText(strcat("Downloading latest version "));

                        url = local.generals.base_url;
                        matlab.net.http.HTTPOptions.VerifyServerName = false;
                        options = weboptions('Timeout',Inf,'RequestMethod','auto');
                        websave(filename,url,options);

                        %% Unzip lasted version
                        pause(1);
                        jObj.setBusyText('Unpacking version...');
                        disp(strcat("-->> Unpacking version..."));

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

                        jObj.stop;
                        jObj.setBusyText('All done!');
                        disp(strcat("-->> All done!"));
                        pause(2);
                        delete(f);

                        disp('-->> The project is already update with the latest version.');
                        disp('-->> The process was stoped to refresh all file');
                        disp('-->> Please configure the app properties file, before restart the process.');

                        uiwait(msgbox(["The latest version was downloaded successfully!";...
                            "Matlab will be closed to complete the installation."],"Success","modal"));
                        exit;

                    case 'No'
                        return;
                    case ''
                        return;   
                end   
            else
                disp('=================================================================');
                cprintf('_red','BC-V -->> There is a BC-VARETA new version available.');fprintf('\n');
                disp(strcat("-->> Name:",online.generals.name));
                disp(strcat("-->> Version:",online.generals.version));
                disp(strcat("-->> Version date:",online.generals.version_date));
                disp("-----------------------------------------------------------------");
                disp("-->> Download the latest version running:");
                disp("-->> Main update");
                disp('=================================================================');
            end
        else
            disp('-->> Nothing to update');
            disp('=================================================================');
        end
    end
catch
    return;
end
end