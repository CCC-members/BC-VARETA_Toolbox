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
                        update_version();
                    case 'No'
                        return;
                    case ''
                        return;   
                end   
            else
                disp('=================================================================');
                cprintf('Keywords','There is a BC-VARETA new version available.');fprintf('\n');
                disp("-----------------------------------------------------------------");
                cprintf('Keywords',strcat("-->> Name:",online.generals.name));fprintf('\n');
                cprintf('Keywords',strcat("-->> Version:",online.generals.version));fprintf('\n');
                cprintf('Keywords',strcat("-->> Version date:",online.generals.version_date));fprintf('\n');
                disp("-----------------------------------------------------------------");
                cprintf('Keywords',"-->> Download the latest version running: "); 
                cprintf('Text',"bcvareta ");cprintf('SystemCommands',"update");fprintf('\n');
                cprintf('Keywords',"     OR");fprintf('\n');
                cprintf('Keywords',"-->> Download the latest version directly from: ");
                cprintf('Keywords','<a href="https://github.com/CCC-members/BC-VARETA_Toolbox">BC-VARETA Toolbox Github</a>');fprintf('\n');
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