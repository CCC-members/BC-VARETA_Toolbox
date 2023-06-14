function Main(varargin)
%% BC-VARETA toolbox v1.0
%%%%%%%%%%%%%%%%%%%%

% Includes the routines of the Brain Connectivity Variable Resolution
% Tomographic Analysis (BC-VARETA), an example for real EEG analysis.
% BC-VARETA toolbox extracts the Source Activity and Connectivity given
% a single frequency component in the Fourier Transform Domain of an
% Individual MEEG Data. See the pdf file "Brief of Theory and Results"
% for an insight to this methodology.

% Authors:
% - Ariosky Areces Gonzalez
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa



%% Preparing WorkSpace
clc;
close all;
restoredefaultpath;
clearvars -except varargin;
disp('-->> Starting process');
disp("=====================================================================");
%restoredefaultpath;
tic
addpath('app');
addpath('bcv_properties');
addpath(genpath('external'));
addpath(genpath('functions'));
addpath(genpath('guide'));
addpath('tools');
% Remove fieldtrip path for override functions 
warning off;
rmpath(genpath(fullfile('external/fieldtrip')));
warning on;


if(isequal(nargin,2))
    idnode = varargin{1};
    count_node = varargin{2};
    if(~isnumeric(idnode) || ~isnumeric(count_node))
        fprintf(2,"\n ->> Error: The selected node and count of nodes have to be numbers \n");
        return;
    end
else
    idnode = 1;
    count_node = 1;
end

setGlobalGuimode(true);
for i=1:length(varargin)
    if(isequal(varargin{i},'nogui'))
      setGlobalGuimode(false);
    end
end

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
    disp('-->> Checking last project version');
    if(isequal(check_version,'updated'))
        return;
    end
end
%%               Upload the actived processes
if(getGlobalGuimode())
    BC_VARETA
else
    BC_VARETA_bash(idnode,count_node);
end

%%
% delete(process_waitbar);
disp('=====================================================================');
disp(strcat("BC-V-->> Process completed on instance: ",num2str(idnode),"."));
hours = fix(toc/3600);
minutes = fix(mod(toc,3600)/60);
disp(strcat("Elapsed time: ", num2str(hours) , " hours with ", num2str(minutes) , " minutes." ));
disp('=====================================================================');
disp(app_properties.generals.name);
disp('=====================================================================');
