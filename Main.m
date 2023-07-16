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
addpath(genpath('app'));
addpath(genpath('bcv_properties'));
addpath(genpath('external'));
addpath(genpath('functions'));
addpath(genpath('guide'));
addpath('tools');
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

%% Init processing
app_properties = init_processing();

%% BC-VARETA processing
if(getGlobalGuimode())
    BC_VARETA
else
    BC_VARETA_bash(idnode,count_node);
end

%%
disp('=====================================================================');
disp(strcat("BC-V-->> Process completed on instance: ",num2str(idnode),"."));
hours = fix(toc/3600);
minutes = fix(mod(toc,3600)/60);
disp(strcat("Elapsed time: ", num2str(hours) , " hours with ", num2str(minutes) , " minutes." ));
disp('=====================================================================');
disp(app_properties.generals.name);
disp('=====================================================================');
