
% Exporting figures with transparent background
addpath(genpath('../tools'));
base_folder = "/mnt/Data/NIM_Data/Figures";
figures = dir(fullfile(base_folder,'*.fig'));
for i=1:length(figures)
    figure = figures(i);
    export_fig(fullfile(figure.folder,figure.name),'-transparent','-png');
end