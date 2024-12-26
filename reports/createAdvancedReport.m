function outputFile = createAdvancedReport(Partic_info, reportPath, type, varargin)
%ADVANCEDREPORTCREATE This function creates the file AdvancedReport.docx
%   First it calls the function "getAdvancedReportData" to get the data
%   that should be shown in the report. Then it loops through the holes in
%   "AdvancedReportTemplate" and fills the holes with the data.
%   Copyright 2016 - 2016 The MathWorks, Inc.

for k = 4 : nargin
    eval([inputname(k) '=  varargin{k-3};']);
end

part = 1;

% Filling presentation document
disp("-->> Processing presentation subdocument");
switch type
    case 'Subject'
        templates = dir(fullfile(pwd,'reports','templates','subject','*.dotx'));
    case 'Dataset'
        templates = dir(fullfile(pwd,'reports','templates','medium','*.dotx'));
    case 'Expert'
        templates = dir(fullfile(pwd,'reports','templates','expert','*.dotx'));
    otherwise
        templates = dir(fullfile(pwd,'reports','templates','*.dotx'));
end
template = templates(find(ismember({templates.name},'01_presentation.dotx'),1));
template_file = fullfile(template.folder,template.name);
report_file = getReportName(reportPath,Partic_info, part);
reportData = getAdvancedReportData('Presentation',Partic_info);
FillingDocHoles(template_file, report_file, reportData);
part = part + 1;

% Filling Structural document
if(isfield(prop,'general_params'))
    disp("-->> Processing Structural subdocument");
    template = templates(find(ismember({templates.name},'02_structural.dotx'),1));
    template_file = fullfile(template.folder,template.name);
    report_file = getReportName(reportPath,Partic_info, part);
    reportData = getAdvancedReportData('Structural',Partic_info, reportPath);
    FillingDocHoles(template_file, report_file, reportData);
    part = part + 1;
end


% Filling Sensor level Analisis document 
if(isfield(prop,'sensor_params'))
    disp("-->> Processing Sensor Level Analisis subdocument");
    template = templates(find(ismember({templates.name},'03_sensor_level_analysis_psd.dotx'),1));
    template_file = fullfile(template.folder,template.name);
    report_file = getReportName(reportPath,Partic_info, part);
    reportData = getAdvancedReportData('Sensor_PSD',Partic_info,reportPath);
    FillingDocHoles(template_file, report_file, reportData);
    part = part + 1;

    % Filling Sensor level Analisis Topography document
    if(prop.general_params.run_frequency_bin.value)
        template = templates(find(ismember({templates.name},'03_sensor_level_analysis_topo_bin.dotx'),1));
        template_file = fullfile(template.folder,template.name);
        report_file = getReportName(reportPath,Partic_info, part);
        reportData = getAdvancedReportData('Sensor_Topography_bin',Partic_info,reportPath);
        FillingDocHoles(template_file, report_file, reportData,prop);
        part = part + 1;
    else
        nfreqs = length(prop.sensor_params.frequencies);
        template = templates(find(ismember({templates.name},strcat('03_sensor_level_analysis_topo_',num2str(nfreqs),'_band.dotx')),1));
        template_file = fullfile(template.folder,template.name);
        report_file = getReportName(reportPath,Partic_info, part);
        reportData = getAdvancedReportData('Sensor_Topography_band',Partic_info,reportPath,prop);
        FillingDocHoles(template_file, report_file, reportData,prop,nfreqs);
        part = part + 1;
    end
end

% Filling Activation level Analisis document 
if(isfield(prop,'activation_params'))
    disp("-->> Processing Activity level subdocument");
    nfreqs = length(prop.sensor_params.frequencies);
    disp(strcat("---->> Processing Activity "));
    template = templates(find(ismember({templates.name},'04_source_level_analysis_topo_5_band.dotx'),1));
    template_file = fullfile(template.folder,template.name);
    report_file = getReportName(reportPath,Partic_info, part);
    reportData = getAdvancedReportData('Activation_Level',Partic_info,reportPath,prop);
    FillingDocHoles(template_file, report_file, reportData,prop,nfreqs);
    part = part + 1;
end

% Filling Connectivity level Analisis document 
if(isfield(prop,'connectivity_params'))

end    
%%
%% merging report parts
%%
disp("-->> Joinning report subdocument");
pdffiles = dir(fullfile(reportPath,'*.pdf'));
pdf_filenames = fullfile({pdffiles.folder},{pdffiles.name});
outputFileName = strcat(Partic_info.SubID, '_Report.pdf');
mergePdfs(pdf_filenames,  fullfile(reportPath,outputFileName));
if(~isempty(find(contains(pdf_filenames,fullfile(reportPath,outputFileName)),1)))
    pdf_filenames(find(contains(pdf_filenames,fullfile(reportPath,outputFileName)),1)) = [];
end
delete(pdf_filenames{:});
outputFile = fullfile(reportPath,outputFileName);
end

%%
%% Get report file
%%
function report_file = getReportName(reportPath,Partic_info, part)
if (part<10)
    report_file = fullfile(reportPath,strcat(Partic_info.SubID,'_report_part_0',num2str(part)));
else
    report_file = fullfile(reportPath,strcat(Partic_info.SubID,'_report_part_',num2str(part)));
end
end

%%
%% Get report file
%%
function sessionName = getSessionName(number)
if (number<10)
    sessionName = strcat('0',num2str(number));
else
    sessionName = strcat(num2str(number));
end
end


