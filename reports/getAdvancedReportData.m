function data = getAdvancedReportData(section,varargin)
%ADVANCEDREPORTDATA Defines an Interface between the report creation
%   function and the data that shall be displayed in the report
%   The function fills a struct where every structure field corresponds
%   to a hole name in the Word template.
%   Copyright 2016 - 2016 The MathWorks, Inc.

%% Inputs
%  section, Partic_info, EEG
%  section, Partic_info, EEG, Session
%  section, Partic_info, EEG, reportPath,'figures', Session
%  section, Partic_info, EEG, reportPath,'figures', Session1, Session2

for i=1:length(varargin)
    eval([inputname(i+1) '= varargin{i};']);
end

%%
%% Defining import data
%%
switch section
    case 'Presentation'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;        
        data.Gender = Partic_info.Gender;
        data.Condition = Partic_info.Condition;
    case 'Structural'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;       
            data.Gender = Partic_info.Gender;            
        data.Condition = Partic_info.Condition;

        data.Image_headmodel_front = fullfile(reportPath,'structural','headmodel_front.png');
        data.Image_headmodel_top = fullfile(reportPath,'structural','headmodel_top.png');
        data.Image_headmodel_right = fullfile(reportPath,'structural','headmodel_right.png');
        data.Image_headmodel_left = fullfile(reportPath,'structural','headmodel_left.png');

        data.Image_sourcemodel_top = fullfile(reportPath,'structural','sourcemodel_top.png');
        data.Image_sourcemodel_bottom = fullfile(reportPath,'structural','sourcemodel_bottom.png');
        data.Image_sourcemodel_right = fullfile(reportPath,'structural','sourcemodel_right.png');
        data.Image_sourcemodel_left = fullfile(reportPath,'structural','sourcemodel_left.png');

        data.Image_leadfield_top = fullfile(reportPath,'structural','leadfield_top.png');
        data.Image_leadfield_front = fullfile(reportPath,'structural','leadfield_front.png');
        data.Image_leadfield_right = fullfile(reportPath,'structural','leadfield_right.png');
        data.Image_leadfield_left = fullfile(reportPath,'structural','leadfield_left.png');

        data.Image_channels_front = fullfile(reportPath,'structural','channels_front.png');
        data.Image_channels_back = fullfile(reportPath,'structural','channels_back.png');
        data.Image_channels_right = fullfile(reportPath,'structural','channels_right.png');
        data.Image_channels_left = fullfile(reportPath,'structural','channels_left.png');
    case 'Sensor_PSD'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;       
            data.Gender = Partic_info.Gender;           
        data.Condition = Partic_info.Condition;
        data.Image_sensorPSD = fullfile(reportPath,'sensor','sensorPSD.png');
    case 'Sensor_Topography_band'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
        data.Gender = Partic_info.Gender;
        data.Condition = Partic_info.Condition;
        for i=1:length(prop.sensor_params.frequencies)
            band = prop.sensor_params.frequencies(i);
            data.(strcat('name_band_',band.name)) = band.str_band;
            data.(strcat('Image_bandtopography_band_',band.name)) = fullfile(reportPath,'sensor',strcat('sensorTopography_',band.str_band,'.png'));
        end

    case 'Sensor_Topography_bin'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
       data.Gender = Partic_info.Gender;
       data.Condition = Partic_info.Condition;
        for i=1:length(prop.sensor_params.frequencies)
            data.Image_sensorPSD = fullfile(reportPath,'sensor','sensorPSD.png');
        end

    case 'Activation_Level'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
       data.Gender = Partic_info.Gender;
       data.Condition = Partic_info.Condition;
        data.view_right = 'Right view';
        data.view_top = 'Top view';
        data.view_left= 'Left view';
        for i=1:length(prop.sensor_params.frequencies)
            band = prop.sensor_params.frequencies(i);
            data.(strcat('band_name_',num2str(i))) = band.str_band;
            data.(strcat('Image_source_band_',band.name,'_right')) = fullfile(reportPath,'activation',strcat('sourceTopography_',band.str_band,'_right.png'));
            data.(strcat('Image_source_band_',band.name,'_top')) = fullfile(reportPath,'activation',strcat('sourceTopography_',band.str_band,'_top.png'));
            data.(strcat('Image_source_band_',band.name,'_left')) = fullfile(reportPath,'activation',strcat('sourceTopography_',band.str_band,'_left.png'));
        end
    case 'SLA_Band_Topography'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;       
            data.Gender = Partic_info.Gender;             
        data.Condition = Partic_info.Condition;

        data.SLA_band_topography = 'Band Topography';
        data.Session = strcat("Session ", Session);
        data.Image_bandtopography_band_delta = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-bandtopography_band-delta.png'));
        data.Image_bandtopography_band_theta = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-bandtopography_band-theta.png'));
        data.Image_bandtopography_band_alpha = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-bandtopography_band-alpha.png'));
        data.Image_bandtopography_band_beta = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-bandtopography_band-beta.png'));
    case 'AlphaThetaRation_1'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
       data.Gender = Partic_info.Gender;
        data.Condition = Partic_info.Condition;

        data.SLA_Alpha_Theta_Ratio = 'Alpha/Theta Ratio';
        data.Session = strcat("Session ", Session);
        data.Image_alphathetaration1 = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-alphathetaration.png'));
    case 'AlphaThetaRation_2'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
        data.Gender = Partic_info.Gender;
        data.Condition = Partic_info.Condition;

        data.SLA_Alpha_Theta_Ratio = 'Alpha/Theta Ratio';
        data.Session1 = strcat("Session ", Session1);
        data.Session2 = strcat("Session ", Session2);
        data.Image_alphathetaration1 = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session1,'_desc-alphathetaration.png'));
        data.Image_alphathetaration2 = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session2,'_desc-alphathetaration.png'));
    case 'SLA_AlphaPeakPosition'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
        data.Gender = Partic_info.Gender;       
        data.Condition = Partic_info.Condition;

        data.SLA_Alpha_Peak_Position = 'Alpha Peak Position';
        data.Session = strcat("Session ", Session);
        data.Image_alphapeakpositionPSD = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-alphapeakpositionPSD.png'));
        data.Image_alphapeakposition = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-alphapeakposition.png'));
    case 'SLA_Normative'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
        data.Gender = Partic_info.Gender;       
        data.Condition = Partic_info.Condition;

        data.SLA_Normative = 'Normative Distribution of Alpha Peak';
        data.Image_normativedist = fullfile(reportPath,'figures','sensor_level',strcat(Partic_info.SubID,'_desc-normativedist.png'));
    case 'SoLA'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
        data.Gender = Partic_info.Gender;       
        data.Condition = Partic_info.Condition;

    case 'SoLA_Activity_Connectivity'
        data.Name = Partic_info.Name;
        data.SubID = Partic_info.SubID;
        data.Age   = Partic_info.Age;
        data.Gender = Partic_info.Gender;       
        data.Condition = Partic_info.Condition;

        data.Session = strcat("Session ", Session);
        data.SoLA_Delta_band = 'Delta band';
        data.Image_deltacoronal = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-deltacoronal.png'));
        data.Image_deltalateral = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-deltalateral.png'));
        data.SoLA_Theta_band = 'Theta band';
        data.Image_thetacoronal = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-thetacoronal.png'));
        data.Image_thetalateral = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-thetalateral.png'));
        data.SoLA_Alpha_band = 'Alpha band';
        data.Image_alphacoronal = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-alphacoronal.png'));
        data.Image_alphalateral = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-alphalateral.png'));
        data.SoLA_Beta_band = 'Beta band';
        data.Image_betacoronal = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-betacoronal.png'));
        data.Image_betalateral = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-betalateral.png'));

        data.SoLA_Activity_Connectivity = 'Activity/Connectivity';
        data.Image_activityconnectivity_part_01 = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-activityconnectivity_part-01.png'));
        data.Image_activityconnectivity_part_02 = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-activityconnectivity_part-02.png'));
        data.Image_activityconnectivity_part_03 = fullfile(reportPath,'figures','source_level',strcat(Partic_info.SubID,'_sess-',Session,'_desc-activityconnectivity_part-03.png'));
end


% data.Image    = fullfile(pwd, 'Image.png');

% data.SimpleTable   = getTableData;
% data.AdvancedTable = getTableData;

% data.Text = evalc('why');

end

function table = getTableData
%This function uses the MATLAB function "why", which returns a different
%random sentence everytime it is called. This sentence is then transformed
%into a 3-by-n cell array, to represent our table data

numCol = 3;
%Why has no return value, so we need to capture the output with evalc
randomString = evalc('why');
%Remove the last carriage return
randomString(end) = [];
%Split the string at each space to get a cell array
cellString = strsplit(randomString, ' ');
le = length(cellString);
if mod(le, numCol) ~= 0
    %If the length is not a multiple of three, we need to add empty elements
    numElementsToAdd = numCol*(floor(le/numCol)+1) - le;
    cellString = [cellString repmat({''},1,numElementsToAdd)];
end
%Reshape the cell array into a 3-by-n cell array
table = reshape(cellString, [numCol, length(cellString)/3])';
end

