function FillingDocHoles(template_file, report_file, reportData, varargin)

import mlreportgen.dom.*

for i=1:length(varargin)
    eval([inputname(i+3) '= varargin{i};']);
end

if(~exist('nfreqs','var'))
    nfreqs = [];
end

%Create new object of class Document, based on AdvancedReportTemplate
doc = Document(report_file, 'docx', fullfile(template_file));
%Move to the first hole
currentHoleId = moveToNextHole(doc);
%While we have not reached the end, we loop through all the holes
while ~strcmp(currentHoleId, '#end#')
    %For some holes we perform a special handling, all "standard" holes
    %are filled with processStandardHole in the otherwise case
    try
        clasif = split(currentHoleId,'_');
        switch clasif{1}
            case 'Image'
                processImage(doc, reportData, currentHoleId, nfreqs);
            case 'SimpleTable'
                processSimpleTable(doc, reportData, currentHoleId);
            case 'AdvancedTable'
                processAdvancedTable(doc, reportData, currentHoleId);
            otherwise
                processStandardHole(doc, reportData, currentHoleId);
        end
        currentHoleId = moveToNextHole(doc);
    catch
        currentHoleId = moveToNextHole(doc);
    end
end
%Close the document and write the result to disc
close(doc);
%Show the result
if ismac
    % Code to run on Mac platform
elseif isunix
    sub_path = fileparts(report_file);
    %% Cerate libreoffice alias in linux
    system(strcat("libreoffice24.8 --headless --convert-to pdf ",report_file,".docx"," --outdir ",sub_path));
elseif ispc
    rptview(report_file, 'docx');
    docview(report_file,'convertdocxtopdf','closeapp');
    % [status,message] = doc2pdf(strcat(report_file,'.docx'),strcat(report_file,'.pdf'));
else
    disp('Platform not supported')
end
delete(strcat(report_file, '.docx'));

end

