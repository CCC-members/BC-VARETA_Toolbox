function [status,message] = doc2pdf(docFile,pdfFile)
    %% DOC2PDF. Save Microsoft Word file into PDF programatically.
    %
    %  [status,message] = doc2pdf(docFile,pdfFile);
    %  Any file that can be opened with MS Word (doc, rtf, txt...) can be saved into PDF.
    %
    %  FILE NAMES MUST BE FULL PATH NAMES i.e.:  C:\Users\Desktop\Daniel\example.doc
    %  
    %  status is an output variable indicating successful conversion (true) or not (false).
    %  message appears if status is false, indicating the cause of the error.
    %
    %  docFile is the Microsoft Word file to be converted, specifying full path
    %  pdfFile is the PDF file into which doc should be converted, specifying full path
    %
    %  Example:
    %
    %   [OK,msg] = doc2pdf([pwd filesep 'myWord.docx'],fullfile(pwd,'PDFs/myWordConverted.pdf'));
    %
    %
    % To convert form Excel to PDF, use the xls2pdf tool:
    %   http://www.mathworks.com/matlabcentral/fileexchange/47308-doc2pdf-m
    %
    % To generate full paths to files, you can use GetFullPath from Matlab
    % Central File Exchange: http://www.mathworks.com/matlabcentral/fileexchange/28249-getfullpath
    
    status  = false;                              % Initialize status
    message = '';                                 % Initialize message
    
    try
        Word    = actxserver('Word.Application'); % Open MS Word
        Docu    = Word.Documents.Open(docFile);   % Open docFile with Word
        Docu.SaveAs2(pdfFile,17);                 % Save as (17 is the code for wdFormatPDF, 18 is for wdExportFormatXPS) specifying the name of the new file.
        % http://social.msdn.microsoft.com/Forums/en-US/d525b79d-ab45-4173-98d9-17f916c35ed4/save-word-document-as-pdf-using-vba?forum=isvvba
        % http://msdn.microsoft.com/en-us/library/bb243311%28v=office.12%29.aspx

        [~, ~, ext] = fileparts(docFile);         % Get file extension
        Docu.Close(false);                        % Close the document without saving
        Word.Quit;                                % Quit the application
        if strcmpi(ext,'.docx')
            delete(Word);                         % Delete the matlab object
        else
            Word.delete;
        end
        status  = true;                           % When all actions are completed, set the status to true
    catch exc
        message = sprintf('%s\n%s',exc.identifier,exc.message);% If error happens, show the message output.
    end