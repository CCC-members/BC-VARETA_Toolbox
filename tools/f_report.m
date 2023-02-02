function report = f_report(action,varargin)

%% Creating report function
%
%
%  Usage:    >> report = f_report(action, value1, value2, ... );
%       action='New'            -> Create a new report
%               report = f_report('New');
%
%       action='clean'            -> Clean the history reports folder
%               report = f_report('clean');
%
%       action='Header'          -> Add a header to the report
%               report = f_report('Header',title,text);
% 
%       action='Title'          -> Add a title to the report
%               report = f_report('Title',title);
% 
%       action='Index'          -> Add a topic index to the report
%               report = f_report('Index',text);
% 
%       action='Sub-Index'      -> Add a sub topic index to the report
%               report = f_report('Sub-Index',text);
% 
%       action='Info'           -> Add a info to the report
%               report = f_report('Info',information); information='some text' OR cell array {'info 1','info 2', ..... , 'info n'}
% 
%       action='Snapshot'      -> Add a figure to the report
%               report = f_report('Snapshot', figures); : figures (fig OR {fig1, fig2, .... , figN})
%               report = f_report('Snapshot', figures, header_text)
%               report = f_report('Snapshot', figures, header_text, description)
%               report = f_report('Snapshot', figures, header_text, description, img_dimention=[left_position bottom_position width high] )
%
%       action='Block'          -> Add a text description to the report
%               report = f_report('Block',title, text);
%
%       action='Table'          -> Add a table to the report
%               report = f_report('Table', table); (Table, Matrix OR Cell array)
%               report = f_report('Table', table, title); table=(Table Or Matrix), header=cell array
%               report = f_report('Table', table, title, colhead); table=(Table Or Matrix), header=cell array
%               report = f_report('Table', table, title, colhead, rowhead); table=(Table Or Matrix), header=cell array
%               report = f_report('Table', table, title, colhead, rowhead, description); table=(Table, Matrix, Cell array, Struct array), header=cell array
%
%       action='Footer'         -> Add a header to the report
%               report = f_report('Footer', title, text, key1, value1,........ , key3, value3);
%               report = f_report('Footer', title, text);
%               report = f_report('Footer', title, text, 'ref', references);
%               report = f_report('Footer', title, text, 'copyright', copyright);
%               report = f_report('Footer', title, text, 'contactus', contact);
%               report = f_report('Footer', title, text, 'ref', references, 'copyright', copyright);
%               report = f_report('Footer', title, text, 'ref', references, 'copyright', copyright, 'contactus', contact);
% 
%       action='Ref'            -> Add references to the report
%               report = f_report('Ref',references);
%
%       action='Export'         -> Export the current report  
%               html = f_report('Export')
%               html = f_report('Export', report_name)
%               html = f_report('Export', report_name, FileFormat)
%
%
% Authors:
% - Ariosky Areces Gonzalez
% - Deirel Paz Linares

% Date: Nov-2020

%%
tmp_path = fullfile(pwd,'tmp');
if(~isfolder(tmp_path))
    mkdir(tmp_path);
end
current_report = [];
reports = dir(tmp_path);
reports(ismember( {reports.name}, {'.', '..'})) = [];  %remove . and ..
for i=1:length(reports)
    try
        load(fullfile(reports(i).folder,reports(i).name));
        if(exist('report','var') && report.iscurrent)
            current_report = report;
            break;
        else
            delete(fullfile(reports(i).folder,reports(i).name));
        end
    catch
        continue;
    end
end
switch lower(action)
    case 'new'
        if(~isfolder(tmp_path))
            mkdir(tmp_path);
        end
        if(~isempty(current_report))
            current_report.iscurrent = false;
            report = current_report;
            save(fullfile(tmp_path,strcat(report.name,'.mat')),'report');
        end
        report = empty_report(tmp_path);
    case 'clean'
        if(isfolder(tmp_path))
            reports = dir(tmp_path);
            reports(ismember( {reports.name}, {'.', '..'})) = [];  %remove . and ..
            for i=1:length(reports)
                report = reports(i);
                delete(fullfile(report.folder,report.name));
            end
        end
    case 'title'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Title(current_report,tmp_path,varargin);
    case 'header'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Header(current_report,tmp_path,varargin);
    case 'info'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Info(current_report,tmp_path,varargin);
    case 'snapshot'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Snapshot(current_report,tmp_path,varargin);
    case 'export'
        if(isempty(current_report))
            fprintf(2,"\n ->> Error: there is not any report to export. \n");
            return;
        end
        report = Export(current_report,varargin);
    case 'index'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Index(current_report,tmp_path,varargin);
    case 'sub-index'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Sub_Index(current_report,tmp_path,varargin);
    case 'block'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Block(current_report,tmp_path,varargin);
    case 'table'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Table(current_report,tmp_path,varargin);
    case 'ref'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Ref(current_report,tmp_path,varargin);
    case 'footer'
        if(isempty(current_report))
            current_report = empty_report(tmp_path);
        end
        report = Footer(current_report,tmp_path,varargin);
    otherwise
        fprintf(2,"\n ->> Error: The action selected is not correct. \n");
        return;
end
end

%%
%% Creating a empty report
%%
function report = empty_report(base_path)
report              = struct;
report.name         = strcat('Report_',num2str(fix(posixtime(datetime()))));
report.iscurrent    = true;
report.time         = datetime();
metadata(1)         = struct('key', 'type', 'value', 'header');
metadata(2)         = struct('key', 'title', 'value', strcat('Report_',num2str(fix(posixtime(datetime())))));
metadata(3)         = struct('key', 'time', 'value', datestr(datetime()));
report.header       = struct('action', 'info', 'value', '', 'metadata', metadata);
report.body         = struct('action', 'info', 'value', '', 'metadata', []);
report.footer       = struct('action', 'footer', 'value', '', 'metadata', []);
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add a title to the current report
%%
% USAGE:  report = Title(current_report,base_path,title)
function report = Title(report,base_path,varargin)
varargin = varargin{1};
if(~isempty(varargin))
    metadata = report.header.metadata;
    metadata(2) = struct('key', 'title', 'value', strcat(convertStringsToChars(varargin{1})));
    report.header.metadata = metadata;
    save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end
end

%%
%% Add a header to the current report
%%
% USAGE:  report = Header(current_report,base_path,title,text)
function report = Header(report,base_path,varargin)
varargin = varargin{1};
metadata = report.header.metadata;
metadata(1) = struct('key', 'type', 'value', 'header');
if(~isempty(varargin))
    title = varargin{1};
    metadata(2) = struct('key', 'title', 'value', title);
end
if(length(varargin)>1)
    text  = varargin{2}; 
    report.header = struct('action', 'header', 'value', convertStringsToChars(text), 'metadata', metadata);
end

save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add a index to the current report
%%
% USAGE:  report = Index(current_report,base_path,title)
function report = Index(report,base_path,varargin)
varargin = varargin{1};
report.body(end+1) = struct('action', 'index', 'value', convertStringsToChars(varargin{1}), 'metadata', struct('key','type','value','text'));
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add a sub-index to the current report
%%
% USAGE:  report = Sub_Index(current_report,base_path,title)
function report = Sub_Index(report,base_path,varargin)
varargin = varargin{1};
report.body(end+1) = struct('action', 'sub-index', 'value', convertStringsToChars(varargin{1}), 'metadata', struct('key','type','value','text'));
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add a Info to the current report
%%
% USAGE:  report = Info(current_report, base_path, text)
function report = Info(report,base_path,varargin)
varargin = varargin{1};
if(iscell(varargin{1}))
    lines = varargin{1};
    for k=1:length(lines)
        report.body(end+1) = struct('action', 'info', 'value', convertStringsToChars(lines{k}), 'metadata', struct('key','type','value','text'));
    end
else
    report.body(end+1) = struct('action', 'info', 'value', convertStringsToChars(varargin{1}), 'metadata', struct('key','type','value','text'));
end
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add a text block to the current report
%%
% USAGE:  report = Block(current_report,base_path,title,text)
function report = Block(report,base_path,varargin)
varargin = varargin{1};
title = varargin{1};
text  = varargin{2};
metadata(1) = struct('key', 'type', 'value', 'block');
metadata(2) = struct('key', 'title', 'value', title);
report.body(end+1) = struct('action', 'block', 'value', convertStringsToChars(text), 'metadata', metadata);
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add a Snapshot to the current report
%%
% USAGE:    report = Snapshot(current_report, base_path, figures): figures (fig OR {fig1, fig2, .... , figN})
%           report = Snapshot(current_report, base_path, figures, header_text)
%           report = Snapshot(current_report, base_path, figures, header_text, description)
%           report = Snapshot(current_report, base_path, figures, header_text, description, img_dimention=[left_position bottom_position width high] )
function report = Snapshot(report, base_path, varargin)
varargin = varargin{1};
winPos = [200 200 875 550];
if(length(varargin)>3)
    winPos = varargin{4};
    if(~isempty(winPos))
        metadata(3) = struct('key', 'dim', 'value', winPos);
    end
end
if(length(varargin)>2)
    description = varargin{3};
    if(~isempty(description))
        metadata(2) = struct('key', 'description', 'value', convertStringsToChars(description));
    end
end
if(~isempty(varargin))
    metadata(1) = struct('key', 'type', 'value', 'image');
    figures = varargin{1};
    if(iscell(figures))
        for i=1:length(figures)
            fig = figures{i};
            report = Snapshot(report, base_path, fig, [], winPos);
        end
    else
        img = get_img(figures, winPos);
    end
end
if(length(varargin)>1)
    metadata(4) = struct('key', 'text', 'value', convertStringsToChars(varargin{2}));
end
report.body(end+1) = struct('action', 'snapshot', 'value', img,  'metadata', metadata);
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add a Table to the current report
%%
% USAGE:    report = Table(current_report, base_path, table) 
%           report = Table(current_report, base_path, table, title)
%           report = Table(current_report, base_path, table, title, colhead)
%           report = Table(current_report, base_path, table, title, colhead, rowhead)
%           report = Table(current_report, base_path, table, title, colhead, rowhead, description)
function report = Table(report, base_path, varargin)
varargin = varargin{1};
if(~isempty(varargin))
    table = varargin{1};
    if(isstruct(table))
        table = struct2table(table);
    end
    if(istable(table))
        table = table2cell(table);
    end
    if(ismatrix(table))
        table = num2cell(table);
    end
    metadata(1) = struct('key', 'type', 'value', 'table');
end
if(length(varargin)>1)
    table_title = varargin{2};
    metadata(2) = struct('key', 'title', 'value', table_title);
end
if(length(varargin)>2)
    header = varargin{3};
    if(~isempty(header))
        metadata(3) = struct('key', 'colhead', 'value', {header});
    end
end
if(length(varargin)>3)
    header = varargin{4};
    if(~isempty(header))
        metadata(4) = struct('key', 'rowhead', 'value', {header});
    end
end
if(length(varargin)>4)
    description = varargin{5};
    if(~isempty(description))
        metadata(5) = struct('key', 'description', 'value', description);
    end
end
report.body(end+1) = struct('action', 'table', 'value', {table},  'metadata', metadata);
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add the references to the current report
%%
% USAGE:  report = Ref(current_report,base_path,references): references = cell array
function report = Ref(report,base_path,varargin)
varargin = varargin{1};
references = varargin{1};
metadata(1) = struct('key', 'type', 'value', 'reference');
report.body(end+1) = struct('action', 'ref', 'value', {references}, 'metadata', metadata);
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Add a footer to the current report
%%
% USAGE:    report = Footer(current_report, base_path, title, text, key1, value1,........ , key3, value3);
%           report = Footer(current_report, base_path, title, text);
%           report = Footer(current_report, base_path, title, text, 'ref', references);
%           report = Footer(current_report, base_path, title, text, 'copyright', copyright);
%           report = Footer(current_report, base_path, title, text, 'contactus', contact);
%           report = Footer(current_report, base_path, title, text, 'ref', references, 'copyright', copyright);
%           report = Footer(current_report, base_path, title, text, 'ref', references, 'copyright', copyright, 'contactus', contact);
function report = Footer(report,base_path,varargin)
varargin = varargin{1};
title = varargin{1};
text  = varargin{2};
metadata(1) = struct('key', 'type', 'value', 'footer');
metadata(2) = struct('key', 'title', 'value', title);


ref_ind = find(strcmp(varargin,'ref'),1);
if(~isempty(ref_ind))
    metadata(end+1) = struct('key', 'ref', 'value', varargin(ref_ind+1));
end
copy_ind = find(strcmp(varargin,'copyright'),1);
if(~isempty(copy_ind))
    metadata(end+1) = struct('key', 'copyright', 'value', varargin{copy_ind+1});
end
contactus_ind = find(strcmp(varargin,'contactus'),1);
if(~isempty(contactus_ind))
    metadata(end+1) = struct('key', 'contactus', 'value', varargin{contactus_ind+1});
end

report.footer(end) = struct('action', 'footer', 'value', convertStringsToChars(text), 'metadata', metadata);
save(fullfile(base_path,strcat(report.name,'.mat')),'report');
end

%%
%% Export the current report to the spesific format
%%
% USAGE:    html = Export(report,varargin)
%           html = Export(report)
%           html = Export(report, report_name)
%           html = Export(report, report_name, FileFormat)
function html = Export(report,varargin)
varargin = varargin{1};
if(~isempty(varargin))
    output_path = varargin{1};
    report_name = report.name;
end
if(length(varargin)>1)
    report_name = varargin{2};
end
% Get output format (HTML or PNG)

FileFormat = 'HTML';
if(length(varargin)>2)
    FileFormat = varargin{3};
end

switch lower(FileFormat)
    case 'html'
        file_name = fullfile(output_path, strcat(report_name, '.html'));
        html = GetHtml(report);
        save(file_name);
        fid = fopen(file_name,'w');
        fprintf(fid,html);
        fclose(fid);
    case 'pdf'
end
end

function img = get_img(hFig, winPos)
drawnow;
% Set figure size
if ~isempty(winPos)
    set(hFig, 'Position', winPos);
end
% Redraw figure
figure(hFig);
drawnow;
% Capture window contents
% Matlab function getframe() was finally fixed in R2014b
frameGfx = getframe(hFig);
img = frameGfx.cdata;
end


%% ===== Get HTML File =====
function html = GetHtml(report)
html = '';

% If no errors, no warnings: no output
if isempty(report.body)
    disp("-->> Nothing to export in the current report.");
    return;
end
% ===== HEADER ====

% HTML header
html = ['<HTML>' 10 ...
    '<STYLE type="text/css">' 10 ...
    'h2 {font-size: 14px; font-weight: bold; padding-top: 12px; padding-bottom: 6px;}' 10 ...
    'li {font-size: 14px;}' 10 ...
    'ul {font-size: 14px;}' 10 ...
    'td {padding: 2px 2px 2px 2px; }' 10 ...
    '.bord td { border-width: 1px; border-style: solid; border-color: #bbbbbb; background-color: #f2f2f2; font: normal 8px Verdana, Arial, Helvetica, sans-serif; }' 10 ...
    '.link {text-decoration: none; color: #0000a0;}' 10 ...
    'p {}'...
    'table, th, td {border: 0.5px solid black; border-collapse: collapse; border-spacing: 0; padding: 10px; }'...
    'table {margin-left: auto; margin-right: auto;  border-radius: 0.25em; font: normal 14px Verdana, Arial, Helvetica;} '...
    'img {margin-left: auto; margin-right: auto;} '...
    '</STYLE>' 10 10 ...
    '<BODY style="padding:15px; margin: 5px 10px 10px 10px; font: normal 14px Verdana, Arial, Helvetica, sans-serif; background-color: #e8e8e8; width: 900px;">' 10 ...
    '<TITLE style="font: normal 14px Verdana, Arial, Helvetica;">' report.header.metadata(2).value '</TITLE>' 10];
html = [html '<div >' 10];
% Time

data_header = report.header;
metadata = data_header.metadata;
text = data_header.value;
html = [html '<H4 style="text-align:right;"> Date: ' metadata(3).value ' </H4><hr>' 10];
html = [html '<header style="text-align:center;"><H1>' metadata(2).value '</H1>'];
html = [html '<div>'];
if(~isempty(text))
   html = [html '<p>' text '</p>']; 
end
html = [html '</div></header>']; 

% ===== BODY =====
for i=1:length(report.body)
    data_row = report.body(i);
    
    switch lower(data_row.action)
        case 'info'           
            if(iscell(data_row.value))
                 html = [html '<div>'];
                for j=1:length(data_row.value)
                    html = [html '<p> ' convertStringsToChars(data_row.value(j)) ' </p>' 10];
                end
                 html = [html '</div>'];
            else
            html = [html '<div><p> ' convertStringsToChars(data_row.value) ' </p></div>' 10];
            end
        case 'snapshot'
%             html = [html '<div><H2> Snapshots </H2><hr>' 10];
            imgRgb = data_row.value;
            metadata = data_row.metadata;
            % Create Base64 encoder
            encoder = sun.misc.BASE64Encoder();
            row = find(strcmp({metadata.key},'text'),1);
            Comment   = metadata(row).value;
            % Convert RGB in 255 to Java color integer
            sz = size( imgRgb );
            imgRgb = double(imgRgb) / 255;
            imgRgb = transpose( reshape( permute( imgRgb, [3 2 1] ), [sz(3) sz(1)*sz(2)] ) );
            imgInt = typecast(bitshift( uint32( 255 ), 24 ), 'int32') + ...
                typecast(bitshift( uint32( 255*imgRgb(:,1) ), 16 ), 'int32') + ...
                typecast(bitshift( uint32( 255*imgRgb(:,2) ), 8 ), 'int32') + ...
                typecast(bitshift( uint32( 255*imgRgb(:,3) ), 0 ), 'int32');
            % Create an image (BufferedImage)
            jImage = java.awt.image.BufferedImage(sz(2), sz(1), java.awt.image.BufferedImage.TYPE_INT_ARGB);
            jImage.setRGB(0, 0, sz(2), sz(1), imgInt(:), 0, sz(2));
            % Convert image to PNG
            jByteStream = java.io.ByteArrayOutputStream();
            javax.imageio.ImageIO.write(jImage, 'gif', jByteStream);
            % Encode PNG image in Base64
            jStringImage = encoder.encode(jByteStream.toByteArray());
            % Display image in HTML
            html = [html, '<H3>', Comment, '</H3><hr>'];
            html = [html, '<BR>'];
            html = [html, '<IMG src="data:image/gif;base64,' char(jStringImage) '" /><BR><BR>'];
            row = find(strcmp({metadata.key},'description'),1);
            if(~isempty(row))
                description   = metadata(row).value;
                html = [html, '<p>', description, '</p>'];
            end
                html = [html, '</div>'];
            
        case 'index'
             html = [html '<ul><li><H1> ' data_row.value ' </H1></li></ul>' 10];
        case 'sub-index'
            html = [html '<ul><li><H2> &#9;' data_row.value ' </H2></li></ul>' 10];
        case 'block'
            metadata = data_row.metadata;            
            title_ind = find(strcmp({metadata.key},'title'),1);
            title = metadata(title_ind).value;
            html = [html '<div> <hr><H2> ' convertStringsToChars(title) ' </H2>'];
            html = [html  data_row.value ' </div>' 10];
        case 'table'
            metadata = data_row.metadata;
            data = data_row.value;
            title_ind = find(strcmp({metadata.key},'title'),1);
            if(~isempty(title_ind))
                title = metadata(title_ind).value;
                html = [html '<div><hr><H2> ' convertStringsToChars(title) ' </H2>' 10];
            end
            
            html = [html ' <table> '];
            colhead_ind = find(strcmp({metadata.key},'colhead'),1);
            rowhead_ind = find(strcmp({metadata.key},'rowhead'),1);
            if(~isempty(colhead_ind))
                html = [html '<thead> <tr> '];
                if(~isempty(rowhead_ind))
                    html = [html '<th>' '</th>'];
                    rowheads = metadata(rowhead_ind).value;
                end
                colheads = metadata(colhead_ind).value;
                for m=1:length(colheads)
                    colhead = colheads{m};
                    if(isnumeric(colhead))
                        colhead = str2num(colhead);
                    end
                    html = [html '<th>' convertStringsToChars(colhead) '</th>'];
                end
                html = [html '</tr></thead>'];
            end
            html = [html '<tbody><tr> '];   
            
            for m=1:size(data,1)
                html = [html '<tr> '];
                if(~isempty(rowhead_ind))
                    rowhead = rowheads{m};
                    if(isnumeric(rowhead))
                        rowhead = str2num(rowhead);
                    end
                    html = [html '<th>' convertStringsToChars(rowhead) '</th>'];
                end
                for h=1:size(data,2) 
                       element = data{m,h}{:};
                    if(isnumeric(element))
                        element = num2str(element);
                    end
                    if(isstruct(element))
                       element = 'struct'; 
                    end
                    if(iscell(element))
                       element = 'cell'; 
                    end
                    html = [html '<td>' convertStringsToChars(element) '</td>'];
                end
                 html = [html '</tr>'];
            end
                html = [html '</tr></tbody>'];
            html = [html '</table> <BR><BR>'];
            row = find(strcmp({metadata.key},'description'),1);
            if(~isempty(row))
                description   = metadata(row).value;
                html = [html, '<p>' convertStringsToChars(description) '</p>'];
            end
            html = [html '</div> '];   
        case 'ref'
            html = [html '<div><hr> '];
            html = [html '<H2> References </H2>'];
            html = [html '<ul>'];
            references = data_row.value;
            for h=1:length(references)
                ref = references{h};
                html = [html '<li><a href="' convertStringsToChars(ref) '">' convertStringsToChars(ref) '</a></li>'];
            end
            html = [html  '</ul>' 10];
            html = [html '</div> '];  
    end    
end
% HTML footer
html = [html '</div><div><footer style="text-align: center;"> <hr>'];

data_footer = report.footer;
metadata = data_footer.metadata;
text = data_footer.value;
if(~isempty(metadata))
    title_ind = find(strcmp({metadata.key},'title'),1);
    ref_ind = find(strcmp({metadata.key},'ref'),1);
    copy_ind = find(strcmp({metadata.key},'copyright'),1);
    contactus_ind = find(strcmp({metadata.key},'contactus'),1);
    if(~isempty(title_ind))
        title = metadata(title_ind).value;
        html = [html '<H2>' convertStringsToChars(title) '</H2>'];
    end
    if(~isempty(text))
        html = [html '<p>' convertStringsToChars(text) '</p>'];
    end
    if(~isempty(ref_ind))
        html = [html '<H3> References </H3>'];
        references = metadata(ref_ind).value;
        for h=1:length(references)
            ref = references{h};
            html = [html '<a href="' ref '">' convertStringsToChars(ref) '</a> | '];
        end
    end
    if(~isempty(contactus_ind))
        contact = metadata(contactus_ind).value;
        html = [html '<BR><a href = "mailto: ' convertStringsToChars(contact) '">Contact us</a>'];
    end
    if(~isempty(copy_ind))
        copy   = metadata(copy_ind).value;
        html = [html '<p>Copyright &copy; ' convertStringsToChars(copy) '</p>'];
    end
end

html = [html '</footer></div></BODY></HTML>'];
% HTML footer

end


