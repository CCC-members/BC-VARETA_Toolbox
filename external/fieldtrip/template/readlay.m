function layout = readlay(filename)
if ~exist(filename, 'file')
    ft_error('could not open layout file "%s"', filename);
end
fid=fopen(filename);
lay_string=fread(fid,inf,'char=>char')';
fclose(fid);

% pattern to match is integer, than 4 numeric values followed by a
% string that can contain whitespaces and plus characters, followed by
% newline
integer='(\d+)';
float='([\d\.-]+)';
space='\s+';
channel_label='([\w \t\r\f\v\-\+]+)';
single_newline='\n';

pat=[integer, space, ...
    float, space, ...
    float, space, ...
    float, space, ...
    float, space, ...
    channel_label, single_newline];

matches=regexp(sprintf('%s\n',lay_string),pat,'tokens');

% convert to (nchannel x 6) matrix
layout_matrix=cat(1,matches{:});

% convert values in first five columns to numeric
num_values_cell=layout_matrix(:,1:5)';
str_values=sprintf('%s %s %s %s %s; ', num_values_cell{:});
num_values=str2num(str_values);

% store layout information (omit channel number in first column)
layout.pos    = num_values(:,2:3);
layout.width  = num_values(:,4);
layout.height = num_values(:,5);

% trim whitespace around channel names
label=layout_matrix(:,6);
label=regexprep(label,'^\s*','');
label=regexprep(label,'\s*$','');
layout.label  = label;
end

