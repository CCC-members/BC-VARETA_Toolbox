function [] = saveJSON(data,output_file)
%SAVEJSON Summary of this function goes here
%   Detailed explanation goes here
 
data = jsonencode(data);
data = strrep(data, ',', sprintf(',\r\t'));
data = strrep(data, '[', sprintf('[\r\t'));
data = strrep(data, '{', sprintf('{\r\t'));
data = strrep(data, '}', sprintf('\r}'));
data = strrep(data, ']', sprintf('\r]'));
fid = fopen(output_file, 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, data, 'char');
fclose(fid);
end
