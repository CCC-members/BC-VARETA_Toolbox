function parsave(file_name,varargin)
variables = cell(1,nargin-1);
for k = 2 : nargin
    eval([inputname(k) '=  varargin{k-1};']);
    variables(1,k-1) = {strcat(inputname(k))};
end
disp("-->> Saving file");
save(file_name,variables{:},'-v7.3','-nocompression');
end

