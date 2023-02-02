function fig_out = merge_figures(fig_name, fig_title, figures, varargin)

%% Merge figures function
%
%
% Authors
% - Ariosky Areces Gonzalez
% - Deirel Paz Linares
%
% Inputs
% - fig_name    :Output figure name
% - fig_tile    :Output figure title
% - figures     :Figures list in array or cell format
%   Optional input params
% - width       :Output width figure ('width',value)
% - height      :Output height figure ('height',value)
% - rows        :Rows subplots for the output figure
% - cols        :Columns subplots for the output figure
%
%
%

for i=1:2:length(varargin)
    eval([varargin{i} '=  varargin{(i+1)};'])
end
Position = {[0 0.5 0.5 0.5],[0.5 0.5 0.5 0.5],[0 0 0.5 0.5],[0.5 0 0.5 0.5]};
if(~exist('position','var'))
    position = 'absolute';
end
if(~exist('rows','var'))
    rows = 1;
end
if(~exist('cols','var'))
    cols = length(figures);
end
figures_path = cell(1,length(figures));
tmp_path = fullfile(pwd,'tmp','figures');
if(~isfolder(tmp_path))
   mkdir(tmp_path); 
end
for i=1:length(figures)
    if(iscell(figures))
        fig = figures{i};
    else
        fig = figures(i);
    end
    if(isempty(fig))
        fig_path = [];
    else
        fig_path = fullfile(tmp_path,['fig_',num2str(i),'.fig']);
        savefig(fig,fig_path);
    end
    figures_path{i} = fig_path;
end
if(exist('width','var') && exist('height','var'))
   fig_Position = [200,200,width,height];
else
    fig_Position = [200,200,900,700];   
end
fig_out = figure('Name', fig_name, 'NumberTitle', 'off', 'Position', fig_Position);
axis off;
for i=1:rows
    for j=1:cols
        p = (i-1)*cols+j;
        if(p>length(figures_path))
            break;
        end
        sub_p=subplot(rows,cols,p);   
        set(sub_p,'XTick',[], 'YTick', []); % all in one  
        if(~isequal(position,'relative'))
            sub_p.Position = Position{p};
        end
        if(isempty(figures_path{p}));continue;end
        f_c = openfig(figures_path{p});        
        clear title;        
        H = findobj(gca, 'Type','Text');
        if(~isempty(H))
            delete(H);
        end        
        % Identify axes to be copied
        axes_to_be_copied = findobj(f_c,'type','axes');
        % Identify the children of this axes
        chilred_to_be_copied = get(axes_to_be_copied,'children');
        % Identify orientation of the axes
        [az,el] = view;
        % Copy the children of the axes
        if(iscell(chilred_to_be_copied))
            copyobj(chilred_to_be_copied{2},sub_p);
        else
            copyobj(chilred_to_be_copied,sub_p);
        end
        % Close the figure
        close(f_c);
        axis vis3d;
        axis equal;
        if(exist('axis_on','var'))
            if(isequal(axis_on{p},'on'));axis on; else;axis off;end
        else
            axis on;
        end
        if(exist('colorbars','var'))
            if(isequal(colorbars{p},'on'));colorbar ; else; colorbar off; end
        else
            colorbar off;
        end
        if(exist('view_orient','var'))
            if(~isempty(view_orient{p}))
                orient = view_orient{p};
                view(orient(1),orient(2));
            end
        end 
        if(exist('subtitles','var'))
            subtitle = subtitles{p};
            title(subtitle);
        end
        pause(1);        
    end
end
if(exist('cmap','var'))
    colormap(cmap);
end
if(isfolder(tmp_path))
   rmdir(tmp_path,'s'); 
end
sgtitle(fig_title);

%% Colspan
% subplot(4,4,[1 6]) % top left
% subplot(4,4,[3 8]) % top right
% subplot(4,4,[9 10]) % bottom left (1)
% subplot(4,4,[13 14]) % bottom left (2)
% subplot(4,4,[11 16]) % bottom right
end

