function  eegplot(EEG)

%%
%% Setup inputs
%%
data = EEG.data;
time = EEG.times;
srate = EEG.srate;
chanlocs = EEG.chanlocs;
nbchan = EEG.nbchan;
pnts = EEG.pnts;
winlength	= 5; 	 % Number of seconds of EEG displayed

%%
%% getting spacing
%%
maxindex = min(1000, pnts);
stds = std(data(:,1:maxindex),[],2);
stds = sort(stds);
if length(stds) > 2
    stds = mean(stds(2:end-1));
else
    stds = mean(stds);
end
spacing = stds*3;
if spacing > 10
    spacing = round(spacing);
end
if spacing  == 0 || isnan(spacing)
    spacing = 1; % default
end

%%
%% Prepare figure and axes
%%
figh = figure('Color',[1 1 1], 'MenuBar','none','Position',[50 50 800 500], ...
    'numbertitle', 'off', 'visible', 'on', 'Units', 'Normalized');

%%
%% Drawing axis
%%
ax1 = axes('parent',figh,...%(when in g, slow down display)
    'Box','on','xgrid', 'off','ygrid', 'off',...
    'gridlinestyle','-',...
    'Ylim',[0 (nbchan+1)*spacing],...
    'YTick',0:spacing:nbchan*spacing,...
    'TickLength',[.005 .005],...
    'Color','none',...
    'XColor','k',...
    'YColor','k');

YLabels = {chanlocs.labels}';
YLabels = strvcat(YLabels);
YLabels = flipud(char(YLabels,' '));
set(ax1,'YTickLabel',YLabels)

hold on

%%
%% Plot EEG Data
%%
lowlim = round(time*srate+1);
highlim = round(min((time+winlength)*srate+2,pnts));
plot(ax1, bsxfun(@plus, data(end:-1:1,lowlim:highlim), spacing*(1:nbchan)')', ...
    'color', 'b', 'clipping','on');
set(ax1, 'Xlim',[1 winlength*srate],...
    'XTick',1:srate:winlength*srate+1);
set(ax1, 'XTickLabel', num2str((time:1:time+winlength)'));

% %% Wizard options
% if p1 == 1
%     time = time-winlength*0.9;     % << subtract one window length
% elseif p1 == 2
%     time = time-fastif(winlength>=5, round(winlength/5), winlength/5);             % < subtract one second
% elseif p1 == 3
%     time = time+fastif(winlength>=5, round(winlength/5), winlength/5);             % > add one second
% elseif p1 == 4
%     time = time+winlength*0.9;     % >> add one window length
% end
% % Update edit box
% % ---------------
% time = max(0,min(time,ceil((pnts-1)/srate)-winlength));
% set(EPosition,'string',num2str(time));
% set(figh, 'userdata', g);



