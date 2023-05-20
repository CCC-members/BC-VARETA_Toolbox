function [outvar1] = eegplot(data, varargin)

outvar1 = [];

try icadefs;
    DEFAULT_FIG_COLOR = BACKCOLOR;
    BUTTON_COLOR = GUIBUTTONCOLOR;
catch
    DEFAULT_FIG_COLOR = [1 1 1];
    BUTTON_COLOR =[0.8 0.8 0.8];
end
DEFAULT_AXIS_COLOR = 'k';         % X-axis, Y-axis Color, text Color
DEFAULT_GRID_SPACING = 1;         % Grid lines every n seconds
DEFAULT_GRID_STYLE = '-';         % Grid line style
SPACING_UNITS_STRING = '';        % '\muV' for microvolt optional units for g.spacingI Ex. uV
%MAXEVENTSTRING = 10;
%DEFAULT_AXES_POSITION = [0.0964286 0.15 0.842 0.75-(MAXEVENTSTRING-5)/100];
% dimensions of main EEG axes
ORIGINAL_POSITION = [50 50 800 500];
matVers = version;
matVers = str2double(matVers(1:3));

% %%%%%%%%%%%%%%%%%%%%%%%%
% Setup inputs
% %%%%%%%%%%%%%%%%%%%%%%%%
if ~ischar(data) % If NOT a 'noui' call or a callback from uicontrols
    try
        options = varargin;
        if ~isempty( varargin )
            for i = 1:2:numel(options)
                g.(options{i}) = options{i+1};
            end
        else g= []; end
    catch
        disp('eegplot() error: calling convention {''key'', value, ... } error'); return;
    end

    % push button: create/remove window
    % ---------------------------------
    defdowncom   = 'eegplot(''defdowncom'',   gcbf);'; % push button: create/remove window
    defmotioncom = 'eegplot(''defmotioncom'', gcbf);'; % motion button: move windows or display current position
    defupcom     = 'eegplot(''defupcom'',     gcbf);';
    defctrldowncom = 'eegplot(''topoplot'',   gcbf);'; % CTRL press and motion -> do nothing by default
    defctrlupcom = ''; % CTRL press and up -> do nothing by default

    try g.srate; 		    catch, g.srate		= 256; 	end
    try g.spacing; 			catch, g.spacing	= 0; 	end
    try g.eloc_file; 		catch, g.eloc_file	= 0; 	end % 0 mean numbered
    try g.winlength; 		catch, g.winlength	= 5; 	end % Number of seconds of EEG displayed
    try g.position; 	    catch, g.position	= ORIGINAL_POSITION; 	end
    try g.title; 		    catch, g.title		= 'Scroll activity -- eegplot()'; 	end
    try g.plottitle; 		catch, g.plottitle	= ''; 	end
    try g.trialstag; 		catch, g.trialstag	= -1; 	end
    try g.winrej; 			catch, g.winrej		= []; 	end
    try g.command; 			catch, g.command	= ''; 	end
    try g.tag; 				catch, g.tag		= 'EEGPLOT'; end
    try g.xgrid;		    catch, g.xgrid		= 'off'; end
    try g.ygrid;		    catch, g.ygrid		= 'off'; end
    try g.color;		    catch, g.color		= 'off'; end
    try g.submean;			catch, g.submean	= 'off'; end
    try g.children;			catch, g.children	= 0; end
    try g.limits;		    catch, g.limits	    = [0 1000*(size(data,2)-1)/g.srate]; end
    try g.freqs;            catch, g.freqs	    = []; end  % Ramon
    try g.freqlimits;	    catch, g.freqlimits	= []; end
    try g.dispchans; 		catch, g.dispchans  = size(data,1); end
    try g.wincolor; 		catch, g.wincolor   = [ 0.7 1 0.9]; end
    try g.butlabel; 		catch, g.butlabel   = 'REJECT'; end
    try g.colmodif; 		catch, g.colmodif   = { g.wincolor }; end
    try g.scale; 		    catch, g.scale      = 'on'; end
    try g.events; 		    catch, g.events      = []; end
    try g.ploteventdur;     catch, g.ploteventdur = 'off'; end
    try g.data2;            catch, g.data2      = []; end
    try g.plotdata2;        catch, g.plotdata2 = 'off'; end
    try g.mocap;		    catch, g.mocap		= 'off'; end % nima
    try g.datastd;          catch, g.datastd = []; end %ozgur
    try g.normed;           catch, g.normed = 0; end %ozgur
    try g.envelope;         catch, g.envelope = 0; end %ozgur
    try g.maxeventstring;   catch, g.maxeventstring = 10; end % JavierLC
    try g.isfreq;           catch, g.isfreq = 0;    end % Ramon
    try g.noui;             catch, g.noui = 'off'; end
    try g.time;             catch, g.time = []; end
    if strcmpi(g.ploteventdur, 'on'), g.ploteventdur = 1; else g.ploteventdur = 0; end
    if ~ismatrix(data)
        g.trialstag = size(	data, 2);
    end

    gfields = fieldnames(g);
    for index=1:length(gfields)
        switch gfields{index}
            case {'spacing', 'srate' 'eloc_file' 'winlength' 'position' 'title' 'plottitle' ...
                    'trialstag'  'winrej' 'command' 'tag' 'xgrid' 'ygrid' 'color' 'colmodif'...
                    'freqs' 'freqlimits' 'submean' 'children' 'limits' 'dispchans' 'wincolor' ...
                    'maxeventstring' 'ploteventdur' 'butlabel' 'scale' 'events' 'data2' 'plotdata2' ...
                    'mocap' 'selectcommand' 'ctrlselectcommand' 'datastd' 'normed' 'envelope' 'isfreq' 'noui' 'time' }
            otherwise, error(['eegplot: unrecognized option: ''' gfields{index} '''' ]);
        end
    end

    % g.data=data; % never used and slows down display dramatically - Ozgur 2010

    if length(g.srate) > 1
        disp('Error: srate must be a single number'); return;
    end
    if length(g.spacing) > 1
        disp('Error: ''spacing'' must be a single number'); return;
    end
    if length(g.winlength) > 1
        disp('Error: winlength must be a single number'); return;
    end
    if ischar(g.title) > 1
        disp('Error: title must be is a string'); return;
    end
    if ischar(g.command) > 1
        disp('Error: command must be is a string'); return;
    end
    if ischar(g.tag) > 1
        disp('Error: tag must be is a string'); return;
    end
    if length(g.position) ~= 4
        disp('Error: position must be is a 4 elements array'); return;
    end
    switch lower(g.xgrid)
  	   case { 'on', 'off' }
  	   otherwise disp('Error: xgrid must be either ''on'' or ''off'''); return;
    end
    switch lower(g.ygrid)
  	   case { 'on', 'off' }
  	   otherwise disp('Error: ygrid must be either ''on'' or ''off'''); return;
    end
    switch lower(g.submean)
  	   case { 'on' 'off' }
  	   otherwise disp('Error: submean must be either ''on'' or ''off'''); return;
    end
    switch lower(g.scale)
  	   case { 'on' 'off' }
  	   otherwise disp('Error: scale must be either ''on'' or ''off'''); return;
    end

    if ~iscell(g.color)
   	   switch lower(g.color)
           case 'on', g.color = { 'k', 'm', 'c', 'b', 'g' };
           case 'off', g.color = { [ 0 0 0.4] };
           otherwise
               disp('Error: color must be either ''on'' or ''off'' or a cell array');
               return;
       end
    end
    if length(g.dispchans) > size(data,1)
  	   g.dispchans = size(data,1);
    end
    if ~iscell(g.colmodif)
        g.colmodif = { g.colmodif };
    end
    if g.maxeventstring>20 % JavierLC
        disp('Error: maxeventstring must be equal or lesser than 20'); return;
    end
    if isempty(g.time)
        g.time = fastif(g.trialstag(1) == -1, 0, 1);
    end

    % max event string;  JavierLC
    % ---------------------------------
    MAXEVENTSTRING = g.maxeventstring;
    DEFAULT_AXES_POSITION = [0.0964286 0.15 0.842 0.75-(MAXEVENTSTRING-5)/100];

    % convert color to modify into array of float
    % -------------------------------------------
    for index = 1:length(g.colmodif)
        if iscell(g.colmodif{index})
            tmpcolmodif{index} = g.colmodif{index}{1} ...
                + g.colmodif{index}{2}*10 ...
                + g.colmodif{index}{3}*100;
        else
            tmpcolmodif{index} = g.colmodif{index}(1) ...
                + g.colmodif{index}(2)*10 ...
                + g.colmodif{index}(3)*100;
        end
    end
    g.colmodif = tmpcolmodif;

    [g.chans,g.frames, tmpnb] = size(data);
    g.frames = g.frames*tmpnb;

    if g.spacing == 0
        maxindex = min(1000, g.frames);
        stds = std(data(:,1:maxindex),[],2);
        g.datastd = stds;
        stds = sort(stds);
        if length(stds) > 2
            stds = mean(stds(2:end-1));
        else
            stds = mean(stds);
        end
        g.spacing = stds*3;
        if g.spacing > 10
            g.spacing = round(g.spacing);
        end
        if g.spacing  == 0 || isnan(g.spacing)
            g.spacing = 1; % default
        end
    end

    % set defaults
    % ------------
    g.incallback = 0;
    g.winstatus = 1;
    g.setelectrode  = 0;
    [g.chans,g.frames,tmpnb] = size(data);
    g.frames = g.frames*tmpnb;
    g.nbdat = 1; % deprecated
    g.elecoffset = 0;

    % %%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare figure and axes
    % %%%%%%%%%%%%%%%%%%%%%%%%

    figh = figure('UserData', g,... % store the settings here
        'Color',DEFAULT_FIG_COLOR, 'name', g.title,...
        'MenuBar','none','tag', g.tag ,'Position',g.position, ...
        'numbertitle', 'off', 'visible', 'off', 'Units', 'Normalized');

    % Background axis
    % ---------------
    axes('tag','backeeg','parent',figh,...
        'Position',DEFAULT_AXES_POSITION,...
        'Box','off','xgrid','off', 'xaxislocation', 'top', 'Units', 'Normalized');

    % Drawing axis
    % ---------------
    YLabels = num2str((1:g.chans)');  % Use numbers as default
    YLabels = flipud(char(YLabels,' '));
    ax1 = axes('Position',DEFAULT_AXES_POSITION,...
        'userdata', data, ...% store the data here
        'tag','eegaxis','parent',figh,...%(when in g, slow down display)
        'Box','on','xgrid', g.xgrid,'ygrid', g.ygrid,...
        'gridlinestyle',DEFAULT_GRID_STYLE,...
        'Ylim',[0 (g.chans+1)*g.spacing],...
        'YTick',0:g.spacing:g.chans*g.spacing,...
        'YTickLabel', YLabels,...
        'TickLength',[.005 .005],...
        'Color','none',...
        'XColor',DEFAULT_AXIS_COLOR,...
        'YColor',DEFAULT_AXIS_COLOR);

    if ischar(g.eloc_file) || isstruct(g.eloc_file)  % Read in electrode name
        if isstruct(g.eloc_file) && length(g.eloc_file) > size(data,1)
            g.eloc_file(end) = []; % common reference channel location
        end
        eegplot('setelect', g.eloc_file, ax1);
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up uicontrols
    % %%%%%%%%%%%%%%%%%%%%%%%%%

    % positions of buttons
    posbut(1,:) = [ 0.0464    0.0254    0.0385    0.0339 ]; % <<
    posbut(2,:) = [ 0.0924    0.0254    0.0288    0.0339 ]; % <
    posbut(3,:) = [ 0.1924    0.0254    0.0299    0.0339 ]; % >
    posbut(4,:) = [ 0.2297    0.0254    0.0385    0.0339 ]; % >>
    posbut(5,:) = [ 0.1287    0.0203    0.0561    0.0390 ]; % Eposition
    posbut(7,:) = [ 0.2762    0.01    0.0582    0.0390 ]; % elec
    posbut(:,1) = posbut(:,1)+0.2;

    % Five move buttons: << < text > >>

    u(1) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position', posbut(1,:), ...
        'Tag','Pushbutton1',...
        'string','<<',...
        'Callback',['global in_callback;', ...
        'if isempty(in_callback);in_callback=1;', ...
        '    try eegplot(''drawp'',1);', ...
        '        clear global in_callback;', ...
        '    catch error_struct;', ...
        '        clear global in_callback;', ...
        '        throw(error_struct);', ...
        '    end;', ...
        'else;return;end;']);%James Desjardins 2013/Jan/22
    u(2) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position', posbut(2,:), ...
        'Tag','Pushbutton2',...
        'string','<',...
        'Callback',['global in_callback;', ...
        'if isempty(in_callback);in_callback=1;', ...
        '    try eegplot(''drawp'',2);', ...
        '        clear global in_callback;', ...
        '    catch error_struct;', ...
        '        clear global in_callback;', ...
        '        throw(error_struct);', ...
        '    end;', ...
        'else;return;end;']);%James Desjardins 2013/Jan/22
    u(5) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',[1 1 1], ...
        'Position', posbut(5,:), ...
        'Style','edit', ...
        'Tag','EPosition',...
        'string', num2str(g.time),...
        'Callback', 'eegplot(''drawp'',0);' );
    u(3) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(3,:), ...
        'Tag','Pushbutton3',...
        'string','>',...
        'Callback',['global in_callback;', ...
        'if isempty(in_callback);in_callback=1;', ...
        '    try eegplot(''drawp'',3);', ...
        '        clear global in_callback;', ...
        '    catch error_struct;', ...
        '        clear global in_callback;', ...
        '        throw(error_struct);', ...
        '    end;', ...
        'else;return;end;']);%James Desjardins 2013/Jan/22
    u(4) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(4,:), ...
        'Tag','Pushbutton4',...
        'string','>>',...
        'Callback',['global in_callback;', ...
        'if isempty(in_callback);in_callback=1;', ...
        '    try eegplot(''drawp'',4);', ...
        '        clear global in_callback;', ...
        '    catch error_struct;', ...
        '        clear global in_callback;', ...
        '        error(error_struct);', ...
        '    end;', ...
        'else;return;end;']);%James Desjardins 2013/Jan/22

    
    if ~isempty(g.events)
        u(17) = uicontrol('Parent',figh, ...
            'Units', 'normalized', ...
            'Position',posbut(17,:), ...
            'string', 'Event types', 'callback', 'eegplot(''drawlegend'', gcbf)');
    end

    for i = 1: length(u) % Matlab 2014b compatibility
        if isprop(eval(['u(' num2str(i) ')']),'Style')
            set(u(i),'Units','Normalized');
        end
    end

    % plot durations
    % --------------
    if g.ploteventdur && isfield(g.events, 'duration')
        disp(['Use menu "Display > Hide event duration" to hide colored regions ' ...
            'representing event duration']);
    end
    if isfield(g.events, 'duration')
        uimenu('Parent',m(1),'Label',fastif(g.ploteventdur, 'Hide event duration', 'Plot event duration'),'Callback', ...
            ['g = get(gcbf, ''userdata'');' ...
            'if ~g.ploteventdur' ...
            '  set(gcbo, ''label'', ''Hide event duration'');' ...
            'else' ...
            '  set(gcbo, ''label'', ''Show event duration'');' ...
            'end;' ...
            'g.ploteventdur = ~g.ploteventdur;' ...
            'set(gcbf, ''userdata'', g);' ...
            'eegplot(''drawb''); clear g;'] )
    end
   
    

    set(figh, 'userdata', g);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot EEG Data
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(ax1)
    hold on

    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot Spacing I
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    YLim = get(ax1,'Ylim');
    A = DEFAULT_AXES_POSITION;
    axes('Position',[A(1)+A(3) A(2) 1-A(1)-A(3) A(4)],'Visible','off','Ylim',YLim,'tag','eyeaxes')
    axis manual
    
    eegplot('drawp', 0);
    if g.dispchans ~= g.chans
    	   eegplot('zoom', gcf);
    end
    
    h = findobj(gcf, 'style', 'pushbutton');
    set(h, 'backgroundcolor', BUTTON_COLOR);
    h = findobj(gcf, 'tag', 'eegslider');
    set(h, 'backgroundcolor', BUTTON_COLOR);
    set(figh, 'visible', 'on');

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End Main Function
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    try
        p1 = varargin{1}; p2 = varargin{2}; p3 = varargin{3}; 
    catch
    end
    switch data
        case 'drawp' % Redraw EEG and change position

            % this test help to couple eegplot windows
            if exist('p3', 'var')
                figh = p3;
                figure(p3);
            else
                figh = gcf;                          % figure handle
            end

            if strcmp(get(figh,'tag'),'dialog')
                figh = get(figh,'UserData');
            end
            ax0 = findobj('tag','backeeg','parent',figh); % axes handle
            ax1 = findobj('tag','eegaxis','parent',figh); % axes handle

            g = get(figh,'UserData');
            data = get(ax1,'UserData');
            ESpacing = findobj('tag','ESpacing','parent',figh);   % ui handle
            EPosition = findobj('tag','EPosition','parent',figh); % ui handle
            if ~isempty(EPosition) && ~isempty(ESpacing)
                if g.trialstag(1) == -1
                    g.time    = str2num(get(EPosition,'string'));
                else
                    g.time    = str2num(get(EPosition,'string'));
                    g.time    = g.time - 1;
                end
                g.spacing = str2num(get(ESpacing,'string'));
            end

            if p1 == 1
                g.time = g.time-g.winlength*0.9;     % << subtract one window length
            elseif p1 == 2
                g.time = g.time-fastif(g.winlength>=5, round(g.winlength/5), g.winlength/5);             % < subtract one second
            elseif p1 == 3
                g.time = g.time+fastif(g.winlength>=5, round(g.winlength/5), g.winlength/5);             % > add one second
            elseif p1 == 4
                g.time = g.time+g.winlength*0.9;     % >> add one window length
            end

          	 if g.trialstag ~= -1 % time in second or in trials
          		multiplier = g.trialstag;
        	 else
          		multiplier = g.srate;
             end

             % Update edit box
             % ---------------
             g.time = max(0,min(g.time,ceil((g.frames-1)/multiplier)-g.winlength));
             if g.trialstag(1) == -1
                 set(EPosition,'string',num2str(g.time));
             else
                 set(EPosition,'string',num2str(g.time+1));
             end
             set(figh, 'userdata', g);

             lowlim = round(g.time*multiplier+1);
             highlim = round(min((g.time+g.winlength)*multiplier+2,g.frames));

             % Plot data and update axes
             % -------------------------
             if ~isempty(g.data2)
                 switch lower(g.submean) % subtract the mean ?
                     case 'on'
                         meandata = mean(g.data2(:,lowlim:highlim)');
                         if any(isnan(meandata))
                             meandata = nan_mean(g.data2(:,lowlim:highlim)');
                         end
                     otherwise, meandata = zeros(1,g.chans);
                 end
             else
                 switch lower(g.submean) % subtract the mean ?
                     case 'on'
                         meandata = mean(data(:,lowlim:highlim)');
                         if any(isnan(meandata))
                             meandata = nan_mean(data(:,lowlim:highlim)');
                         end
                     otherwise, meandata = zeros(1,g.chans);
                 end
             end
             if strcmpi(g.plotdata2, 'on')
                 hold on;
             else
                 cla(ax1);
             end

             oldspacing = g.spacing;
             if g.envelope
                 g.spacing = 0;
             end

             % plot channels whose "badchan" field is set to 1.
             % Bad channels are plotted first so that they appear behind the good
             % channels in the eegplot figure window.
             for i = 1:g.chans
                 if strcmpi(g.plotdata2, 'on')
                     tmpcolor = [ 1 0 0 ];
                 else tmpcolor = g.color{mod(i-1,length(g.color))+1};
                 end

                 if isfield(g, 'eloc_file') && isfield(g.eloc_file, 'badchan') && g.eloc_file(g.chans-i+1).badchan
                     tmpcolor = [ .85 .85 .85 ];
                     plot(ax1, data(g.chans-i+1,lowlim:highlim) -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), ...
                         'color', tmpcolor, 'clipping','on')
                     plot(ax1, 1,mean(data(g.chans-i+1,lowlim:highlim) -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing),2),'<r','MarkerFaceColor','r','MarkerSize',6);
                 end

             end

             % plot good channels on top of bad channels (if g.eloc_file(i).badchan = 0... or there is no bad channel information)
             if strcmpi(g.plotdata2, 'on')
                 tmpcolor = [ 1 0 0 ];
             else tmpcolor = g.color{mod(g.chans-i,length(g.color))+1};
             end

             %        keyboard;
             if (isfield(g, 'eloc_file') && isfield(g.eloc_file, 'badchan') && ~g.eloc_file(g.chans-i+1).badchan) || ...
                     (~isfield(g, 'eloc_file')) || (~isfield(g.eloc_file, 'badchan'))
                 plot(ax1, bsxfun(@plus, data(end:-1:1,lowlim:highlim), g.spacing*[1:g.chans]'-meandata(end:-1:1)')' + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), ...
                     'color', tmpcolor, 'clipping','on');                
             end

                         
             g.spacing = oldspacing;
             set(ax1, 'Xlim',[1 g.winlength*multiplier],...
                 'XTick',1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1);
            
             set(ax1, 'XTickLabel', num2str((g.time:DEFAULT_GRID_SPACING:g.time+g.winlength)'));
             

             % ordinates: even if all elec are plotted, some may be hidden
             set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );             
             eegplot('drawb'); % draw background first
             
        case 'drawb' % Draw background ******************************************************
            % Redraw EEG and change position

            ax0 = findobj('tag','backeeg','parent',gcf); % axes handle
            ax1 = findobj('tag','eegaxis','parent',gcf); % axes handle
            ylims=ylim(ax0);

            g = get(gcf,'UserData');  % Data (Note: this could also be global)

            % Plot data and update axes
            cla(ax0);
            hold(ax0, 'on');            
            % ordinates: even if all elec are plotted, some may be hidden
            set(ax0, 'ylim',ylims );
            set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );

        case 'loadelect' % load channels
            [inputname,inputpath] = uigetfile('*','Channel locations file');
            if inputname == 0
                return
            end
            if ~exist([ inputpath inputname ])
                error('no such file');
            end

            AXH0 = findobj('tag','eegaxis','parent',gcf);
            eegplot('setelect',[ inputpath inputname ],AXH0);
            return;

        case 'setelect'
            % Set channels
            eloc_file = p1;
            axeshand = p2;
            outvar1 = 1;
            if isempty(eloc_file)
                outvar1 = 0;
                return
            end
            tmplocs = readlocs(eloc_file);
            YLabels = { tmplocs.labels };
            YLabels = strvcat(YLabels);

            YLabels = flipud(char(YLabels,' '));
            set(axeshand,'YTickLabel',YLabels)        
        otherwise
            error(['Error - invalid eegplot() parameter: ',data])
    end
end
