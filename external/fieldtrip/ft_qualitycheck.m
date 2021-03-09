function [varargout] = ft_qualitycheck(cfg)

% FT_QUALITYCHECK performs a quality inspection of a given MEG/EEG dataset,
% stores (.mat), and visualizes the result (.png and .pdf).
%
% This function segments the data into 10-second pieces and performs the
% following analyses:
%  1) reads the properties of the dataset
%  2) computes the headpositions and distance covered from recording onset (CTF only)
%  3) computes the mean, max, min, and range of the signal amplitude
%  4) detects trigger events
%  5) detects jump artifacts
%  6) computes the powerspectrum
%  7) estimates the low-frequency (<2 Hz) and line noise (~50 Hz)
%
% Use as
%   [info, timelock, freq, summary, headpos] = ft_qualitycheck(cfg)
% where info contains the dataset properties, timelock the timelocked data,
% freq the powerspectra, summary the mean descriptives, and headpos the
% headpositions throughout the recording
%
% The configuration should contain:
%   cfg.dataset = a string (e.g. 'dataset.ds')
%
% The following parameters can be used:
%   cfg.analyze   = string, 'yes' or 'no' to analyze the dataset (default = 'yes')
%   cfg.savemat   = string, 'yes' or 'no' to save the analysis (default = 'yes')
%   cfg.matfile   = string, filename (e.g. 'previousoutput.mat'), preferably in combination
%                    with analyze = 'no'
%   cfg.visualize = string, 'yes' or 'no' to visualize the analysis (default = 'yes')
%   cfg.saveplot  = string, 'yes' or 'no' to save the visualization (default = 'yes')
%   cfg.linefreq  = scalar, frequency of power line (default = 50)
%   cfg.plotunit  = scalar, the length of time to be plotted in one panel (default = 3600)
%
% See also FT_PREPROCESSING, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT

% Copyright (C) 2010-2011, Arjen Stolk, Bram Daams, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id:%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.analyze   = ft_getopt(cfg, 'analyze',   'yes');
cfg.savemat   = ft_getopt(cfg, 'savemat',   'yes');
cfg.matfile   = ft_getopt(cfg, 'matfile',   []);
cfg.visualize = ft_getopt(cfg, 'visualize', 'yes');
cfg.saveplot  = ft_getopt(cfg, 'saveplot',  'yes');
cfg.linefreq  = ft_getopt(cfg, 'linefreq',  50);
cfg.plotunit  = ft_getopt(cfg, 'plotunit',  3600);

%% ANALYSIS
if strcmp(cfg.analyze,'yes')
  tic
  
  % checks
  cfg   = ft_checkconfig(cfg, 'dataset2files', 'yes'); % translate into datafile+headerfile
  
  % these will be replaced by more appropriate values
  info.filename    = cfg.dataset;
  info.datasetname = 'unknown';
  info.starttime   = 'unknown';
  info.startdate   = 'unknown';
  info.stoptime    = 'unknown';
  info.stopdate    = 'unknown';
  
  % the exportname is also used in the cron job
  exportname = qualitycheck_exportname(cfg.dataset);
  [iseeg, ismeg, isctf, fltp] = filetyper(cfg.dataset);
  if isctf
    try
      % update the info fields
      info = read_ctf_hist(cfg.dataset);
    end
  end
  
  % add info
  info.event                  = ft_read_event(cfg.dataset);
  info.hdr                    = ft_read_header(cfg.dataset);
  info.filetype               = fltp;
  
  % trial definition
  cfgdef                      = [];
  cfgdef.dataset              = cfg.dataset;
  cfgdef.trialdef.triallength = 10;
  cfgdef.continuous           = 'yes';
  cfgdef                      = ft_definetrial(cfgdef);
  ntrials                     = size(cfgdef.trl,1)-1; % remove last trial
  timeunit                    = cfgdef.trialdef.triallength;
  
  % channelselection for jump detection (all) and for FFT (brain)
  if ismeg
    allchans                   = ft_channelselection({'MEG','MEGREF'}, info.hdr.label);
    chans                      = ft_channelselection('MEG', info.hdr.label); % brain
    allchanindx                = match_str(info.hdr.label, allchans);
    chanindx                   = match_str(chans, allchans);
    jumpthreshold              = 1e-10;
  elseif iseeg
    allchans                   = ft_channelselection('EEG', info.hdr.label);
    if isempty(allchans)
      % some EEG systems and data files use non-standard channel names that are not detected automatically
      ft_warning('no EEG channels detected, selecting all channels');
      allchans = info.hdr.label;
    end
    chans                      = allchans;  % brain
    allchanindx                = match_str(info.hdr.label, allchans);
    chanindx                   = match_str(chans, allchans);
    jumpthreshold              = 1e4;
  end
  
  % find headcoil channels
  if isctf % this fails for older CTF data sets
    Nx = strmatch('HLC0011', info.hdr.label); % x nasion coil
    Ny = strmatch('HLC0012', info.hdr.label); % y nasion
    Nz = strmatch('HLC0013', info.hdr.label); % z nasion
    Lx = strmatch('HLC0021', info.hdr.label); % x left coil
    Ly = strmatch('HLC0022', info.hdr.label); % y left
    Lz = strmatch('HLC0023', info.hdr.label); % z left
    Rx = strmatch('HLC0031', info.hdr.label); % x right coil
    Ry = strmatch('HLC0032', info.hdr.label); % y right
    Rz = strmatch('HLC0033', info.hdr.label); % z right
    headpos.dimord = 'chan_time';
    headpos.time   = (timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2);
    headpos.label  = {'Nx';'Ny';'Nz';'Lx';'Ly';'Lz';'Rx';'Ry';'Rz'};
    headpos.avg    = NaN(length(headpos.label), ntrials);
    headpos.grad   = info.hdr.grad;
    
    if numel(cat(1,Nx,Ny,Nz,Lx,Ly,Lz,Rx,Ry,Rz))==9
      hasheadpos = true;
    else
      hasheadpos = false;
    end
    
  end % if
  
  % analysis settings
  cfgredef             = [];
  cfgredef.length      = 1;
  cfgredef.overlap     = 0;
  
  cfgfreq              = [];
  cfgfreq.output       = 'pow';
  cfgfreq.channel      = allchans;
  cfgfreq.method       = 'mtmfft';
  cfgfreq.taper        = 'hanning';
  cfgfreq.keeptrials   = 'no';
  cfgfreq.foilim       = [0 min(info.hdr.Fs/2, 400)];
  
  % output variables
  timelock.dimord = 'chan_time';
  timelock.label  = allchans;
  timelock.time   = (timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2);
  timelock.avg    = NaN(length(allchans), ntrials); % updated in loop
  timelock.median = NaN(length(allchans), ntrials); % updated in loop
  timelock.jumps  = NaN(length(allchans), ntrials); % updated in loop
  timelock.range  = NaN(length(allchans), ntrials); % updated in loop
  timelock.min    = NaN(length(allchans), ntrials); % updated in loop
  timelock.max    = NaN(length(allchans), ntrials); % updated in loop
  
  freq.dimord     = 'chan_freq_time';
  freq.label      = allchans;
  freq.freq       = (cfgfreq.foilim(1):cfgfreq.foilim(2));
  freq.time       = (timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2);
  freq.powspctrm  = NaN(length(allchans), length(freq.freq), ntrials); % updated in loop
  
  summary.dimord  = 'chan_time';
  summary.time    = (timeunit-timeunit/2:timeunit:timeunit*ntrials-timeunit/2);
  summary.label   = {'Mean';'Median';'Min';'Max';'Range';'HmotionN';'HmotionL';'HmotionR';'LowFreqPower';'LineFreqPower';'Jumps'};
  summary.avg     = NaN(length(summary.label), ntrials); % updated in loop
  
  % try add gradiometer info
  if isfield(info.hdr, 'grad')
    timelock.grad = info.hdr.grad;
    freq.grad     = info.hdr.grad;
    summary.grad  = info.hdr.grad;
  end
  
  
  % process trial by trial
  for t = 1:ntrials
    fprintf('analyzing trial %s of %s \n', num2str(t), num2str(ntrials));
    
    % preprocess
    cfgpreproc     = cfgdef;
    cfgpreproc.trl = cfgdef.trl(t,:);
    data           = ft_preprocessing(cfgpreproc); clear cfgpreproc;
    
    % determine headposition
    if isctf && hasheadpos
      headpos.avg(1,t) = mean(data.trial{1,1}(Nx,:) * 100);  % meter to cm
      headpos.avg(2,t) = mean(data.trial{1,1}(Ny,:) * 100);
      headpos.avg(3,t) = mean(data.trial{1,1}(Nz,:) * 100);
      headpos.avg(4,t) = mean(data.trial{1,1}(Lx,:) * 100);
      headpos.avg(5,t) = mean(data.trial{1,1}(Ly,:) * 100);
      headpos.avg(6,t) = mean(data.trial{1,1}(Lz,:) * 100);
      headpos.avg(7,t) = mean(data.trial{1,1}(Rx,:) * 100);
      headpos.avg(8,t) = mean(data.trial{1,1}(Ry,:) * 100);
      headpos.avg(9,t) = mean(data.trial{1,1}(Rz,:) * 100);
    end
    
    % update values
    timelock.avg(:,t)    = mean(data.trial{1}(allchanindx,:),2);
    timelock.median(:,t) = median(data.trial{1}(allchanindx,:),2);
    timelock.range(:,t)  = max(data.trial{1}(allchanindx,:),[],2) - min(data.trial{1}(allchanindx,:),[],2);
    timelock.min(:,t)    = min(data.trial{1}(allchanindx,:),[],2);
    timelock.max(:,t)    = max(data.trial{1}(allchanindx,:),[],2);
    
    % detect jumps
    for c = 1:size(data.trial{1}(allchanindx,:),1)
      timelock.jumps(c,t) = length(find(diff(data.trial{1,1}(allchanindx(c),:)) > jumpthreshold));
    end
    
    % FFT and noise estimation
    redef                 = ft_redefinetrial(cfgredef, data); clear data;
    FFT                   = ft_freqanalysis(cfgfreq, redef); clear redef;
    freq.powspctrm(:,:,t) = FFT.powspctrm;
    summary.avg(9,t)      = mean(mean(findpower(0,  2,  FFT, chanindx))); % Low Freq Power
    summary.avg(10,t)     = mean(mean(findpower(cfg.linefreq-1, cfg.linefreq+1, FFT, chanindx))); clear FFT; % Line Freq Power
    
    toc
  end % end of trial loop
  
  % determine headmotion: distance from initial trial (in cm)
  if isctf && hasheadpos
    summary.avg(6,:) = sqrt(sum((headpos.avg(1:3,:)-repmat(headpos.avg(1:3,1),1,size(headpos.avg,2))).^2,1)); % N
    summary.avg(7,:) = sqrt(sum((headpos.avg(4:6,:)-repmat(headpos.avg(4:6,1),1,size(headpos.avg,2))).^2,1)); % L
    summary.avg(8,:) = sqrt(sum((headpos.avg(7:9,:)-repmat(headpos.avg(7:9,1),1,size(headpos.avg,2))).^2,1)); % R
  end
  
  % summarize/mean and store variables of brain info only
  summary.avg(1,:)   = mean(timelock.avg(chanindx,:),1);
  summary.avg(2,:)   = mean(timelock.median(chanindx,:),1);
  summary.avg(3,:)   = mean(timelock.min(chanindx,:),1);
  summary.avg(4,:)   = mean(timelock.max(chanindx,:),1);
  summary.avg(5,:)   = mean(timelock.range(chanindx,:),1);
  summary.avg(11,:)  = mean(timelock.jumps(chanindx,:),1);
  
  % save to .mat
  if strcmp(cfg.savemat, 'yes')
    if isctf && hasheadpos
      headpos.cfg        = cfg;
      save(exportname, 'info','timelock','freq','summary','headpos');
    else
      save(exportname, 'info','timelock','freq','summary');
    end
  end
  
end % end of analysis

%% VISUALIZATION
if strcmp(cfg.visualize, 'yes')
  
  % load data
  if strcmp(cfg.analyze, 'no')
    if ~isempty(cfg.matfile)
      exportname = cfg.matfile;
    else
      exportname = qualitycheck_exportname(cfg.dataset);
    end
    fprintf('loading %s \n', exportname);
    load(exportname);
  end
  
  % determine number of 1-hour plots to be made
  nplots = ceil(length(freq.time)/(cfg.plotunit/10));
  
  % create GUI-like figure(s)
  for p = 1:nplots
    fprintf('visualizing %s of %s \n', num2str(p), num2str(nplots));
    toi = [p*cfg.plotunit-(cfg.plotunit-5) p*cfg.plotunit-5]; % select 1-hour chunks
    
    tmpcfg.latency = toi;
    temp_timelock  = ft_selectdata(tmpcfg, timelock);
    temp_freq      = ft_selectdata(tmpcfg, freq);
    temp_summary   = ft_selectdata(tmpcfg, summary);
    if exist('headpos','var')
      temp_headpos  = ft_selectdata(tmpcfg, headpos);
      draw_figure(info, temp_timelock, temp_freq, temp_summary, temp_headpos, toi);
      clear temp_timelock; clear temp_freq; clear temp_summary; clear temp_headpos; clear toi;
    else
      draw_figure(info, temp_timelock, temp_freq, temp_summary, toi);
      clear temp_timelock; clear temp_freq; clear temp_summary; clear toi;
    end
    
    % export to .PNG and .PDF
    if strcmp(cfg.saveplot, 'yes')
      [pathstr,name,extr] = fileparts(exportname);
      if p == 1
        exportfilename = name;
      else
        exportfilename = strcat(name,'_pt',num2str(p));
      end
      fprintf('exporting %s of %s \n', num2str(p), num2str(nplots));
      set(gcf, 'PaperType', 'a4');
      print(gcf, '-dpng', strcat(exportfilename,'.png'));
      orient landscape;
      print(gcf, '-dpdf', strcat(exportfilename,'.pdf'));
      close
    end
  end % end of nplots
end % end of visualization

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble history timelock   % add the input cfg to multiple outputs
ft_postamble history freq       % add the input cfg to multiple outputs
ft_postamble history summary    % add the input cfg to multiple outputs

%% VARARGOUT
if nargout>0
  mOutputArgs{1} = info;
  mOutputArgs{2} = timelock;
  mOutputArgs{3} = freq;
  mOutputArgs{4} = summary;
  try
    mOutputArgs{5} = headpos;
  end
  [varargout{1:nargout}] = mOutputArgs{:};
  clearvars -except varargout
else
  clear
end

%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = clipat(x, v, v2)
v = [v v2]; % clip between value v and v2
if length(v) == 1
  x(x>v) = v;
elseif length(v) == 2
  x(x<v(1)) = v(1);
  x(x>v(2)) = v(2);
end

%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function [iseeg, ismeg, isctf, fltp] = filetyper(dataset)
fltp = ft_filetype(dataset);
iseeg = ft_filetype(dataset,'brainvision_eeg') | ...
  ft_filetype(dataset,'ns_eeg') | ...
  ft_filetype(dataset,'bci2000_dat') | ...
  ft_filetype(dataset,'neuroprax_eeg') | ...
  ft_filetype(dataset,'egi_sbin') | ...
  ft_filetype(dataset,'biosemi_bdf');
ismeg = ft_filetype(dataset,'ctf_ds') | ...
  ft_filetype(dataset,'4d') | ...
  ft_filetype(dataset,'neuromag_fif') | ...
  ft_filetype(dataset,'itab_raw');
isctf = ft_filetype(dataset, 'ctf_ds');
if ~ismeg && ~iseeg % if none found, try less strict checks
  [p, f, ext] = fileparts(dataset);
  if strcmp(ext, '.eeg')
    fltp = 'brainvision_eeg';
    iseeg = 1;
  elseif strcmp(ext, '.bdf')
    fltp = 'biosemi_bdf';
    iseeg = 1;
  elseif strcmp(ext, '.ds')
    fltp = 'ctf_ds';
    ismeg = 1;
  else % otherwise use eeg settings for stability reasons
    iseeg = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function [power, freq] = findpower(low, high, freqinput, chans)
% replace value with the index of the nearest bin
xmin  = nearest(getsubfield(freqinput, 'freq'), low);
xmax  = nearest(getsubfield(freqinput, 'freq'), high);
% select the freq range
power = freqinput.powspctrm(chans, xmin:xmax);
freq  = freqinput.freq(:, xmin:xmax);

%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function draw_figure(varargin)
% deal with input
if nargin == 6
  info     = varargin{1};
  timelock = varargin{2};
  freq     = varargin{3};
  summary  = varargin{4};
  headpos  = varargin{5};
  toi      = varargin{6};
elseif nargin == 5
  info     = varargin{1};
  timelock = varargin{2};
  freq     = varargin{3};
  summary  = varargin{4};
  toi      = varargin{5};
end

% determine whether it is EEG or MEG
if isfield(info, 'filename') % supported as of January 2018
  filename = info.filename;
elseif isfield(timelock.cfg, 'dataset')
  filename = timelock.cfg.dataset;
elseif isfield(info.hdr.orig, 'FileName')
  filename = info.hdr.orig.FileName;
elseif exist('headpos','var') && isfield(headpos.cfg.previous, 'dataset')
  filename = headpos.cfg.previous.dataset;
else
  error('could not determine the filename');
end
[iseeg, ismeg, isctf, fltp] = filetyper(filename);

if ismeg
  scaling = 1e15; % assuming data is in T and needs to become fT
  powscaling = scaling^2;
  ylab = 'fT';
elseif iseeg
  scaling = 1e0; % assuming data is in muV already
  powscaling = scaling^2;
  ylab = '\muV';
end

% PARENT FIGURE
h.MainFigure = figure(...
  'MenuBar','none',...
  'Name','ft_qualitycheck',...
  'Units','normalized',...
  'color','white',...
  'Position',[0.01 0.01 .99 .99]); % nearly fullscreen

if strcmp(info.startdate,'unknown')
  tmp = 'unknown';
else
  [d,w] = weekday(info.startdate);
  tmp = [w ' ' info.startdate];
end

h.MainText = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',tmp,...
  'Backgroundcolor','white',...
  'Position',[.06 .96 .15 .02]);

h.MainText2 = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String','Jump artefacts',...
  'Backgroundcolor','white',...
  'Position',[.08 .46 .12 .02]);

h.MainText3 = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String','Mean powerspectrum',...
  'Backgroundcolor','white',...
  'Position',[.4 .3 .15 .02]);

h.MainText4 = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String','Timecourses',...
  'Backgroundcolor','white',...
  'Position',[.5 .96 .11 .02]);

h.MainText5 = uicontrol(...
  'Parent',h.MainFigure,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String','Events',...
  'Backgroundcolor','white',...
  'Position',[.81 .3 .06 .02]);

% HEADMOTION PANEL
h.HmotionPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.01 .5 .25 .47]);

h.DataText = uicontrol(...
  'Parent',h.HmotionPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',info.datasetname,...
  'Backgroundcolor','white',...
  'Position',[.01 .85 .99 .1]);

h.TimeText = uicontrol(...
  'Parent',h.HmotionPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',[info.starttime ' - ' info.stoptime],...
  'Backgroundcolor','white',...
  'Position',[.01 .78 .99 .1]);

if ismeg
  allchans = ft_senslabel(ft_senstype(timelock));
  misschans = setdiff(ft_channelselection('MEG', info.hdr.label), allchans);
  nchans = num2str(size(ft_channelselection('MEG', info.hdr.label),1));
else
  misschans = '';
  nchans = num2str(size(ft_channelselection('EEG', info.hdr.label),1));
end

h.DataText2 = uicontrol(...
  'Parent',h.HmotionPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',[fltp ', fs: ' num2str(info.hdr.Fs) ', nchans: ' nchans],...
  'Backgroundcolor','white',...
  'Position',[.01 .71 .99 .1]);

h.DataText3 = uicontrol(...
  'Parent',h.HmotionPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',['missing chans: ' misschans'],...
  'Backgroundcolor','white',...
  'Position',[.01 .64 .99 .1]);

% boxplot headmotion (*10; cm-> mm) per coil
if exist('headpos','var')
  h.HmotionAxes = axes(...
    'Parent',h.HmotionPanel,...
    'Units','normalized',...
    'color','white',...
    'Position',[.05 .08 .9 .52]);
  
  hmotions = ([summary.avg(8,:)' summary.avg(7,:)' summary.avg(6,:)'])*10;
  boxplot(h.HmotionAxes, hmotions, 'orientation', 'horizontal', 'notch', 'on');
  set(h.HmotionAxes,'YTick',1:3);
  set(h.HmotionAxes,'YTickLabel',{'R','L','N'});
  xlim(h.HmotionAxes, [0 10]);
  xlabel(h.HmotionAxes, 'Headmotion from start [mm]');
end

% TIMECOURSE PANEL
h.SignalPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.28 .34 .71 .63]);

h.SignalAxes = axes(...
  'Parent',h.SignalPanel,...
  'Units','normalized',...
  'color','white',...
  'Position',[.08 .36 .89 .3]);

h.LinenoiseAxes = axes(...
  'Parent',h.SignalPanel,...
  'Units','normalized',...
  'color','white',...
  'Position',[.08 .23 .89 .1]);

h.LowfreqnoiseAxes = axes(...
  'Parent',h.SignalPanel,...
  'Units','normalized',...
  'color','white',...
  'Position',[.08 .1 .89 .1]);

% plot hmotion timecourses per coil (*10; cm-> mm)
if exist('headpos','var')
  h.HmotionTimecourseAxes = axes(...
    'Parent',h.SignalPanel,...
    'Units','normalized',...
    'color','white',...
    'Position',[.08 .73 .89 .22]);
  
  plot(h.HmotionTimecourseAxes, summary.time, clipat(summary.avg(6,:)*10, 0, 10), ...
    summary.time, clipat(summary.avg(7,:)*10, 0, 10), ...
    summary.time, clipat(summary.avg(8,:)*10, 0, 10), 'LineWidth',2);
  ylim(h.HmotionTimecourseAxes,[0 10]);
  ylabel(h.HmotionTimecourseAxes, 'Coil distance [mm]');
  xlim(h.HmotionTimecourseAxes,toi);
  grid(h.HmotionTimecourseAxes,'on');
  legend(h.HmotionTimecourseAxes, 'N','L','R');
end

% plot mean and range of the raw signal
plot(h.SignalAxes, summary.time, summary.avg(5,:)*scaling, summary.time, summary.avg(1,:)*scaling, 'LineWidth', 2);
set(h.SignalAxes,'Nextplot','add');
plot(h.SignalAxes, summary.time, summary.avg(3,:)*scaling, summary.time, summary.avg(4,:)*scaling, 'LineWidth', 1, 'Color', [255/255 127/255 39/255]);
grid(h.SignalAxes,'on');
ylabel(h.SignalAxes, ['Amplitude [' ylab ']']);
xlim(h.SignalAxes,toi);
legend(h.SignalAxes,'Range','Mean','Min','Max');
set(h.SignalAxes,'XTickLabel','');

% plot linenoise
semilogy(h.LinenoiseAxes, summary.time, clipat(summary.avg(10,:)*powscaling, 1e2, 1e4), 'LineWidth',2);
grid(h.LinenoiseAxes,'on');
legend(h.LinenoiseAxes, ['LineFreq [' ylab '^2/Hz]']);
set(h.LinenoiseAxes,'XTickLabel','');
xlim(h.LinenoiseAxes,toi);
ylim(h.LinenoiseAxes,[1e2 1e4]); % before april 28th this was 1e0 - 1e3

% plot lowfreqnoise
semilogy(h.LowfreqnoiseAxes, summary.time, clipat(summary.avg(9,:)*powscaling, 1e10, 1e12), 'LineWidth',2);
grid(h.LowfreqnoiseAxes,'on');
xlim(h.LowfreqnoiseAxes,toi);
ylim(h.LowfreqnoiseAxes,[1e10 1e12]);
legend(h.LowfreqnoiseAxes, ['LowFreq [' ylab '^2/Hz]']);
xlabel(h.LowfreqnoiseAxes, 'Time [seconds]'); % before april 28th this was 1e0 - 1e10

% EVENT PANEL
h.EventPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.7 .01 .29 .3]);

% event details
eventtypes    = {};
eventtriggers = {};
eventvalues   = {};

if ~isempty(info.event)
  [a,b,c] = unique({info.event.type});
  for j=1:length(a)
    eventtypes{j,1}    = a{j};
    eventtriggers{j,1} = sum(c==j);
    eventvalues{j,1}   = length(unique([info.event(c==j).value]));
  end
end
if isempty(eventtypes)
  eventtypes{1,1} = 'no triggers found';
end

h.EventText = uicontrol(...
  'Parent',h.EventPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',['Types'; ' '; eventtypes],...
  'Backgroundcolor','white',...
  'Position',[.05 .05 .4 .85]);

h.EventText2 = uicontrol(...
  'Parent',h.EventPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',['Triggers'; ' '; eventtriggers],...
  'Backgroundcolor','white',...
  'Position',[.55 .05 .2 .85]);

h.EventText3 = uicontrol(...
  'Parent',h.EventPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',['Values'; ' '; eventvalues],...
  'Backgroundcolor','white',...
  'Position',[.8 .05 .15 .85]);

% POWERSPECTRUM PANEL
h.SpectrumPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.28 .01 .4 .3]);

h.SpectrumAxes = axes(...
  'Parent',h.SpectrumPanel,...
  'Units','normalized',...
  'color','white',...
  'Position',[.15 .2 .8 .7]);

% plot powerspectrum
loglog(h.SpectrumAxes, freq.freq, mean(mean(freq.powspctrm,1),3)*powscaling,'r','LineWidth',2);
xlabel(h.SpectrumAxes, 'Frequency [Hz]');
ylabel(h.SpectrumAxes, ['Power [' ylab '^2/Hz]']);

% ARTEFACT PANEL
h.JumpPanel = uipanel(...
  'Parent',h.MainFigure,...
  'Units','normalized',...
  'Backgroundcolor','white',...
  'Position',[.01 .01 .25 .46]);

% jump details
jumpchans  = {};
jumpcounts = {};
[jumps,i] = find(timelock.jumps>0); % find all jumps
[a,b,c] = unique(jumps);
for j=1:length(a)
  jumpchans{j,1}  = timelock.label{a(j)};
  jumpcounts{j,1} = sum(c==j);
end
if isempty(jumpchans)
  jumpchans{1,1} = 'no jumps detected';
end

h.JumpText = uicontrol(...
  'Parent',h.JumpPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',jumpchans,...
  'Backgroundcolor','white',...
  'Position',[.15 .5 .25 .4]);

h.JumpText2 = uicontrol(...
  'Parent',h.JumpPanel,...
  'Style','text',...
  'Units','normalized',...
  'FontSize',10,...
  'String',jumpcounts,...
  'Backgroundcolor','white',...
  'Position',[.65 .5 .2 .4]);

% plot jumps on the dewar sensors
if ismeg
  h.TopoMEG = axes(...
    'Parent',h.JumpPanel,...
    'color','white',...
    'Units','normalized',...
    'Position',[0.4 0.05 0.55 0.4]);
  
  MEGchans                 = ft_channelselection('MEG', timelock.label);
  MEGchanindx              = match_str(timelock.label, MEGchans);
  cfgtopo                  = [];
  cfgtopo.marker           = 'off';
  cfgtopo.colorbar         = 'no';
  cfgtopo.comment          = 'no';
  cfgtopo.style            = 'blank';
  cfgtopo.layout           = ft_prepare_layout([], timelock);
  cfgtopo.highlight        = 'on';
  cfgtopo.highlightsymbol  = '.';
  cfgtopo.highlightsize    = 14;
  cfgtopo.highlightchannel = find(sum(timelock.jumps(MEGchanindx,:),2)>0);
  data.label               = MEGchans;
  data.powspctrm           = sum(timelock.jumps(MEGchanindx,:),2);
  data.dimord              = 'chan_freq';
  data.freq                = 1;
  axes(h.TopoMEG);
  ft_topoplotTFR(cfgtopo, data); % FIXME should this not be ft_topoplotER?
  clear data
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ready longer than the total length requested
        begsample  = cfg.trl(i,1);
        endsample  = cfg.trl(i,2);
        offset     = cfg.trl(i,3);
        begpadding = 0;
        endpadding = 0;
        if padding > 0
          ft_warning('no padding applied because the padding duration is shorter than the trial');
        end
      else
        switch cfg.paddir
          case 'both'
            % begpadding+nsamples+endpadding = total length of raw data that will be read
            begpadding = ceil((padding-nsamples)/2);
            endpadding = floor((padding-nsamples)/2);
          case 'left'
            begpadding = padding-nsamples;
            endpadding = 0;
          case 'right'
            begpadding = 0;
            endpadding = padding-nsamples;
          otherwise
            ft_error('unsupported requested direction of padding');
        end
        
        if strcmp(cfg.padtype, 'data')
          begsample  = cfg.trl(i,1) - begpadding;
          endsample  = cfg.trl(i,2) + endpadding;
        else
          % padding will be done below
          begsample  = cfg.trl(i,1);
          endsample  = cfg.trl(i,2);
        end
        if begsample<1
          ft_warning('cannot apply enough padding at begin of file');
          begpadding = begpadding - (1 - begsample);
          begsample  = 1;
        end
        if endsample>(hdr.nSamples*hdr.nTrials)
          ft_warning('cannot apply enough padding at end of file');
          endpadding = endpadding - (endsample - hdr.nSamples*hdr.nTrials);
          endsample  = hdr.nSamples*hdr.nTrials;
        end
        offset = cfg.trl(i,3) - begpadding;
      end
      
      % read the raw data with padding on both sides of the trial - this
      % includes datapadding
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', rawindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
      
      % convert the data to another numeric precision, i.e. double, single or int32
      if ~isempty(cfg.precision)
        dat = cast(dat, cfg.precision);
      end
      
      % pad in case of no datapadding
      if ~strcmp(cfg.padtype, 'data')
        dat = ft_preproc_padding(dat, cfg.padtype, begpadding, endpadding);
        tim = offset2time(offset+begpadding, hdr.Fs, size(dat,2));
      else
        tim = offset2time(offset, hdr.Fs, size(dat,2));
      end
      
      % do the preprocessing on the padded trial data and remove the padding after filtering
      [cutdat{i}, label, time{i}, cfg] = preproc(dat, hdr.label(rawindx), tim, cfg, begpadding, endpadding);
      
      if isfield(cfg, 'export') && ~isempty(cfg.export)
        % write the processed data to an original manufacturer format file
        newhdr        = [];
        newhdr.Fs     = hdr.Fs;
        newhdr.label  = label;
        newhdr.nChans = length(newhdr.label);
        % only append for the second and consecutive trials
        ft_write_data(cfg.export.dataset, cutdat{i}, 'dataformat', cfg.export.dataformat, 'header', newhdr, 'append', i~=1);
        if nargout==0
          % don't keep the processed data in memory
          cutdat(i) = [];
        end
      end
      
    end % for all trials
    ft_progress('close');
    
    % don't keep hdr.orig in the output data if it is too large
    % hdr.orig can be large when caching data from specific file formats, such as bci2000_dat and mega_neurone
    if isfield(hdr, 'orig')
      s = hdr.orig;
      s = whos('s');
      if s.bytes>10240
        hdr = rmfield(hdr, 'orig');
      end
    end
    
    dataout                    = [];
    dataout.hdr                = hdr;                  % header details of the datafile
    dataout.label              = label;                % labels of channels that have been read, can be different from labels in file due to montage
    dataout.time               = time;                 % vector with the timeaxis for each individual trial
    dataout.trial              = cutdat;
    dataout.fsample            = hdr.Fs;
    dataout.sampleinfo         = cfg.trl(:,1:2);
    if size(cfg.trl,2) > 3
      dataout.trialinfo      = cfg.trl(:,4:end);
    end
    if isfield(hdr, 'grad')
      dataout.grad             = hdr.grad;             % MEG gradiometer information in header (f.e. headerformat = 'ctf_ds')
    end
    if isfield(hdr, 'elec')
      dataout.elec             = hdr.elec;             % EEG electrode information in header (f.e. headerformat = 'neuromag_fif')
    end
    if isfield(hdr, 'opto')
      dataout.opto             = hdr.opto;             % NIRS optode information in header (f.e. headerformat = 'artinis')
    end
    
  end % for all channel groups
  
end % if hasdata

if strcmp(cfg.updatesens, 'yes')
  % updating the sensor descriptions can be done on basis of the montage or the rereference settings
  if ~isempty(cfg.montage) && ~isequal(cfg.montage, 'no')
    montage = cfg.montage;
  elseif strcmp(cfg.reref, 'yes')
    if strcmp(cfg.refmethod, 'bipolar') || strcmp(cfg.refmethod, 'avg')
      tmpcfg = keepfields(cfg, {'refmethod', 'implicitref', 'refchannel', 'channel'});
      tmpcfg.showcallinfo = 'no';
      montage = ft_prepare_montage(tmpcfg, data);
    else
      % do not update the sensor description
      montage = [];
    end
  else
    % do not update the sensor description
    montage = [];
  end
  
  if ~isempty(montage)
    % apply the linear projection also to the sensor description
    if issubfield(montage, 'type')
      bname = montage.type;
    else
      bname = 'preproc';
    end
    if isfield(dataout, 'grad')
      ft_info('applying the montage to the grad structure\n');
      dataout.grad = ft_apply_montage(dataout.grad, montage, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
    end
    if isfield(dataout, 'elec')
      ft_info('applying the montage to the elec structure\n');
      dataout.elec = ft_apply_montage(dataout.elec, montage, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
    end
    if isfield(dataout, 'opto')
      ft_info('applying the montage to the opto structure\n');
      dataout.opto = ft_apply_montage(dataout.opto, montage, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
    end
  end
end % if updatesens

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = dataout;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

                              axis\n', cfg.symmetry);
  % expand the number of parameters from one (3) to two dipoles (6)
  sourcemodel.pos = sourcemodel.pos(:,expand) .* repmat(mirror, size(sourcemodel.pos,1), 1);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance sourcemodel
ft_postamble history    sourcemodel
ft_postamble savevar    sourcemodel

%--------------------------------------------------------------
% helper function for basedonmri method to determine the inside
% returns a boolean vector
function inside = getinside(pos, mask)

% it might be that the box with the points does not completely fit into the
% mask
dim = size(mask);
sel = find(pos(:,1)<1 |  pos(:,1)>dim(1) | ...
  pos(:,2)<1 |  pos(:,2)>dim(2) | ...
  pos(:,3)<1 |  pos(:,3)>dim(3));
if isempty(sel)
  % use the efficient implementation
  inside = mask(sub2ind(dim, pos(:,1), pos(:,2), pos(:,3)));
else
  % only loop over the points that can be dealt with
  inside = false(size(pos,1), 1);
  for i=setdiff(1:size(pos,1), sel(:)')
    inside(i) = mask(pos(i,1), pos(i,2), pos(i,3));
  end
end

%--------------------------------------------------------------------------
% helper function to return the fieldnames of the boolean fields in a
% segmentation, should work both for volumetric and for source
function fn = booleanfields(mri)

fn = fieldnames(mri);
isboolean = false(1,numel(fn));
for i=1:numel(fn)
  if islogical(mri.(fn{i})) && isequal(numel(mri.(fn{i})),prod(mri.dim))
    isboolean(i) = true;
  end
end
fn  = fn(isboolean);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
      case 'q'
        again = 0;

      otherwise
        ft_warning('invalid button (%d)', k);
    end
  end % while again
  % remember this set of polygons as the outline
  outline = polygon;

  % convert the sensor positions into a layout structure
  layout.pos = pos;
  nchans = size(pos,1);
  for i=1:nchans
    layout.label{i,1} = num2str(i);
  end
  % compute the width and height for multiplotting
  d = dist(pos');
  for i=1:nchans
    d(i,i) = inf; % exclude the diagonal
  end
  mindist = min(d(:));
  layout.width  = ones(nchans,1) * mindist * 0.8;
  layout.height = ones(nchans,1) * mindist * 0.6;
  % add the polygons that describe the mask and outline
  layout.mask    = mask;
  layout.outline = outline;

  finalhelp = [ ...
    '-----------------------------------------------------------------------------------\n' ...
    'you should update the channel labels, and check the width and height in the output layout\n' ...
    ];
  fprintf(finalhelp);
  fprintf('\n');

else
  ft_error('no layout detected, please specify cfg.layout')
end

% make the subset as specified in cfg.channel
cfg.channel = ft_channelselection(cfg.channel, setdiff(layout.label, {'COMNT', 'SCALE'}));  % COMNT and SCALE are not really channels
chansel = match_str(layout.label, cat(1, cfg.channel(:), 'COMNT', 'SCALE'));                % include COMNT and SCALE, keep all channels in the order of the layout
% return the layout for the subset of channels
layout.pos    = layout.pos(chansel,:);
layout.label  = layout.label(chansel);
if strcmpi(cfg.style, '2d')
  % width and height only apply to the 2D layout
  layout.width  = layout.width(chansel);
  layout.height = layout.height(chansel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overrule the width and height when specified by the user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.width)
  layout.width(:) = cfg.width;
end
if ~isempty(cfg.height)
  layout.height(:) = cfg.height;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether the outline and mask are available, create them if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~isfield(layout, 'outline') || ~isfield(layout, 'mask')) && ~strcmpi(cfg.style, '3d')
  % the reason to check for style=3d rather than 2d is that cfg.style is also an option in ft_topoplotER and ft_topoplotTFR
  % the style option of that function easily "leaks" into here, causing the default 2d not to be selected at the top

  if strcmp(cfg.outline, 'circle') || strcmp(cfg.mask, 'circle')
    % Scale the electrode positions to fit within a unit circle, i.e. electrode radius = 0.45
    ind_scale = find(strcmp('SCALE', layout.label));
    ind_comnt = find(strcmp('COMNT', layout.label));
    sel = setdiff(1:length(layout.label), [ind_scale ind_comnt]); % these are excluded for scaling
    x = layout.pos(sel,1);
    y = layout.pos(sel,2);
    if istrue(cfg.center)
      % the following centers all electrodes around zero
      xrange = range(x);
      yrange = range(y);
    else
      % the following prevent topography distortion in case electrodes are not evenly distributed over the whole head
      xrange = 2*( max(max(x),abs(min(x)) ));
      yrange = 2*( max(max(y),abs(min(y)) ));
    end
    if xrange==0
      xrange = 1;
    end
    if yrange==0
      yrange = 1;
    end
    % First scale the width and height of the box for multiplotting
    layout.width  = layout.width./xrange;
    layout.height = layout.height./yrange;
    % Then shift and scale the electrode positions
    layout.pos(:,1) = 0.8*((layout.pos(:,1)-min(x))/xrange-0.5);
    layout.pos(:,2) = 0.8*((layout.pos(:,2)-min(y))/yrange-0.5);
  end

  if ~isfield(layout, 'outline') && ischar(cfg.outline)
    switch cfg.outline
      case 'circle'
        layout.outline = outline_circle();
      case 'convex'
        layout.outline = outline_convex(layout);
      case 'square'
        layout.outline = outline_square(layout);
      case {'headshape', 'mri'}
        % the configuration should contain the headshape or mri
        % the (segmented) mri will be converted into a headshape on the fly
        hsoutline = outline_headshape(cfg, sens); % used for mask if possible
        layout.outline = hsoutline;
      otherwise
        layout.outline = {};
    end
  end

  if ~isfield(layout, 'mask') && ischar(cfg.mask)
    switch cfg.mask
      case 'circle'
        layout.mask = outline_circle();
        layout.mask = layout.mask(1); % the first is the circle, the others are nose and ears
      case 'convex'
        layout.mask = outline_convex(layout);
      case 'square'
        layout.mask = outline_square(layout);
      case {'headshape', 'mri'}
        % the configuration should contain the headshape or mri
        % the (segmented) mri will be converted into a headshape on the fly
        if isequal(cfg.mask,cfg.outline) && exist('hsoutline','var')
          layout.mask = hsoutline;
        else
          layout.mask = outline_headshape(cfg, sens);
        end
      otherwise
        layout.mask = {};
    end
  end

end % create outline if style=2d


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the montage, e.g. convert from monopolar to bipolar channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(cfg.montage, 'no')
  Nold = length(cfg.montage.labelold);
  Nnew = length(cfg.montage.labelnew);
  for i=1:Nnew
    cfg.montage.tra(i,:) = abs(cfg.montage.tra(i,:));
    cfg.montage.tra(i,:) = cfg.montage.tra(i,:) ./ sum(cfg.montage.tra(i,:));
  end
  % pretend it is a sensor structure, this achieves averaging after channel matching
  tmp.tra   = layout.pos;
  tmp.label = layout.label;
  new = ft_apply_montage(tmp, cfg.montage);
  layout.pos   = new.tra;
  layout.label = new.label;
  % do the same for the width and height
  tmp.tra = layout.width(:);
  new = ft_apply_montage(tmp, cfg.montage);
  layout.width = new.tra;
  tmp.tra = layout.height(:);
  new = ft_apply_montage(tmp, cfg.montage);
  layout.height = new.tra;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add axes positions for comments and scale information if required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if skipcomnt || ~isequal(cfg.commentpos, 'layout')
  % remove the comnt entry
  sel = find(strcmp('COMNT', layout.label));
  layout.label(sel)  = [];
  layout.pos(sel,:)  = [];
  layout.width(sel)  = [];
  layout.height(sel) = [];
end

if skipscale || ~isequal(cfg.scalepos, 'layout')
  % remove the scale entry
  sel = find(strcmp('SCALE', layout.label));
  layout.label(sel)  = [];
  layout.pos(sel,:)  = [];
  layout.width(sel)  = [];
  layout.height(sel) = [];
end

if (~skipcomnt || ~skipscale) && ~strcmpi(cfg.style, '3d')
  % this is used for the placement of the comment and scale
  pos = layout.pos;
  if isfield(layout, 'outline')
    pos = cat(1, pos, layout.outline{:});
  end
  if isfield(layout, 'mask')
    pos = cat(1, pos, layout.mask{:});
  end
  width  = mean(layout.width);
  height = mean(layout.height);
  middle = @(x) min(x) + (max(x)-min(x))/2;
end

if ~skipcomnt && ~any(strcmp('COMNT', layout.label)) && ~strcmpi(cfg.style, '3d') && ~isequal(cfg.commentpos, 'title')
  % add a placeholder for the comment in the desired location
  if strcmp(cfg.commentpos, 'layout')
    cfg.commentpos = 'leftbottom'; % set the default position
  end
  if strcmp(cfg.commentpos, 'lefttop')
    layout.pos(end+1,:) = [min(pos(:,1))-width/2 max(pos(:,2))+height/2];
  elseif strcmp(cfg.commentpos, 'leftbottom')
    layout.pos(end+1,:) = [min(pos(:,1))-width/2 min(pos(:,2))-height/2];
  elseif strcmp(cfg.commentpos, 'middletop')
    layout.pos(end+1,:) = [middle(pos(:,1)) max(pos(:,2))+height/2];
  elseif strcmp(cfg.commentpos, 'middlebottom')
    layout.pos(end+1,:) = [middle(pos(:,1)) min(pos(:,2))-height/2];
  elseif strcmp(cfg.commentpos, 'righttop')
    layout.pos(end+1,:) = [max(pos(:,1))+width/2 max(pos(:,2))+height/2];
  elseif strcmp(cfg.commentpos, 'rightbottom')
    layout.pos(end+1,:) = [max(pos(:,1))+width/2 min(pos(:,2))-height/2];
  elseif isnumeric(cfg.commentpos)
    layout.pos(end+1,:) = cfg.commentpos;
  else
    ft_error('invalid specification of cfg.commentpos');
  end
  layout.label{end+1}  = 'COMNT';
  layout.width(end+1)  = width;
  layout.height(end+1) = height;
end

if ~skipscale && ~any(strcmp('SCALE', layout.label)) && ~strcmpi(cfg.style, '3d')
  % add a placeholder for the scale in the desired location
  if strcmp(cfg.scalepos, 'layout')
    cfg.scalepos = 'rightbottom'; % set the default position
  end
  if strcmp(cfg.scalepos, 'lefttop')
    layout.pos(end+1,:) = [min(pos(:,1))-width/2 max(pos(:,2))+height/2];
  elseif strcmp(cfg.scalepos, 'leftbottom')
    layout.pos(end+1,:) = [min(pos(:,1))-width/2 min(pos(:,2))-height/2];
  elseif strcmp(cfg.scalepos, 'middletop')
    layout.pos(end+1,:) = [middle(pos(:,1)) max(pos(:,2))+height/2];
  elseif strcmp(cfg.scalepos, 'middlebottom')
    layout.pos(end+1,:) = [middle(pos(:,1)) min(pos(:,2))-height/2];
  elseif strcmp(cfg.scalepos, 'righttop')
    layout.pos(end+1,:) = [max(pos(:,1))+width/2 max(pos(:,2))+height/2];
  elseif strcmp(cfg.scalepos, 'rightbottom')
    layout.pos(end+1,:) = [max(pos(:,1))+width/2 min(pos(:,2))-height/2];
  elseif isnumeric(cfg.scalepos)
    layout.pos(end+1,:) = cfg.scalepos;
  else
    ft_error('invalid specification of cfg.scalepos');
  end
  layout.label{end+1}  = 'SCALE';
  layout.width(end+1)  = width;
  layout.height(end+1) = height;
end

% these should be represented in a column vector (see bug 1909 -roevdmei)
layout.label = layout.label(:);
% the width and height are not present in a 3D layout as used in SPM
if ~strcmpi(cfg.style, '3d')
  layout.width  = layout.width(:);
  layout.height = layout.height(:);
end

% to plot the layout for debugging, you can use this code snippet
if strcmp(cfg.feedback, 'yes') && ~strcmpi(cfg.style, '3d')
  tmpcfg = [];
  tmpcfg.layout = layout;
  ft_layoutplot(tmpcfg); % FIXME this should use ft_plot_layout
end

% to write the layout to a .mat or text file, you can use this code snippet
if ~isempty(cfg.output) && ~strcmpi(cfg.style, '3d')
  ft_info('writing layout to ''%s''\n', cfg.output);
  if strcmpi(cfg.output((end-3):end), '.mat')
    save(cfg.output,'layout');
  else
    fid = fopen(cfg.output, 'wt');
    for i=1:numel(layout.label)
      fprintf(fid, '%d %f %f %f %f %s\n', i, layout.pos(i,1), layout.pos(i,2), ...
        layout.width(i), layout.height(i), layout.label{i});
    end
    fclose(fid);
  end
elseif ~isempty(cfg.output) && strcmpi(cfg.style, '3d')
  % the layout file format does not support 3D positions, furthermore for
  % a 3D layout the width and height are currently set to NaN
  ft_error('writing a 3D layout to an output file is not supported');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble history layout
ft_postamble savevar layout


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% read the layout information from the ascii file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
return % function readlay


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% convert 3D electrode positions into 2D layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = sens2lay(sens, rotatez, projmethod, style, overlap, viewpoint, boxchannel)

% remove the balancing from the sensor definition, e.g. 3rd order gradients, PCA-cleaned data or ICA projections
% this not only removed the linear projections, but also ensures that the channel labels are correctly named

if isfield(sens, 'chanposold')
  chanposold = sens.chanposold;
else
  chanposold = [];
end
if isfield(sens, 'balance') && ~strcmp(sens.balance.current, 'none')
  sens = undobalancing(sens);
  if size(chanposold, 1) == numel(sens.label)
    sens.chanpos = chanposold;
  end
  % In case not all the locations have NaNs it might still be useful to plot them
  % But perhaps it'd be better to have any
elseif any(all(isnan(sens.chanpos)))
  [sel1, sel2] = match_str(sens.label, sens.labelold);
  sens.chanpos = chanposold(sel2, :);
  sens.label   = sens.labelold(sel2);
end

ft_info('creating layout for %s system\n', ft_senstype(sens));

% apply rotation, but only if viewpoint is not used specifically
if isempty(viewpoint)
  if isempty(rotatez)
    switch ft_senstype(sens)
      case {'ctf151', 'ctf275', 'bti148', 'bti248', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar', 'yokogawa160', 'yokogawa160_planar', 'yokogawa64', 'yokogawa64_planar', 'yokogawa440', 'yokogawa440_planar', 'magnetometer', 'meg'}
        rotatez = 90;
      case {'neuromag122', 'neuromag306'}
        rotatez = 0;
      case 'electrode'
        rotatez = 90;
      otherwise
        rotatez = 0;
    end
  end
  sens.chanpos = ft_warp_apply(rotate([0 0 rotatez]), sens.chanpos, 'homogenous');
end

% determine the 3D channel positions
pos   = sens.chanpos;
label = sens.label;

if strcmpi(style, '3d')
  layout.pos   = pos;
  layout.label = label;

else
  if isempty(viewpoint)
    % projection other than viewpoint-specific orthographic projection is requested, use elproj
    prj = elproj(pos, projmethod);

  else
    % apply viewpoint-specific orthographic projection

    % determine auto view if requested
    if strcmp(viewpoint, 'auto')
      % simple automatic determination of the 'ideal' viewpoint
      % first, depth or not: if Xvar (l/r axis) is bigger than both Yvar (post/ant axis) and Zvar (top/bottom axis), it's a depth
      % if yes, superior (screw inferior) is more appriorate if Yvar > Zvar, otherwise posterior (screw anterior)
      % if no, it's left/right, sign of mean(X) indicates which side the grid is on (note, for interhemispheric grids, both left/right (doenst) work)
      posvar = var(pos);
      if (posvar(1)>posvar(2)) && (posvar(1)>posvar(3)) % if they're roughly equal, it's likely a diagonal depth, and any view would (not) work
        if posvar(2)>posvar(3)
          viewpoint = 'superior';
        else
          viewpoint = 'posterior';
        end
      else
        if sign(mean(pos(:,1))) == -1
          viewpoint = 'left';
        else
          viewpoint = 'right';
        end
      end
    end

    % 3D to 2D
    prj = getorthoviewpos(pos, sens.coordsys, viewpoint);
  end

  % this copy will be used to determine the minimum distance between channels
  % we need a copy because prj retains the original positions, and
  % prjForDist might need to be changed if the user wants to keep
  % overlapping channels
  % also subselect channels for computing width/height if requested by boxchannel
  boxchannel = ft_channelselection(boxchannel, label);
  boxchansel = match_str(label, boxchannel);
  prjForDist = prj(boxchansel,:);

  % check whether many channels occupy identical positions, if so shift them around if requested
  if size(unique(prjForDist,'rows'),1) / size(prjForDist,1) < 0.8
    if strcmp(overlap, 'shift')
      ft_warning('the specified sensor configuration has many overlapping channels, creating a layout by shifting them around - use a template layout for better control over the positioning');
      prj = shiftxy(prj', 0.2)';
      prjForDist = prj(boxchansel,:);
    elseif strcmp(overlap, 'no')
      ft_error('the specified sensor configuration has many overlapping channels, you specified not to allow that');
    elseif strcmp(overlap, 'keep')
      prjForDist = unique(prj(boxchansel,:), 'rows');
    else
      ft_error('unknown value for cfg.overlap = ''%s''', overlap);
    end
  end

  d = dist(prjForDist');
  d(logical(eye(size(d)))) = inf;

  % This is a fix for .sfp files, containing positions of 'fiducial
  % electrodes'. Their presence determines the minimum distance between
  % projected electrode positions, leading to very small boxes.
  % This problem has been detected and reported by Matt Mollison.
  % FIXME: consider changing the box-size being determined by mindist
  % by a mean distance or so; this leads to overlapping boxes, but that
  % also exists for some .lay files.
  if any(strmatch('Fid', label(boxchansel)))
    tmpsel = strmatch('Fid', label(boxchansel));
    d(tmpsel, :) = inf;
    d(:, tmpsel) = inf;
  end

  if any(isfinite(d(:)))
    % take mindist as the median of the first quartile of closest channel pairs with non-zero distance
    mindist = min(d); % get closest neighbour for all channels
    mindist = sort(mindist(mindist>1e-6),'ascend');
    mindist = mindist(1:round(numel(label(boxchansel))/4));
    mindist = median(mindist);
    %%% /workaround - make a safe guess to detect iEEG until a better solution is found
    if any(strcmp(ft_senstype(sens),{'eeg','unknown'})) && ~isempty(viewpoint) && (numel(boxchannel) ~= numel(label))
      mindist = min(d);
      mindist = mindist(mindist>1e-6); % allows for substantially more overlap than the above
      mindist = median(mindist);
    end
    %%% \workaround - make a safe guess to detect iEEG until a better solution is found
  else
    mindist = eps; % not sure this is a good value but it's just to prevent crashes when
    % the EEG sensor definition is meaningless
  end

  X = prj(:,1);
  Y = prj(:,2);
  Width  = ones(size(X)) * mindist * 0.8;
  Height = ones(size(X)) * mindist * 0.6;
  layout.pos    = [X Y];
  layout.width  = Width;
  layout.height = Height;
  layout.pos    = prj;
  layout.label  = label;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% convert 2D optode positions into 2D layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout = opto2lay(opto, label, rotatez)

if isempty(rotatez)
  rotatez = 90;
end

layout = [];
layout.pos    = [];
layout.label  = {};
layout.width  = [];
layout.height = [];

% NIRS channels are named 'RxY - TxZ [wavelength]'
[rxnames, rem] = strtok(label, {'-', ' '});
[txnames, rem] = strtok(rem,   {'-', ' '});

for i=1:numel(label)
  % create positions halfway between transmitter and receiver
  rxid = ismember(opto.optolabel, rxnames(i));
  txid = ismember(opto.optolabel, txnames(i));
  layout.pos(i, :) = opto.optopos(rxid, :)/2 + opto.optopos(txid, :)/2;
end

layout.label  = label;
layout.width  = ones(numel(label),1);
layout.height = ones(numel(label),1);

% apply the rotation around the z-axis
layout.pos = ft_warp_apply(rotate([0 0 rotatez]), layout.pos, 'homogenous');

% prevent the circle-with-ears-and-nose to be added
layout.outline = {};

% construct a mask for topographic interpolation
pos1 = layout.pos; pos1(:,1) = pos1(:,1)-layout.width;
pos2 = layout.pos; pos2(:,1) = pos2(:,1)+layout.width;
pos3 = layout.pos; pos3(:,2) = pos3(:,2)-layout.height;
pos4 = layout.pos; pos4(:,2) = pos4(:,2)+layout.height;
pos = [pos1; pos2; pos3; pos4];
indx = convhull(pos);
layout.mask{1} = pos(indx,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% shift 2D positions around so that the minimum distance between any pair
% is mindist
%
% Credit for this code goes to Laurence Hunt at UCL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xy = shiftxy(xy, mindist)

x = xy(1,:);
y = xy(2,:);

l=1;
i=1; % filler
mindist = mindist/0.999; % limits the number of loops
while (~isempty(i) && l<50)
  xdiff = repmat(x,length(x),1) - repmat(x',1,length(x));
  ydiff = repmat(y,length(y),1) - repmat(y',1,length(y));
  xydist= sqrt(xdiff.^2 + ydiff.^2); %euclidean distance between all sensor pairs

  [i,j] = find(xydist<mindist*0.999);
  rm=(i<=j); i(rm)=[]; j(rm)=[]; %only look at i>j

  for m = 1:length(i)
    if (xydist(i(m),j(m)) == 0)
      diffvec = [mindist./sqrt(2) mindist./sqrt(2)];
    else
      xydiff = [xdiff(i(m),j(m)) ydiff(i(m),j(m))];
      diffvec = xydiff.*mindist./xydist(i(m),j(m)) - xydiff;
    end
    x(i(m)) = x(i(m)) - diffvec(1)/2;
    y(i(m)) = y(i(m)) - diffvec(2)/2;
    x(j(m)) = x(j(m)) + diffvec(1)/2;
    y(j(m)) = y(j(m)) + diffvec(2)/2;
  end
  l = l+1;
end

xy = [x; y];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to obtain XY pos from XYZ pos as orthographic projections depending on
% the viewpoint and coordsys. See also ELPROJ and COORDSYS2LABEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getorthoviewpos(pos, coordsys, viewpoint)
% see also

if size(pos,2)~=3
  ft_error('XYZ coordinates are required to obtain the orthographic projections based on a viewpoint')
end

% create view(az,el) transformation matrix
switch coordsys
  case {'ras' 'itab' 'neuromag' 'acpc' 'spm' 'mni' 'tal'}
    switch viewpoint
      case 'left'
        transmat = viewmtx(-90, 0);
      case 'right'
        transmat = viewmtx(90, 0);
      case 'topleft'
        transmat = viewmtx(-90, 45);
      case 'topright'
        transmat = viewmtx(90, 45);
      case 'superior'
        transmat = viewmtx(0, 90);
      case 'inferior'
        transmat = viewmtx(180, -90);
      case 'posterior'
        transmat = viewmtx(0, 0);
      case 'anterior'
        transmat = viewmtx(180, 0);
      otherwise
        ft_error('orthographic projection using viewpoint "%s" is not supported', viewpoint)
    end % switch viewpoint
  case {'als' 'ctf' '4d' 'bti'}
    switch viewpoint
      case 'left'
        transmat = viewmtx(180, 0);
      case 'right'
        transmat = viewmtx(0, 0);
      case 'topleft'
        transmat = viewmtx(180, 45);
      case 'topright'
        transmat = viewmtx(0, 45);
      case 'superior'
        transmat = viewmtx(-90, 90);
      case 'inferior'
        transmat = viewmtx(90, -90);
      case 'posterior'
        transmat = viewmtx(-90, 0);
      case 'anterior'
        transmat = viewmtx(90, 0);
      otherwise
        ft_error('orthographic projection using viewpoint "%s" is not supported', viewpoint)
    end % switch viewpoint
  otherwise
    ft_error('orthographic projection using coordinate system "%s" is not supported', coordsys)
end % switch coordsys

% extract xy
pos      = ft_warp_apply(transmat, pos, 'homogenous');
pos      = pos(:,[1 2]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate a square outline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_square(layout)

% get index of all relevant channels
ind1 = strcmp(layout.label,'COMNT');
ind2 = strcmp(layout.label,'SCALE');
ind  = ~(ind1 | ind2);

x = layout.pos(ind,1);
y = layout.pos(ind,2);
w = layout.width(ind);
h = layout.height(ind);

% determine the bounding box and add a little space around it
space = min(mean(w), mean(h))/4;
xmin = min(x - 0.5*w - space);
xmax = max(x + 0.5*w + space);
ymin = min(y - 0.5*h - space);
ymax = max(y + 0.5*h + space);

% construct the outline
outline = {[
  xmin ymax
  xmax ymax
  xmax ymin
  xmin ymin
  xmin ymax
  ]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate an outline from the boundary/convex hull of pos+width/height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_circle()
% create the default "circle with triangle" to resemble the head
% note that the electrode positions should be scaled accordingly
rmax  = 0.5;
l     = 0:2*pi/100:2*pi;
HeadX = cos(l).*rmax;
HeadY = sin(l).*rmax;
NoseX = [0.18*rmax 0 -0.18*rmax];
NoseY = [rmax-.004 rmax*1.15 rmax-.004];
EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
% Define the outline of the head, ears and nose
outline{1} = [HeadX(:) HeadY(:)];
outline{2} = [NoseX(:) NoseY(:)];
outline{3} = [ EarX(:)  EarY(:)];
outline{4} = [-EarX(:)  EarY(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate an outline from the boundary/convex hull of pos+width/height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_convex(layout)

% get index of all relevant channels
ind1 = strcmp(layout.label,'COMNT');
ind2 = strcmp(layout.label,'SCALE');
ind  = ~(ind1 | ind2);
% get position of all channel 'boxes' and draw boundary around it, or, in older matlabs, a convex hull
x = layout.pos(ind,1);
y = layout.pos(ind,2);
w = layout.width(ind);
h = layout.height(ind);
boxpos = [
  x - (w/2) y - (h/2);  % lb
  x - (w/2) y + (h/2);  % lt
  x + (w/2) y - (h/2);  % rb
  x + (w/2) y + (h/2);  % rt
  ];
if ft_platform_supports('boundary')
  k = boundary(boxpos,.2);
  outline{1} = boxpos(k,:);
else
  outline{1} = boxpos(convhull(boxpos(:,1),boxpos(:,2)),:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION generate an outline from the headshape or anatomical mri
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline = outline_headshape(cfg, sens)

if ~isempty(cfg.headshape)
  if ischar(cfg.headshape) && exist(cfg.headshape, 'file')
    ft_info('reading headshape from file %s\n', cfg.headshape);
    outlbase = ft_read_headshape(cfg.headshape);
  elseif isstruct(cfg.headshape)
    outlbase = cfg.headshape;
  else
    ft_error('incorrect specification of cfg.headshape')
  end
elseif ~isempty(cfg.mri)
  if ischar(cfg.mri) && exist(cfg.mri, 'file')
    ft_info('reading MRI from file %s\n', cfg.mri);
    outlbase = ft_read_mri(cfg.mri);
  elseif ft_datatype(cfg.mri, 'volume')
    outlbase = cfg.mri;
  else
    ft_error('incorrect specification of cfg.mri')
  end
  % create mesh from anatomical field, this will be used as headshape below
  cfgpm = [];
  cfgpm.method      = 'projectmesh';
  cfgpm.tissue      = 'brain';
  cfgpm.numvertices = 1e5;
  outlbase = ft_prepare_mesh(cfgpm, outlbase);
end

% check that we have the right data in outlbase
assert(isfield(outlbase, 'pos'), 'the headshape does not contain any vertices')

% check coordinate system of outlbase
assert(isfield(outlbase, 'coordsys'), 'no coordsys field found in headshape/mri, use ft_determine_coordsys')
assert(isfield(sens, 'coordsys'), 'no coordsys field found in sensor structure, use ft_determine_coordsys')
assert(isequal(outlbase.coordsys, sens.coordsys), 'the coordinate system of headshape/mri does not match that of sensors')

% match head geometry units with that of the sensors
outlbase = ft_convert_units(outlbase, sens.unit);

% there can be multiple meshes, e.g. left and right hemispheres
outline = cell(size(outlbase));
for i=1:numel(outlbase)
  % generate outline based on matlab version
  if ft_platform_supports('boundary')

    % extract points indicating brain
    braincoords = outlbase(i).pos;
    % apply projection and extract XY
    if isequal(cfg.projection,'orthographic') && ~isempty(cfg.viewpoint)
      braincoords = getorthoviewpos(braincoords, outlbase(i).coordsys, cfg.viewpoint);
    else
      % project identically as in sens2lay using cfg.rotate and elproj (this should be kept identical to sens2lay)
      if isempty(cfg.rotate)
        switch ft_senstype(sens)
          case {'ctf151', 'ctf275', 'bti148', 'bti248', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar', 'yokogawa160', 'yokogawa160_planar', 'yokogawa64', 'yokogawa64_planar', 'yokogawa440', 'yokogawa440_planar', 'magnetometer', 'meg'}
            rotatez = 90;
          case {'neuromag122', 'neuromag306'}
            rotatez = 0;
          case 'electrode'
            rotatez = 90;
          otherwise
            rotatez = 0;
        end
      end
      braincoords = ft_warp_apply(rotate([0 0 rotatez]), braincoords, 'homogenous');
      braincoords = elproj(braincoords, cfg.projection);
    end

    % get outline
    k = boundary(braincoords,.8);
    outline{i} = braincoords(k,:);

  else % matlab version fallback

    % plot mesh in rotated view, rotate, screencap, and trace frame to generate outline
    if isequal(cfg.projection,'orthographic') && ~isempty(cfg.viewpoint)
      outlbase(i).pos = getorthoviewpos(outlbase(i).pos, outlbase(i).coordsys, cfg.viewpoint);
    else
      % project identically as in sens2lay using cfg.rotate and elproj (this should be kept identical to sens2lay)
      if isempty(cfg.rotate)
        switch ft_senstype(sens)
          case {'ctf151', 'ctf275', 'bti148', 'bti248', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar', 'yokogawa160', 'yokogawa160_planar', 'yokogawa64', 'yokogawa64_planar', 'yokogawa440', 'yokogawa440_planar', 'magnetomete