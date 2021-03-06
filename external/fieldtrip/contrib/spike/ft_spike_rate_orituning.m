function [Stat] = ft_spike_rate_orituning(cfg, varargin)

% FT_SPIKE_RATE_ORITUNING computes a model of the firing rate as a function
% of orientation or direction.
%
% Use as
%   [stat] = ft_spike_rate_tuning(cfg, rate1, rate2, ... rateN)
%
% The inputs RATE should be the output from FT_SPIKE_RATE. 
%
% Configurations:
%   cfg.stimuli  = should be an 1 x nConditions array of orientations or
%                  directions in radians
%                  varargin{i} corresponds to cfg.stimuli(i)
%   cfg.method   = model to apply, implemented are 'orientation' and 'direction'
%
% Outputs:
%   stat.ang       = mean angle of orientation / direction (1 x nUnits)
%   stat.osi       = orientation selectivity index (Womelsdorf et al., 2012,
%                    PNAS), that is resultant length.
%                    if cfg.method = 'orientation', then orientations are
%                    first projected on the unit circle.
%   stat.di        = direction index, 1 - min/max response

% FIXME: models for contrast etc.

% Copyright (C) 2010, Martin Vinck
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
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance varargin
ft_preamble trackconfig

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'stimuli', 'method'});
cfg = ft_checkopt(cfg, 'stimuli', 'doublevector');
cfg = ft_checkopt(cfg, 'method', 'char', {'orientation', 'direction'});

if length(varargin)<2, error('can only compute orituning if multiple inputs are specified'); end 
  
% check whether trials were kept in the rate function
for k = 1:length(varargin)
  try
    varargin{k} = rmfield(varargin{k}, 'trial');
  end
end

Tune        = ft_appendtimelock([],varargin{:});
Tune.avg    = Tune.trial;

% get the unique stimuli and the number of directions, these should match
stimuli     = cfg.stimuli(:);
nConditions = length(varargin);
nStimuli    = length(stimuli);
if nConditions~=nStimuli, error('Length of cfg.stimuli should match number of data inputs'); end

nUnits      = size(Tune.avg,2);

% change the directions so it is a circular variable with range 2*pi
if strcmp(cfg.method,'orientation')
  if (max(stimuli)-min(stimuli))>pi
    error('If cfg.tuningtype is "orientation", cfg.stimuli should have range of pi');
  end
  stimuli = stimuli*2; % convert to make it circular (see Womelsdorf et al. 2012, PNAS).
  
  if (max(stimuli)-min(stimuli))<0.5*pi
    warning('Orientations have a range < 1/2 pi. Are you sure this is correct?. Stats will be biased');
  end
  
  % compute the directionality index
  Stat.di = 1 - min(Tune.avg)./max(Tune.avg);
  
  % transform the data into complex numbers to compute resultant length
  z = exp(1i*stimuli(:)*ones(1,nUnits));
  sumAvg = sum(Tune.avg);
  z = Tune.avg.*z./(sumAvg(ones(nStimuli,1),:));
  Stat.osi  = abs(sum(z));
  
  % make preferred angle modulo 2pi, convert it back to range pi and convert to rad or deg
  prefAngle           = angle(nansum(z));
  prefAngle           = mod(prefAngle,2*pi);
  Stat.ang      = prefAngle/2;
elseif strcmp(cfg.method,'direction')
  if (max(stimuli)-min(stimuli))>2*pi
    error('Directions has a range > 2*pi. Are you sure this is radians and not degrees?')
  end
  if (max(stimuli)-min(stimuli))<0.5*pi
    warning('Directions have a range < 1/2 pi. Are you sure this is correct?');
  end
  
  % compute the directionality index
  Stat.di = 1 - min(Tune.avg)./max(Tune.avg);
  
  % transform the data into complex numbers to compute resultant length
  z = exp(1i*stimuli(:)*ones(1,nUnits));
  sumAvg = sum(Tune.avg);
  z = Tune.avg.*z./(sumAvg(ones(nStimuli,1),:));
  Stat.osi  = abs(sum(z));
  
  % make preferred angle modulo 2pi, convert it back to range pi and convert to rad or deg
  prefAngle      = angle(nansum(z));
  Stat.ang = mod(prefAngle,2*pi);
end

Stat.label   = Tune.label;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   Tune
ft_postamble provenance Stat
ft_postamble history    Stat

                                                                                                                                                                                                        Latency,endTrialLatency);

% check which trials will be used based on the latency
overlaps      = endTrialLatency>(cfg.latency(1)) & begTrialLatency<(cfg.latency(2));
hasWindow     = true(nTrials,1);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
  startsLater    = begTrialLatency > (cfg.latency(1) + 0.0001);
  endsEarlier    = endTrialLatency < (cfg.latency(2) - 0.0001);
  hasWindow      = ~(startsLater | endsEarlier); % it should not start later or end earlier
  trialDur       = ones(nTrials,1)*(cfg.latency(2)-cfg.latency(1));
elseif strcmp(cfg.vartriallen,'yes')
  winBeg      = max([begTrialLatency(:) cfg.latency(1)*ones(nTrials,1)],[],2); 
  winEnd      = min([endTrialLatency(:) cfg.latency(2)*ones(nTrials,1)],[],2);
  trialDur    = winEnd-winBeg; % the effective trial duration  
end
cfg.trials         = cfg.trials(overlaps(:) & hasWindow(:));
trialDur           = trialDur(overlaps(:) & hasWindow(:)); % select the durations, we will need this later
nTrials            = length(cfg.trials);
if isempty(cfg.trials), warning('No trials were selected in the end'); end

% preallocate before computing the psth
keepTrials = strcmp(cfg.keeptrials,'yes');
if  keepTrials, singleTrials = NaN(nTrials,nUnits); end % preallocate single trials with NaNs
[s,ss]   = deal(zeros(nUnits,1));
dof      = nTrials*ones(nUnits,1); % compute degrees of freedom
for iUnit = 1:nUnits
  unitIndx   = spikesel(iUnit);
  ts         = spike.time{unitIndx}(:); % get the times
  latencySel = ts>=cfg.latency(1) & ts<=cfg.latency(2); % select timestamps within latency
  trialSel   = ismember(spike.trial{unitIndx},cfg.trials);
  trialNums  = spike.trial{unitIndx}(latencySel(:) & trialSel(:));
  
  % use the fact that trial numbers are integers >=1 apart, so we can use histc
  trialBins   = sort([cfg.trials-0.5; cfg.trials+0.5]);
  trialRate   = histc(trialNums(:),trialBins);
  trialRate   = trialRate(1:2:end-1); % the uneven bins correspond to the trial integers
  if isempty(trialRate), trialRate = zeros(nTrials,1); end
  
  % convert to firing rates if requested
  if strcmp(cfg.outputunit,'rate') ,trialRate = trialRate(:)./trialDur(:); end
  
  % store and compute sum, just fill the single trials up with nans
  if keepTrials, singleTrials(:,iUnit) = trialRate(:); end
  s(iUnit)  = sum(trialRate);
  ss(iUnit) = sum(trialRate.^2);
end

% compute the average rate
rate.avg       = s ./ dof;
rate.var       = (ss - s.^2./dof)./(dof); % since sumrate.^2 ./ dof = dof .* (sumrate/dof).^2
rate.var(dof<2)  = 0;

% gather the rest of the results
rate.dof       = dof;
rate.label     = spike.label(spikesel);
rate.dimord    = 'chan_time';
rate.time      = mean(cfg.latency);
if (strcmp(cfg.keeptrials,'yes'))
  rate.trial = singleTrials;
  rate.dimord = 'rpt_chan_time';
end
if isfield(spike, 'trialinfo'), rate.trialinfo = spike.trialinfo(cfg.trials,:); end
  
% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   spike
ft_postamble provenance rate
ft_postamble history    rate


%%%%%%%%% SUB FUNCTIONS %%%%%%%%%
function [cfg] = latencyselection(cfg,begTrialLatency,endTrialLatency)

if strcmp(cfg.latency,'minperiod')
  cfg.latency = [max(begTrialLatency) min(endTrialLatency)];
elseif strcmp(cfg.latency,'maxperiod')
  cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
  cfg.latency = [min(begTrialLatency) 0];
elseif strcmp(cfg.latency,'poststim')
  cfg.latency = [0 max(endTrialLatency)];
end

function [cfg] = trialselection(cfg,spike)

% get the number of trials or change DATA according to cfg.trials
nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, error('maximum trial number in cfg.trials should not exceed size(spike.trialtime,1)'); end
if isempty(cfg.trials), error('no trials were selected by you'); end

                   f
end

for iTrial = 1:nTrials 
  origTrial = cfg.trials(iTrial);
  if  ~ (allStartEarlier && allEndLater) % select bins and count dof + 1
    binSel = begTrialLatency(iTrial)<=bins(1:end-1) & endTrialLatency(iTrial)>=bins(2:end);
    dof(binSel)      = dof(binSel) + 1;
  else
    binSel           = 1:(nBins-1); % always deselect the last bin
  end
  
  for iUnit = 1:nUnits
    unitIndx      = spikesel(iUnit); % select the unit
    spikesInTrial = (spike.trial{unitIndx}==origTrial); % get spikes in trial
    spikeTimes    = spike.time{unitIndx}(spikesInTrial);
    
    % compute the psth
    trialPsth   = histc(spikeTimes(:), bins); % we deselect the last bin per default
    trialPsth   = trialPsth(:)'; % force into row vec
    if isempty(trialPsth), trialPsth = zeros(1,length(bins)); end
    
    % convert to firing rates if requested, with spikecount do nothing
    if strcmp(cfg.outputunit,'rate'), 
      trialPsth = trialPsth/cfg.binsize; 
    elseif strcmp(cfg.outputunit,'proportion'), 
      trialPsth = trialPsth./nansum(trialPsth); 
    end
    
    % compute the sum and the sum of squares for the var and the mean on the fly
    s(iUnit,binSel)  = s(iUnit,binSel)   + trialPsth(binSel); % last bin is single value
    ss(iUnit,binSel) = ss(iUnit,binSel)  + trialPsth(binSel).^2;
        
    if strcmp(cfg.keeptrials,'yes'),
      singleTrials(iTrial,iUnit,binSel) = trialPsth(binSel);
    end
  end
end

% get the results
dof            = dof(ones(nUnits,1),:);
psth.avg       = s ./ dof;
psth.var       = (ss - s.^2./dof)./(dof-1); % since sumPsth.^2 ./ dof = dof .* (sumPsth/dof).^2
psth.dof       = dof; % combined with psth.var we can get SEM
psth.fsample   = 1/(cfg.binsize);   % might be more compatible with feeding psth in other funcs
psth.time      = bins(1:end-1) + 0.5*cfg.binsize; 
psth.label     = spike.label(spikesel);
if (strcmp(cfg.keeptrials,'yes'))
  psth.trial  = singleTrials;
  psth.dimord = 'rpt_chan_time';
else
  psth.dimord = 'chan_time';
end
if isfield(spike,'sampleinfo'), psth.sampleinfo = spike.sampleinfo(cfg.trials,:); end
if isfield(spike,'trialinfo'),  psth.trialinfo  = spike.trialinfo(cfg.trials,:);  end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   spike
ft_postamble provenance psth
ft_postamble history    psth


%%%%%%%%% SUB FUNCTIONS %%%%%%%%%
function [cfg] = latencyselection(cfg,begTrialLatency,endTrialLatency)

if strcmp(cfg.latency,'minperiod')
  cfg.latency = [max(begTrialLatency) min(endTrialLatency)];
elseif strcmp(cfg.latency,'maxperiod')
  cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
  cfg.latency = [min(begTrialLatency) 0];
elseif strcmp(cfg.latency,'poststim')
  cfg.latency = [0 max(endTrialLatency)];
end
% check whether the time window fits with the data
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency);
  warning('Correcting begin latency because it is before all trial beginnings');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('Correcting end latency because it is after all trial ends');
end


function [cfg] = trialselection(cfg,spike)

nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, 
  error('maximum trial number in cfg.trials should not exceed number of rows of spike.trialtime'); 
end
if isempty(cfg.trials), error('No trials were selected by you, rien ne va plus'); end



                                                                                                                                                                                                                                                                                                                                                                                                                        h
    avgHdl  = feval(cfg.topplotfunc,dataX,dataY');
    
    % set the color for every unit
    for iUnit = 1:nUnits
      set(avgHdl(iUnit),'Color',cfg.cmapneurons(iUnit,:));
    end
    
    % ensure that the data limits are increasing
    mnmax = [nanmin(dataY(:))-eps nanmax(dataY(:))+eps];
    
    % plot the errorbars
    if ~strcmp(cfg.errorbars,'no')
      
      % check if the right error information is there
      if ~isfield(timelock, 'var') || ~isfield(timelock, 'dof') 
        error('timelock should contain field .var and .dof for errorbars'); 
      end
      
      % select the degrees of freedom
      df = timelock.dof(unitIndx,binSel);
      
      % plot the different error bars
      if strcmp(cfg.errorbars, 'sem')
        err = sqrt(timelock.var(unitIndx,binSel)./df);
      elseif strcmp(cfg.errorbars, 'std')
        err = sqrt(timelock.var(unitIndx,binSel));
      elseif strcmp(cfg.errorbars, 'var')
        err = timelock.var(unitIndx,binSel);
      elseif strcmp(cfg.errorbars, 'conf95%')
        tCrit = tinv(0.975,df);
        err   = tCrit.*sqrt(timelock.var(unitIndx,binSel)./df); 
      end
      
      % ensure division by zero does not count, should not happen however
      err(~isfinite(err)) = NaN;
      
      % make a polygon of the error bars that allows easy filling in adobe
      for iUnit = 1:nUnits
        upb  = dataY(iUnit,:) + err(iUnit,:);
        lowb = dataY(iUnit,:) - err(iUnit,:);
        sl   = ~isnan(upb);
        [X,Y] = polygon(dataX(sl),upb(sl)+0.0001,lowb(sl)-0.0001);
        hold on
        hd = plot(X,Y,'--');
        set(hd,'Color', cfg.cmapneurons(iUnit,:));
        hold on
      end
      mnmax = [nanmin(dataY(:) - err(:)) nanmax(dataY(:) + err(:))];
    end
    
    % set the y limits explicitly
    ylim = mnmax;
    d    = ylim(2)-ylim(1);
    yl = ylim(1)-0.05*d;
    if ylim(1)>0
      if yl<0, yl = 0; end
    end
    ylim(1) = yl;

    yl = ylim(2)+0.05*d;
    if ylim(2)<0
      if yl>0, yl = 0; end
    end
    ylim(2) = yl;
  if ylim(2) == ylim(1) %if the plot is empty
    ylim(2)=ylim(1)+1; %be nice to set
  end
    set(gca,'YLim', ylim)
  end    
  
  % store the handle
  cfg.hdl.psth = avgHdl;
  
  % modify the axes
  set(ax(2),'YAxisLocation', 'right') % swap y axis location
  set(ax(2),'XTickLabel', {}) % remove ticks and labels for x axis
  set(ax(2),'XTick', [])
  set(ax(2), 'Box', 'off') % turn the box off
  
  % change the axis settings
  try
    ylabel(timelock.cfg.outputunit)
  catch
    warning('unit of the y-axis is unknown'); 
    ylabel('Signal intensity')
  end
end

% set the limits for the axis
set(ax,'XLim', [cfg.latency])
if nTrialsShown==0; nTrialsShown = 1; end %
set(ax(1), 'YLim', [0.5 nTrialsShown+0.5]); % number of trials
set(ax,'TickDir','out') % put the tickmarks outside

% now link the axes, constrain zooming and keep ticks intact
limX       = [cfg.latency];
limY       = get(ax,'YLim');
if ~iscell(limY), limY = {limY}; end

% constrain the zooming and zoom psth together with the jpsth, remove
% ticklabels jpsth
if strcmp(cfg.interactive,'yes')
  set(zoom,'ActionPostCallback',{@mypostcallback,ax,limX,limY});
  set(pan,'ActionPostCallback',{@mypostcallback,ax,limX,limY});
end

% pass positions and axis handles so downstream
% functions can respect the positions set here.
cfg.pos.posRaster = posRaster;
cfg.hdl.axRaster = ax(1);
if doTopData
 cfg.pos.posTopPlot = posTopPlot;
 cfg.hdl.axTopPlot = ax(2);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous spike
ft_postamble provenance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = mypostcallback(fig,evd,ax,lim,ylim)
currentAxes = evd.Axes;
indx = find(currentAxes==ax);

% keep the y axis within boundaries
origLim = ylim{indx};
ylim = get(ax(indx), 'YLim');
if origLim(1)>ylim(1), ylim(1) = origLim(1); end
if origLim(2)<ylim(2), ylim(2) = origLim(2); end
set(ax(indx), 'YLim', ylim)

% keep the x axis within boundaries and reset for both
xlim = get(ax(indx), 'XLim');
if lim(1)>xlim(1), xlim(1) = lim(1); end
if lim(2)<xlim(2), xlim(2) = lim(2); end
set(ax,'XLim', xlim)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [colorspace] = colormap_cgbprb(N)
% COLORMAP_SEPARATION returns a color map with well separated colors that does not include
% white, yellow or orange since these are not well visible in white plots.
%
% Inputs: N is the number of colors desired
% Outputs: COLORSPACE is an N-by-3 colorspace.
% Chosen colors are green, cyan, blue, purple, red and black
% Script was adapted from varycolors.m from FEX
% Copyright (c) 2009, Martin Vinck

if nargin<1
  error('specify the number of colors that you want')
end

% create a well visible color space
%        green     cyan  blue    purple    red     black
colors = [[0 1 0]; [0 1 1] ; [0 0 1] ; [1 0 1]; [1 0 0]; [0 0 0]];
order  = [1 5 3 2 4 6]; % our preference order when plotting different lines
rm{2} = [2:5];
rm{3} = [2 4 5];
rm{4} = [2 4];
rm{5} = [4];

if N==1
  colorspace = colors(end,:);
elseif N>1&&N<6
  colors(rm{N},:) = [];
  colorspace = colors;
else
  n      = floor(N/5)*ones(1,5);
  rest   = mod(N,5);
  order  = [1 5 3 2 4];
  % if we have some rest, add them starting oppositly
  n(order(1:rest)) = n(order(1:rest)) + 1;
  
  colorspace = [];
  for iColor = 1:5
    corStart  = colors(iColor,:);
    corEnd    = colors(iColor+1,:);
    dim       = corStart~=corEnd;
    subColors = corStart(ones(n(iColor)+1,1),:);
    if iColor>1
      subColors(:,dim) = linspace(corStart(dim),corEnd(dim),n(iColor)+1);
      subColors(1,:) = [];
    else
      subColors(1:end-1,dim) =   linspace(corStart(dim),corEnd(dim),n(iColor));
      subColors(end,:) = [];
    end
    colorspace = [colorspace;subColors];
  end
end

function [X,Y] = polygon(x,sm1,sm2,multiplier)

x = x(:)';
if nargin < 4
    multiplier = 1;
end
X = [x x(end:-1:1) x(1)];
up = sm1(:)';
down = sm2(:)';
Y = [down up(end:-1:1) down(1)];
                                                                                                                                                                                                                                                                                                                                                                                                                                                  trialselection(cfg,spike)

% get the number of trials or change DATA according to cfg.trials
nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials, 'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, warning('maximum trial number in cfg.trials should not exceed length of DATA.trial')
end
if isempty(cfg.trials), error('No trials were selected');
end

function m = nansum(x,dim)

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x);
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim);
end






                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     av) + 1i*sum(ones(1,timwinSamples).*sinwav); % fft of unit direct current * taper
    spctrm  = spctrm + (conv_fftbased(dat(:),wavelet) - (beta0*fftDC + beta1.*fftRamp))/(numsmp/2);           
                       % fft                       % mean            %linear ramp      % make magnitude invariant to window length                             
end
spctrm = spctrm./nTapers; % normalize by number of tapers
spctrm = spctrm.*exp(-1i*phaseCor);

% set the part of the spectrum without a valid phase to NaNs
n = (timwinSamples-1)/2;
spctrm(1:n) = NaN; 
spctrm(end-n+1:end) = NaN;

function [taper] = double_dpss(a, b, varargin)
taper = dpss(double(a), double(b), varargin{:});


function [spikelabel, eeglabel] = detectspikechan(data)

maxRate = 1000; % default on what we still consider a neuronal signal

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    hasAllInts    = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:) == round(data.trial{i}(j,:)));
    hasAllPosInts = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:)>=0);
    fr            = nansum(data.trial{i}(j,:),2) ./ (data.time{i}(end)-data.time{i}(1));    
    spikechan(j) = spikechan(j) + double(hasAllInts & hasAllPosInts & fr<=maxRate);
  end
end
spikechan = (spikechan==ntrial);

spikelabel = data.label(spikechan);
eeglabel   = data.label(~spikechan);

% CONVOLUTION: FFT BASED IMPLEMENTATION
function c = conv_fftbased(a, b)

P = numel(a);
Q = numel(b);
L = P + Q - 1;
K = 2^nextpow2(L);

c = ifft(fft(a, K) .* fft(b, K));
c = c(1:L);
toRm1 = [1:(Q-1)/2];
toRm2 = [(1 + length(c) - (Q-1)/2) : length(c)];
toRm  = [toRm1(:); toRm2(:)];
c(toRm) = [];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          