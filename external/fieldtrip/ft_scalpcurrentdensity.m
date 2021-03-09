function [scd] = ft_scalpcurrentdensity(cfg, data)

% FT_SCALPCURRENTDENSITY computes an estimate of the SCD using the
% second-order derivative (the surface Laplacian) of the EEG potential
% distribution
%
% The relation between the surface Laplacian and the SCD is explained
% in more detail on http://tinyurl.com/ptovowl.
%
% Use as
%   [data] = ft_scalpcurrentdensity(cfg, data)
% or
%   [timelock] = ft_scalpcurrentdensity(cfg, timelock)
% where the input data is obtained from FT_PREPROCESSING or from
% FT_TIMELOCKANALYSIS. The output data has the same format as the input
% and can be used in combination with most other FieldTrip functions
% such as FT_FREQNALYSIS or FT_TOPOPLOTER.
%
% The configuration should contain
%   cfg.method       = 'finite' for finite-difference method or
%                      'spline' for spherical spline method
%                      'hjorth' for Hjorth approximation method
%   cfg.elec         = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.trials       = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.feedback     = string, 'no', 'text', 'textbar', 'gui' (default = 'text')
%
% The finite method require the following
%   cfg.conductivity = conductivity of the skin (default = 0.33 S/m)
%
% The spline and finite method require the following
%   cfg.conductivity = conductivity of the skin (default = 0.33 S/m)
%   cfg.lambda       = regularization parameter (default = 1e-05)
%   cfg.order        = order of the splines (default = 4)
%   cfg.degree       = degree of legendre polynomials (default for
%                       <=32 electrodes  = 9,
%                       <=64 electrodes  = 14,
%                       <=128 electrodes = 20,
%                       else             = 32
%
% The hjorth method requires the following
%   cfg.neighbours   = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%
% For the spline method you can specify the following
%   cfg.badchannel      = cell-array, see FT_CHANNELSELECTION for details (default = [])
%
% Note that the skin conductivity, electrode dimensions and the potential
% all have to be expressed in the same SI units, otherwise the units of
% the SCD values are not scaled correctly. The spatial distribution still
% will be correct.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% The 'finite' method implements
%   TF Oostendorp, A van Oosterom; The surface Laplacian of the potential:
%   theory and application. IEEE Trans Biomed Eng, 43(4): 394-405, 1996.
%   G Huiskamp; Difference formulas for the surface Laplacian on a
%   triangulated sphere. Journal of Computational Physics, 2(95): 477-496,
%   1991.
%
% The 'spline' method implements
%   F. Perrin, J. Pernier, O. Bertrand, and J. F. Echallier.
%   Spherical splines for scalp potential and curernt density mapping.
%   Electroencephalogr Clin Neurophysiol, 72:184-187, 1989
% including their corrections in
%   F. Perrin, J. Pernier, O. Bertrand, and J. F. Echallier.
%   Corrigenda: EEG 02274, Electroencephalography and Clinical
%   Neurophysiology 76:565.
%
% The 'hjorth' method implements
%   B. Hjort; An on-line transformation of EEG scalp potentials into
%   orthogonal source derivation. Electroencephalography and Clinical
%   Neurophysiology 39:526-530, 1975.
%
% See also FT_PREPROCESSING, FT_TIMELOCKANALYSIS, FT_FREQNALYSIS, FT_TOPOPLOTER.

% Copyright (C) 2004-2012, Robert Oostenveld
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
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.method          = ft_getopt(cfg, 'method',       'spline');
cfg.conductivity    = ft_getopt(cfg, 'conductivity', 0.33); % in S/m
cfg.trials          = ft_getopt(cfg, 'trials',       'all', 1);
cfg.feedback        = ft_getopt(cfg, 'feedback',     'text');
cfg.badchannel      = ft_getopt(cfg, 'badchannel',     {});

switch cfg.method
  case 'hjorth'
    cfg = ft_checkconfig(cfg, 'required', {'neighbours'});
  case 'spline'
    cfg.lambda  = ft_getopt(cfg, 'lambda', 1e-5);
    cfg.order   = ft_getopt(cfg, 'order', 4);
    cfg.degree  = ft_getopt(cfg, 'degree', []);

    if isempty(cfg.degree) % determines degree of Legendre polynomials bases on number of electrodes
      nchan = numel(data.label);
      if nchan<=32
        cfg.degree = 9;
      elseif nchan<=64
        cfg.degree = 14;
      elseif nchan<=128
        cfg.degree = 20;
      else
        cfg.degree = 32;
      end
    end
  otherwise
    cfg = ft_checkconfig(cfg); % perform a simple consistency check
end

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'ismeg', []);

% get the electrode positions
tmpcfg = cfg;
tmpcfg.senstype = 'EEG';
elec = ft_fetch_sens(tmpcfg, data);

% select channels and trials of interest
tmpcfg = keepfields(cfg, {'trials', 'showcallinfo'});
tmpcfg.channel = elec.label;
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

Ntrials = numel(data.trial);

if isempty(cfg.badchannel)
  % check if the first sample of the first trial contains NaNs; if so treat it as a bad channel
  cfg.badchannel = ft_channelselection(find(isnan(data.trial{1}(:,1))), data.label);
  if ~isempty(cfg.badchannel)
    ft_info('detected channel %s as bad\n', cfg.badchannel{:});
  end
end

% match the order of the data channels with the channel positions, order them according to the data
[datindx, elecindx] = match_str(data.label, elec.label);
[goodindx, tmp]     = match_str(data.label, setdiff(data.label, cfg.badchannel));

allchanpos = elec.chanpos(elecindx,:);    % the position of all channels, ordered according to the data
goodchanpos = allchanpos(goodindx,:);     % the position of good channels
anychanpos  = [0 0 1];                    % random channel, will be used in case there are no bad channels

% compute SCD for each trial
if strcmp(cfg.method, 'spline')
  ft_progress('init', cfg.feedback, 'computing SCD for trial...')
  for trlop=1:Ntrials
    % do compute interpolation
    ft_progress(trlop/Ntrials, 'computing SCD for trial %d of %d', trlop, Ntrials);
    if ~isempty(cfg.badchannel)
      % compute scd for all channels, also for the bad ones
      fprintf('computing scd also at locations of bad channels');
      [V2, L2, L1] = splint(goodchanpos, data.trial{trlop}(goodindx,:), allchanpos, cfg.order, cfg.degree, cfg.lambda);
      scd.trial{trlop} = L2;
    else
      % just compute scd for input channels, specify arbitrary channel to-be-discarded for interpolation to save computation time
      [V2, L2, L1] = splint(goodchanpos, data.trial{trlop}(goodindx,:), anychanpos, cfg.order, cfg.degree, cfg.lambda);
      scd.trial{trlop} = L1;
    end
  end
  ft_progress('close');

elseif strcmp(cfg.method, 'finite')
  if ~isempty(cfg.badchannel)
    ft_error('the method "%s" does not support the specification of bad channels', cfg.method);
  end
  % the finite difference approach requires a triangulation
  prj = elproj(allchanpos);
  tri = delaunay(prj(:,1), prj(:,2));
  % the new electrode montage only needs to be computed once for all trials
  montage.tra = lapcal(allchanpos, tri);
  montage.labelold = data.label;
  montage.labelnew = data.label;
  % apply the montage to the data, also update the electrode definition
  scd  = ft_apply_montage(data, montage);
  elec = ft_apply_montage(elec, montage);

elseif strcmp(cfg.method, 'hjorth')
  if ~isempty(cfg.badchannel)
    ft_error('the method "%s" does not support the specification of bad channels', cfg.method);
  end
  % convert the neighbourhood structure into a montage
  labelnew = {};
  labelold = {};
  for i=1:length(cfg.neighbours)
    labelnew = cat(2, labelnew, cfg.neighbours(i).label);
    labelold = cat(2, labelold, cfg.neighbours(i).neighblabel(:)');
  end
  labelold = cat(2, labelnew, labelold);
  labelold = unique(labelold);
  tra = zeros(length(labelnew), length(labelold));
  for i=1:length(cfg.neighbours)
    thischan   = match_str(labelold, cfg.neighbours(i).label);
    thisneighb = match_str(labelold, cfg.neighbours(i).neighblabel);
    tra(i, thischan) = 1;
    tra(i, thisneighb) = -1/length(thisneighb);
  end
  % combine it in a montage
  montage.tra = tra;
  montage.labelold = labelold;
  montage.labelnew = labelnew;
  % apply the montage to the data, also update the electrode definition
  scd  = ft_apply_montage(data, montage);
  elec = ft_apply_montage(elec, montage);

else
  ft_error('unknown method "%s"', cfg.method);
end

if strcmp(cfg.method, 'spline') || strcmp(cfg.method, 'finite')
  % correct the units
  ft_warning('trying to correct the units, assuming uV and mm');
  for trlop=1:Ntrials
    % The surface laplacian is proportional to potential divided by squared distance which means that, if
    % - input potential is in uV, which is 10^6 too large
    % - units of electrode positions are in mm, which is 10^3 too large
    % these two cancel out against each other. Hence the computed laplacian
    % is in SI units (MKS).
    scd.trial{trlop} = cfg.conductivity * -1 * scd.trial{trlop};
  end
  fprintf('output surface laplacian is in V/m^2\n');
else
  fprintf('output Hjorth filtered potential is in uV\n');
end

% collect the results
scd.elec    = elec;
scd.time    = data.time;
scd.label   = data.label;
scd.fsample = 1/mean(diff(data.time{1}));
if isfield(data, 'sampleinfo')
  scd.sampleinfo = data.sampleinfo;
end
if isfield(data, 'trialinfo')
  scd.trialinfo = data.trialinfo;
end

% convert back to input type if necessary
switch dtype
  case 'timelock'
    scd = ft_checkdata(scd, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = scd;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
                                                                                                                                                                                                                                                                                                                                           the specified new time axes for each trial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ntr = length(data.trial);
  
  ft_progress('init', cfg.feedback, 'resampling data');
  for itr = 1:ntr
    ft_progress(itr/ntr, 'resampling data in trial %d from %d\n', itr, ntr);
    
    olddat = data.trial{itr};
    oldtim = data.time{itr};
    
    % detrending is in general not recommended
    if istrue(cfg.detrend)
      if ~strcmp(cfg.baselinewindow, 'all')
        olddat = ft_preproc_detrend(olddat, nearest(oldtim, cfg.baselinewindow(1)), nearest(oldtim, cfg.baselinewindow(2)));
      else
        olddat = ft_preproc_detrend(olddat);
      end
    end
    
    % always remove the mean to avoid edge effects when there's a strong offset, the cfg.demean option is dealt with below
    if ~strcmp(cfg.baselinewindow, 'all')
      [olddat, bsl] = ft_preproc_baselinecorrect(olddat, nearest(oldtim, cfg.baselinewindow(1)), nearest(oldtim, cfg.baselinewindow(2)));
    else
      [olddat, bsl] = ft_preproc_baselinecorrect(olddat);
    end
    
    % perform the resampling
    newtim = cfg.time{itr};
    if length(oldtim)>1
      newdat = interp1(oldtim', olddat', newtim', cfg.method)';
    else
      newdat = repmat(olddat, [1 numel(newtim)]);
    end
    
    % add back the mean
    if ~strcmp(cfg.demean, 'yes')
      nsmp   = size(newdat, 2);
      newdat = newdat + bsl(:,ones(1,nsmp));
    end
    
    data.trial{itr} = newdat;
    data.time{itr}  = newtim;
    
  end % for itr
  ft_progress('close');
  
  % specify the new sampling frequency in the output
  t1 = cfg.time{1}(1);
  t2 = cfg.time{1}(2);
  data.fsample = 1/(t2-t1);
  
end % if usefsample or usetime

ft_info('original sampling rate = %d Hz\nnew sampling rate = %d Hz\n', cfg.origfs, data.fsample);

% convert back to input type if necessary
switch convert
  case 'timelock'
    data = ft_checkdata(data, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that decimates along the columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = my_decimate(x, varargin)
[n, m] = size(x);
% decimate the first column
y = decimate(x(:,1), varargin{:});
if m>1
  % increase the size of the output matrix
  y(:,m) = 0;
  % decimate the subsequent columns
  for i=2:m
    y(:,i) = decimate(x(:,i), varargin{:});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that does a block-wise average along the columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = my_mean(x, r)
[n, m] = size(x);
n = n - mod(n,r);
x = x(1:n,:);
x = reshape(x, [r n/r m]);
y = mean(x, 1);
y = reshape(y, [n/r m]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that does a block-wise median along the columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = my_median(x, r)
[n, m] = size(x);
n = n - mod(n,r);
x = x(1:n,:);
x = reshape(x, [r n/r m]);
y = median(x, 1);
y = reshape(y, [n/r m]);
                                                                                                                                                                                    fctdef, 'writerej') && ~isempty(cfg.artfctdef.writerej)
  fid = fopen(cfg.artfctdef.writerej, 'w');
  if fid<0
    ft_error('could not open rejection file for writing');
  else
    % determine the begin and end of each rejected period (in samples)
    rejectonset = find(diff([0 rejectall])== 1);
    rejectofset = find(diff([rejectall 0])==-1);
    % determine the begin and end of each rejected period (in seconds)
    rejectonset = (rejectonset-1)/hdr.Fs;
    rejectofset = (rejectofset-1)/hdr.Fs;
    for rejlop=1:length(rejectonset)
      fprintf(fid, '%f-%f\n', rejectonset(rejlop), rejectofset(rejlop));
    end
    fclose(fid);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the trials that (partially) coincide with a rejection mark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(cfg.artfctdef.reject, {'partial', 'complete', 'nan', 'value'}))
  trialok = [];
  
  count_complete_reject = 0;
  count_partial_reject  = 0;
  count_nan             = 0;
  count_value           = 0;
  count_outsidecrit     = 0;
  
  trlCompletelyRemovedInd = [];
  trlPartiallyRemovedInd  = [];
  
  for trial=1:size(trl,1)
    % cut out the part of the rejection-axis corresponding with this trial
    rejecttrial = rejectall(trl(trial,1):trl(trial,2));
    
    if all(not(rejecttrial))
      % the whole trial is good
      trialok = [trialok; trl(trial,:)];
      
    elseif all(rejecttrial) && strcmp(cfg.artfctdef.reject, 'nan')
      % the whole trial is bad, but it is requested to be replaced with nans
      data.trial{trial}(:,rejecttrial) = nan;
      count_nan = count_nan + 1;
      trialok = [trialok; trl(trial,:)]; % Mark the trial as good as nothing will be removed
      
    elseif all(rejecttrial) && strcmp(cfg.artfctdef.reject, 'value')
      % the whole trial is bad, but it is requested to be replaced with a specific value
      data.trial{trial}(:,rejecttrial) = cfg.artfctdef.value;
      count_value = count_value + 1;
      trialok = [trialok; trl(trial,:)]; % Mark the trial as good as nothing will be removed

    elseif all(rejecttrial)
      % the whole trial is bad
      count_complete_reject = count_complete_reject + 1;
      trlCompletelyRemovedInd = [trlCompletelyRemovedInd trial];
      
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'complete')
      % some part of the trial is bad, check if within crittoilim?
      if (checkCritToi)
        critInd = (data.time{trial} >= cfg.artfctdef.crittoilim(trial,1) ...
          & data.time{trial} <= cfg.artfctdef.crittoilim(trial,2));
        if (any(critInd & rejecttrial))
          count_complete_reject = count_complete_reject + 1;
          trlCompletelyRemovedInd = [trlCompletelyRemovedInd trial];
          continue;
        else
          trialok = [trialok; trl(trial,:)];
          count_outsidecrit = count_outsidecrit + 1;
        end
      else % no crittoilim checking required
        count_complete_reject = count_complete_reject + 1;
        trlCompletelyRemovedInd = [trlCompletelyRemovedInd trial];
        continue;
      end
      
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'partial')
      % some part of the trial is bad, reject only the bad part
      trialnew = [];
      rejecttrial = [0 not(rejecttrial) 0];
      % the first column is the begin sample, the second the end sample and the third is the offset
      trialnew(:,1) = find(diff(rejecttrial(:))== 1) + trl(trial,1) - 1;
      trialnew(:,2) = find(diff(rejecttrial(:))==-1) + trl(trial,1) - 2;
      trialnew(:,3) = find(diff(rejecttrial(:))== 1) - 1 + trl(trial,3);
      % some people use additional columns in the trl matrix to store trigger values and/or reaction times
      % these should remain linked to the original trial, i.e. they should be copied for each new fragment
      for i=4:size(trl,2)
        trialnew(:,i) = trl(trial,i);
      end
      minacceptnumsmp = round(cfg.artfctdef.minaccepttim .* hdr.Fs);
      trialnew((trialnew(:,2)-trialnew(:,1))<minacceptnumsmp,:) = [];
      count_partial_reject = count_partial_reject + 1;
      trialok = [trialok; trialnew];
      trlPartiallyRemovedInd = [trlPartiallyRemovedInd trial];
      
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'nan')
      % Some part of the trial is bad, replace bad part with nans
      data.trial{trial}(:,rejecttrial) = nan;
      count_nan = count_nan + 1;
      trialok = [trialok; trl(trial,:)]; % Mark the trial as good as nothing will be removed
   
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'value')
      % Some part of the trial is bad, replace bad part with specified value
      data.trial{trial}(:,rejecttrial) = cfg.artfctdef.value;
      count_value = count_value + 1;
      trialok = [trialok; trl(trial,:)]; % Mark the trial as good as nothing will be removed

    end
  end % for each trial
  
  fprintf('rejected  %3d trials completely\n', count_complete_reject);
  fprintf('rejected  %3d trials partially\n', count_partial_reject);
  fprintf('filled parts of  %3d trials with nans\n', count_nan);
  fprintf('filled parts of  %3d trials with the specified value\n', count_value);
  if (checkCritToi)
    fprintf('retained  %3d trials with artifacts outside critical window\n', count_outsidecrit);
  end
  fprintf('resulting %3d trials\n', size(trialok,1));
  cfg.trlold = trl;      % return the original trial definition in the configuration
  cfg.trl    = trialok;  % return the cleaned trial definition in the configuration
  
  if strcmp(cfg.artfctdef.feedback, 'yes')
    fprintf('the following trials were completely removed: ');
    for k = trlCompletelyRemovedInd
      fprintf('%d ', k);
    end
    fprintf('\nthe following trials were partially removed: ');
    for k = trlPartiallyRemovedInd
      fprintf('%d ', k);
    end
    fprintf('\n');
  end
  
else
  fprintf('not rejecting any data, only marking the artifacts\n');
end

if isempty(cfg.trl)
  ft_error('No trials left after artifact rejection.')
else
  if hasdata && ~strcmp(cfg.artfctdef.reject, 'nan') % Skip this step to avoid removing parts that should be filled with nans
    % apply the updated trial definition on the data
    tmpcfg      = keepfields(cfg, {'trl', 'showcallinfo'});
    data        = removefields(data, {'trialinfo'});
    data        = ft_redefinetrial(tmpcfg, data);
    % restore the provenance information
    [cfg, data] = rollback_provenance(cfg, data);
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
if hasdata
  ft_postamble previous data
  ft_postamble history data
  ft_postamble savevar data
  % return the data, the output variable is called cfg instead of data
  cfg = data;
end
                                                                                                                                                                                                                                                                                                                                                    'Parent',h.MainFigure,...
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