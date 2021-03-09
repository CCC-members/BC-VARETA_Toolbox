function [source] = ft_source2full(source)

% FT_SOURCE2FULL recreates the grid locations outside the brain in the source 
% reconstruction, so that the source volume again describes the full grid.
% This undoes the memory savings that can be achieved using FT_SOURCE2SPARSE
% and makes it possible again to plot  the source volume and save it to an
% external file.
%
% Use as
%   [source] = ft_source2full(source)
%
% See also FT_SOURCE2SPARSE

% Copyright (C) 2004, Robert Oostenveld
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

ft_defaults

if ~isfield(source, 'inside')  || ...
   ~isfield(source, 'outside') || ...
   ~isfield(source, 'dim')
  ft_error('one of the required fields is missing in the source structure');
end

if ~isfield(source, 'pos') && (~isfield(source, 'xgrid') || ~isfield(source, 'ygrid') || ...
                               ~isfield(source, 'zgrid'))
  ft_error('the input data needs at least a ''pos'' field, or ''x/y/zgrid''');
end

if isfield(source, 'xgrid'),
  xgrid = source.xgrid;
  ygrid = source.ygrid;
  zgrid = source.zgrid;
  sparsepos = source.pos;
  
  % recreate the positions of the dipole grid
  [X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
  pos = [X(:) Y(:) Z(:)];
else
  %FIXME this assumes that the voxel data are ordered as if in a regularly spaced 3D grid,
  %but with only the inside voxels present
  ft_warning('assuming the voxel data to be ordered as if in a regularly spaced 3D grid');
  xgrid = 1:source.dim(1);
  ygrid = 1:source.dim(2);
  zgrid = 1:source.dim(3);

  %establish a homogeneous transformation matrix from voxels to headspace based on the sparse positions
  sparsepos = source.pos;
  ok  = 0;
  cnt = 0;
  while ok==0,
    cnt  = cnt+1;
    dpos = sparsepos - sparsepos(cnt*ones(size(sparsepos,1),1),:);
    [srt, indx] = sort(sum(dpos.^2,2));
    srt    = dpos(indx,:);
    tmpsrt = abs(srt(2:7,:));
    csrt   = tmpsrt*tmpsrt';
    sel    = find(sum(csrt==0)>=2);
    if numel(sel)>=3, 
      ok = 1;
    end
  end  
  tmppos  = sparsepos(indx([1 sel(:)'+1]),:);
  tmpdpos = dpos(indx([1 sel(:)'+1]),:);
 
  % FIXME the following is a bit experimental and not fully tested yet it works in general case rotation
  M         = pinv(tmpdpos(2:4,:));
  
  % get rotation such that maxima are on diagonal and positive
  m(1) = find(M(1,:)==max(abs(M(1,:))));
  m(2) = find(M(2,:)==max(abs(M(2,:))));
  m(3) = find(M(3,:)==max(abs(M(3,:))));
  [srt, indx] = sort(m);
  M    = M(indx,:);
  M    = M*diag(sign(diag(M)));
  sparsepos = sparsepos*M;
  
  % translation
  T         = -min(sparsepos,[],1)+1;
  sparsepos = sparsepos + T(ones(size(sparsepos,1),1), :);  

  % recreate the positions of the dipole grid
  [X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
  pos = [X(:) Y(:) Z(:)];
  pos = ft_warp_apply(inv([M T(:);0 0 0 1]), pos);
end

Nsparse = length(source.inside);
siz     = source.dim;
Nfull   = prod(siz);

% determine the size that each slice takes in memory
sx = 1;
sy = siz(1);
sz = siz(1) * siz(2);

if isfield(source, 'inside') && isfield(source, 'outside') && size(source.pos,1)==Nfull
  % it contains all source positions
  inside = source.inside;
  outside = source.outside;
else
  % it only contains the inside source positions, which are all inside the brain
  % reconstruct the original inside and outside grid locations
  inside = zeros(Nsparse,1);
  for i=1:Nsparse
    fx = find(xgrid==sparsepos(i,1));
    fy = find(ygrid==sparsepos(i,2));
    fz = find(zgrid==sparsepos(i,3));
      inside(i) = (fx-1)*sx + (fy-1)*sy + (fz-1)*sz + 1;
  end
  outside = setdiff([1:Nfull]', inside);
end

fprintf('total number of dipoles        : %d\n', length(inside)+length(outside));
fprintf('number of dipoles inside  brain: %d\n', length(inside));
fprintf('number of dipoles outside brain: %d\n', length(outside));

% determine whether the source is old or new style
fnames = fieldnames(source);
if any(~cellfun('isempty', strfind(fnames, 'dimord'))),
  stype = 'new';
else
  stype = 'old';
end

if strcmp(stype, 'old')
  % original code
  % first do the non-trial fields
  source.dim = [1 length(inside) 1]; %to fool parameterselection
  [param]    = parameterselection('all', source);
  trlparam   = strmatch('trial', param);
  sel        = setdiff(1:length(param), trlparam);
  ind=find(ismember(param,'inside')); % find the index of 'inside' field
  % because its position varies with isfield('plvspctrm') vs. 'cohspctrm'
  param      = param(sel(ind));
  
  for j = 1:length(param)
    dat = getsubfield(source, param{j});
    if islogical(dat)
      tmp         = false(1,Nfull); 
      tmp(inside) = dat;
    elseif iscell(dat)
      tmp          = cell(1,Nfull);
      tmp(inside)  = dat;
      %tmp(outside) = nan;
    else
      tmp         = nan(1,Nfull);
      tmp(inside) = dat;   
    end
    source = setsubfield(source, param{j}, tmp);
  end
  
  % then do the trial fields
  if     isfield(source, 'trial' )
    for j = 1:length(source.trial)
      tmpsource     = source.trial(j);
      tmpsource.dim = source.dim; % to fool parameterselection
      tmpparam      = parameterselection('all', tmpsource);
      for k = 1:length(tmpparam)
        dat = getsubfield(tmpsource, tmpparam{k});
        if islogical(dat)
          tmp         = false(1,Nfull); 
          tmp(inside) = dat;
        elseif iscell(dat)
          tmp          = cell(1,Nfull);
          tmp(inside)  = dat;
          %tmp(outside) = nan;
        else
          tmp         = nan(1,Nfull);
          tmp(inside) = dat;   
        end
        tmpsource = setsubfield(tmpsource, tmpparam{k}, tmp);
      end
      tmpsource       = rmfield(tmpsource, 'dim');
      source.trial(j) = tmpsource;
    end   
  elseif isfield(source, 'trialA')
    for j = 1:length(source.trialA)
      tmpsource     = source.trialA(j);
      tmpsource.dim = source.dim; % to fool parameterselection
      tmpparam      = parameterselection('all', tmpsource);
      for k = 1:length(tmpparam)
        dat = getsubfield(tmpsource, tmpparam{k});
        if islogical(dat)
          tmp         = false(1,Nfull); 
          tmp(inside) = dat;
        elseif iscell(dat)
          tmp          = cell(1,Nfull);
          tmp(inside)  = dat;
          %tmp(outside) = nan;
        else
          tmp         = nan(1,Nfull);
          tmp(inside) = dat;   
        end
        tmpsource = setsubfield(tmpsource, tmpparam{k}, tmp);
      end
      tmpsource        = rmfield(tmpsource, 'dim');
      source.trialA(j) = tmpsource;   
    end
  elseif isfield(source, 'trialB')
    for j = 1:length(source.trialB)
      tmpsource     = source.trialB(j);
      tmpsource.dim = source.dim; % to fool parameterselection
      tmpparam      = parameterselection('all', tmpsource);
      for k = 1:length(tmpparam)
        dat = getsubfield(tmpsource, tmpparam{k});
        if islogical(dat)
          tmp         = false(1,Nfull); 
          tmp(inside) = dat;
        elseif iscell(dat)
          tmp          = cell(1,Nfull);
          tmp(inside)  = dat;
          %tmp(outside) = nan;
        else
          tmp         = nan(1,Nfull);
          tmp(inside) = dat;   
        end
        tmpsource = setsubfield(tmpsource, tmpparam{k}, tmp);
      end
      tmpsource        = rmfield(tmpsource, 'dim');
      source.trialB(j) = tmpsource;   
    end
  end
  
  % and finally do the coherence-like matrices (size Nvox X Nvox)
  fn = fieldnames(source);
  for i=1:length(fn)
    d = getfield(source, fn{i});
    m = size(d, 1);
    n = size(d, 2);
    if m==Nsparse && n==Nsparse
      tmp = nan(Nfull,Nfull);
      tmp(inside,inside) = d;
      source = setfield(source, fn{i}, tmp);
    end
  end
  
  % update the inside and outside definitions
  source.inside  = inside;
  source.outside = outside;
  source.pos     = pos;
  source.dim     = siz;
elseif strcmp(stype, 'new')
  % new style conversion
  fn = fieldnames(source);
  for i=1:numel(fn)
    if any(size(source.(fn{i}))==Nsparse)
      if iscell(source.(fn{i}))
        indx = find(size(source.(fn{i}))==Nsparse);
        if all(indx==1)
          tmp            = cell(Nfull,1);
          tmp(inside,1)  = source.(fn{i});
          source.(fn{i}) = tmp;
        elseif all(indx==2)
          tmp            = cell(1,Nfull);
          tmp(1,inside)  = source.(fn{i});
          source.(fn{i}) = tmp;
        else
          ft_warning('sparse to full conversion failed for field %s\n', fn{i});
        end
      else
        indx = find(size(source.(fn{i}))==Nsparse);
        if all(indx==1)
          tmpsiz = [size(source.(fn{i})) 1];
          tmp    = nan([Nfull tmpsiz(2:end)]);
          tmp(inside,:,:,:,:) = source.(fn{i});
        elseif all(indx==2)
          tmpsiz = [size(source.(fn{i})) 1];
          tmp    = nan([tmpsiz(1) Nfull tmpsiz(3:end)]);
          tmp(:,inside,:,:,:) = source.(fn{i});
        elseif all(indx==[1 2])
          % bivariate matrix
          tmpsiz = [size(source.(fn{i})) 1];
          tmp    = nan([Nfull Nfull tmpsiz(3:end)]);
          tmp(inside,inside,:,:,:) = source.(fn{i});
        else
          ft_warning('sparse to full conversion failed for field %s\n', fn{i});
        end
      end
      % nothing to do
    end
  end
  
  % update the inside and outside definitions and pos
  source.inside  = inside;
  source.outside = outside;
  source.pos     = pos;

end
cfg = [];
% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with MATLAB versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = source.cfg; end
% remember the exact configuration details in the output 
source.cfg = cfg;
                                                                                                                                                                                                                                                 req    = istrue(ft_getopt(cfg, 'avgoverfreq', false));
avgovertime    = istrue(ft_getopt(cfg, 'avgovertime', false));
avgoverrpt     = istrue(ft_getopt(cfg, 'avgoverrpt',  false));

% do a sanity check for the averaging options
if avgoverpos,     assert(haspos,     'there are no source positions, so averaging is not possible'); end
if avgoverchan,    assert(haschan,    'there is no channel dimension, so averaging is not possible'); end
if avgoverchancmb, assert(haschancmb, 'there are no channel combinations, so averaging is not possible'); end
if avgoverfreq,    assert(hasfreq,    'there is no frequency dimension, so averaging is not possible'); end
if avgovertime,    assert(hastime,    'there is no time dimension, so averaging over time is not possible'); end
if avgoverrpt,     assert(hasrpt||hasrpttap, 'there are no repetitions, so averaging is not possible'); end

% set averaging function
cfg.nanmean = ft_getopt(cfg, 'nanmean', 'no');
if strcmp(cfg.nanmean, 'yes')
  average = @nanmean;
else
  average = @mean;
end

% by default we keep most of the dimensions in the data structure when averaging over them
keepposdim     = istrue(ft_getopt(cfg, 'keepposdim',  true));
keepchandim    = istrue(ft_getopt(cfg, 'keepchandim', true));
keepchancmbdim = istrue(ft_getopt(cfg, 'keepchancmbdim', true));
keepfreqdim    = istrue(ft_getopt(cfg, 'keepfreqdim', true));
keeptimedim    = istrue(ft_getopt(cfg, 'keeptimedim', true));
keeprptdim     = istrue(ft_getopt(cfg, 'keeprptdim', ~avgoverrpt));

if ~keepposdim,     assert(avgoverpos,     'removing a dimension is only possible when averaging'); end
if ~keepchandim,    assert(avgoverchan,    'removing a dimension is only possible when averaging'); end
if ~keepchancmbdim, assert(avgoverchancmb, 'removing a dimension is only possible when averaging'); end
if ~keepfreqdim,    assert(avgoverfreq,    'removing a dimension is only possible when averaging'); end
if ~keeptimedim,    assert(avgovertime,    'removing a dimension is only possible when averaging'); end
if ~keeprptdim,     assert(avgoverrpt,     'removing a dimension is only possible when averaging'); end

% trim the selection to all inputs, rpt and rpttap are dealt with later
if hasspike,   [selspike,   cfg] = getselection_spike  (cfg, varargin{:}); end
if haspos,     [selpos,     cfg] = getselection_pos    (cfg, varargin{:}, cfg.tolerance, cfg.select); end
if haschan,    [selchan,    cfg] = getselection_chan   (cfg, varargin{:}, cfg.select); end
if haschancmb, [selchancmb, cfg] = getselection_chancmb(cfg, varargin{:}, cfg.select); end
if hasfreq,    [selfreq,    cfg] = getselection_freq   (cfg, varargin{:}, cfg.tolerance, cfg.select); end
if hastime,    [seltime,    cfg] = getselection_time   (cfg, varargin{:}, cfg.tolerance, cfg.select); end

% this is to keep track of all fields that should be retained in the output
keepfield = datfield;

for i=1:numel(varargin)
  
  for j=1:numel(datfield)
    dimtok = tokenize(dimord{j}, '_');
    
    % the rpt selection should only work with a single data argument
    % in case tapers were kept, selrpt~=selrpttap, otherwise selrpt==selrpttap
    [selrpt{i}, dum, rptdim{i}, selrpttap{i}] = getselection_rpt(cfg, varargin{i}, dimord{j});
    
    % check for the presence of each dimension in each datafield
    fieldhasspike   = ismember('spike',   dimtok);
    fieldhaspos     = ismember('pos',     dimtok) || ismember('{pos}', dimtok);
    fieldhaschan    = (ismember('chan',    dimtok) || ismember('{chan}', dimtok)) && isfield(varargin{1}, 'label');
    fieldhaschancmb = ismember('chancmb', dimtok);
    fieldhastime    = ismember('time',    dimtok) && ~hasspike;
    fieldhasfreq    = ismember('freq',    dimtok);
    fieldhasrpt     = ismember('rpt',     dimtok) | ismember('subj', dimtok) | ismember('{rpt}', dimtok);
    fieldhasrpttap  = ismember('rpttap',  dimtok);
    
    % cfg.latency is used to select individual spikes, not to select from a continuously sampled time axis
    
    if fieldhasspike,   varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(strcmp(dimtok,'spike')),             selspike{i},   false,           'intersect', average); end
    if fieldhaspos,     varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(ismember(dimtok, {'pos', '{pos}'})), selpos{i},     avgoverpos,      cfg.select, average);  end
    if fieldhaschan,    varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(ismember(dimtok,{'chan' '{chan}'})), selchan{i},    avgoverchan,     cfg.select, average);  end
    if fieldhaschancmb, varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(strcmp(dimtok,'chancmb')),           selchancmb{i}, avgoverchancmb,  cfg.select, average);  end
    if fieldhastime,    varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(strcmp(dimtok,'time')),              seltime{i},    avgovertime,     cfg.select, average);  end
    if fieldhasfreq,    varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(strcmp(dimtok,'freq')),              selfreq{i},    avgoverfreq,     cfg.select, average);  end
    if fieldhasrpt,     varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, rptdim{i},                                selrpt{i},     avgoverrpt,      'intersect', average); end
    if fieldhasrpttap,  varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, rptdim{i},                                selrpttap{i},  avgoverrpt,      'intersect', average); end
    
    % update the fields that should be kept in the structure as a whole
    % and update the dimord for this specific datfield
    keepdim = true(size(dimtok));
    
    if avgoverchan && ~keepchandim
      keepdim(strcmp(dimtok, 'chan')) = false;
      keepfield = setdiff(keepfield, 'label');
    else
      keepfield = [keepfield 'label'];
    end
    
    if avgoverchancmb && ~keepchancmbdim
      keepdim(strcmp(dimtok, 'chancmb')) = false;
      keepfield = setdiff(keepfield, 'labelcmb');
    else
      keepfield = [keepfield 'labelcmb'];
    end
    
    if avgoverfreq && ~keepfreqdim
      keepdim(strcmp(dimtok, 'freq')) = false;
      keepfield = setdiff(keepfield, 'freq');
    else
      keepfield = [keepfield 'freq'];
    end
    
    if avgovertime && ~keeptimedim
      keepdim(strcmp(dimtok, 'time')) = false;
      keepfield = setdiff(keepfield, 'time');
    else
      keepfield = [keepfield 'time'];
    end
    
    if avgoverpos && ~keepposdim
      keepdim(strcmp(dimtok, 'pos'))   = false;
      keepdim(strcmp(dimtok, '{pos}')) = false;
      keepdim(strcmp(dimtok, 'dim'))   = false;
      keepfield = setdiff(keepfield, {'pos' '{pos}' 'dim'});
    elseif avgoverpos && keepposdim
      keepfield = setdiff(keepfield, {'dim'}); % this should be removed anyway
    else
      keepfield = [keepfield {'pos' '{pos}' 'dim'}];
    end
    
    if avgoverrpt && ~keeprptdim
      keepdim(strcmp(dimtok, 'rpt'))    = false;
      keepdim(strcmp(dimtok, 'rpttap')) = false;
      keepdim(strcmp(dimtok, 'subj'))   = false;
    end
    
    % update the sampleinfo, if possible, and needed
    if strcmp(datfield{j}, 'sampleinfo') && ~isequal(cfg.latency, 'all')
      if iscell(seltime{i}) && numel(seltime{i})==size(varargin{i}.sampleinfo,1)
        for k = 1:numel(seltime{i})
          varargin{i}.sampleinfo(k,:) = varargin{i}.sampleinfo(k,1) - 1 + seltime{i}{k}([1 end]);
        end
      elseif ~iscell(seltime{i}) && ~isempty(seltime{i}) && ~all(isnan(seltime{i}))
        nrpt       = size(varargin{i}.sampleinfo,1);
        seltime{i} = seltime{i}(:)';
        varargin{i}.sampleinfo = varargin{i}.sampleinfo(:,[1 1]) - 1 + repmat(seltime{i}([1 end]),nrpt,1);
      end
    end
    
    varargin{i}.(datfield{j}) = squeezedim(varargin{i}.(datfield{j}), ~keepdim);
    
  end % for datfield
  
  % also update the fields that describe the dimensions, time/freq/pos have been dealt with as data
  if haschan,    varargin{i} = makeselection_chan   (varargin{i}, selchan{i}, avgoverchan); end % update the label field
  if haschancmb, varargin{i} = makeselection_chancmb(varargin{i}, selchancmb{i}, avgoverchancmb); end % update the labelcmb field
  
end % for varargin

if strcmp(cfg.select, 'union')
  % create the union of the descriptive axes
  if haspos,      varargin = makeunion(varargin, 'pos'); end
  if haschan,     varargin = makeunion(varargin, 'label'); end
  if haschancmb,  varargin = makeunion(varargin, 'labelcmb'); end
  if hastime,     varargin = makeunion(varargin, 'time'); end
  if hasfreq,     varargin = makeunion(varargin, 'freq'); end
end

% remove all fields from the data structure that do not pertain to the selection
sel = strcmp(keepfield, '{pos}'); if any(sel), keepfield(sel) = {'pos'}; end
sel = strcmp(keepfield, 'chan');  if any(sel), keepfield(sel) = {'label'}; end
sel = strcmp(keepfield, 'chancmb');  if any(sel), keepfield(sel) = {'labelcmb'}; end

if avgoverrpt
  % these are invalid after averaging
  keepfield = setdiff(keepfield, {'cumsumcnt' 'cumtapcnt' 'trialinfo' 'sampleinfo'});
end

if avgovertime
  % these are invalid after averaging or making a latency selection
  keepfield = setdiff(keepfield, {'sampleinfo'});
end

for i=1:numel(varargin)
  varargin{i} = keepfields(varargin{i}, [keepfield ignorefields('selectdata')']);
end

% restore the original dimord fields in the data
for i=1:length(orgdim1)
  dimtok = tokenize(orgdim2{i}, '_');
  
  % using a setdiff may result in double occurrences of chan and pos to
  % disappear, so this causes problems as per bug 2962
  % if ~keeprptdim, dimtok = setdiff(dimtok, {'rpt' 'rpttap' 'subj'}); end
  % if ~keepposdim, dimtok = setdiff(dimtok, {'pos' '{pos}'}); end
  % if ~keepchandim, dimtok = setdiff(dimtok, {'chan'}); end
  % if ~keepfreqdim, dimtok = setdiff(dimtok, {'freq'}); end
  % if ~keeptimedim, dimtok = setdiff(dimtok, {'time'}); end
  if ~keeprptdim, dimtok = dimtok(~ismember(dimtok, {'rpt' 'rpttap' 'subj'})); end
  if ~keepposdim, dimtok = dimtok(~ismember(dimtok, {'pos' '{pos}'}));         end
  if ~keepchandim, dimtok = dimtok(~ismember(dimtok, {'chan'})); end
  if ~keepfreqdim, dimtok = dimtok(~ismember(dimtok, {'freq'})); end
  if ~keeptimedim, dimtok = dimtok(~ismember(dimtok, {'time'})); end
  
  dimord = sprintf('%s_', dimtok{:});
  dimord = dimord(1:end-1); % remove the trailing _
  for j=1:length(varargin)
    varargin{j}.(orgdim1{i}) = dimord;
  end
end

% restore the source.avg field, this keeps the output reasonably consistent with the
% old-style source representation of the input
if strcmp(dtype, 'source') && ~isempty(restoreavg)
  for i=1:length(varargin)
    varargin{i}.avg = keepfields(varargin{i}, restoreavg);
    varargin{i}     = removefields(varargin{i}, restoreavg);
  end
end

varargout = varargin;

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance varargout
ft_postamble history varargout
ft_postamble savevar varargout

% the varargout variable can be cleared when written to outputfile
if exist('varargout', 'var') && ft_nargout>numel(varargout)
  % also return the input cfg with the combined selection over all input data structures
  varargout{end+1} = cfg;
end

end % function ft_selectdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection(data, datfield, dimtok, seldim, selindx, avgoverdim, selmode, average)

if numel(seldim) > 1
  for k = 1:numel(seldim)
    data = makeselection(data, datfield, dimtok, seldim(k), selindx, avgoverdim, selmode, average);
  end
  return;
end

if isnumeric(data.(datfield))
  if isrow(data.(datfield)) && seldim==1
    if length(dimtok)==1
      seldim = 2; % switch row and column
    end
  elseif iscolumn(data.(datfield)) && seldim==2
    if length(dimtok)==1
      seldim = 1; % switch row and column
    end
  end
elseif iscell(data.(datfield))
  if isrow(data.(datfield){1}) && seldim==2
    if length(dimtok)==2
      seldim = 3; % switch row and column
    end
  elseif iscolumn(data.(datfield){1}) && seldim==3
    if length(dimtok)==2
      seldim = 2; % switch row and column
    end
  end
end

% an empty selindx means that nothing(!) should be selected and hence everything should be removed, which is different than keeping everything
% the selindx value of NaN indicates that it is not needed to make a selection

switch selmode
  case 'intersect'
    if iscell(selindx)
      % there are multiple selections in multipe vectors, the selection is in the matrices contained within the cell-array
      for j=1:numel(selindx)
        if ~isempty(selindx{j}) && all(isnan(selindx{j}))
          % no selection needs to be made
        else
          data.(datfield){j} = cellmatselect(data.(datfield){j}, seldim-1, selindx{j}, numel(dimtok)==1);
        end
      end
      
    else
      % there is a single selection in a single vector
      if ~isempty(selindx) && all(isnan(selindx))
        % no selection needs to be made
      else
        data.(datfield) = cellmatselect(data.(datfield), seldim, selindx, numel(dimtok)==1);
      end
    end
    
    if avgoverdim
      data.(datfield) = cellmatmean(data.(datfield), seldim, average);
    end
    
  case 'union'
    if ~isempty(selindx) && all(isnan(selindx))
      % no selection needs to be made
    else
      tmp = data.(datfield);
      siz = size(tmp);
      siz(seldim) = numel(selindx);
      data.(datfield) = nan(siz);
      sel = isfinite(selindx);
      switch seldim
        case 1
          data.(datfield)(sel,:,:,:,:,:) = tmp(selindx(sel),:,:,:,:,:);
        case 2
          data.(datfield)(:,sel,:,:,:,:) = tmp(:,selindx(sel),:,:,:,:);
        case 3
          data.(datfield)(:,:,sel,:,:,:) = tmp(:,:,selindx(sel),:,:,:);
        case 4
          data.(datfield)(:,:,:,sel,:,:) = tmp(:,:,:,selindx(sel),:,:);
        case 5
          data.(datfield)(:,:,:,:,sel,:) = tmp(:,:,:,:,selindx(sel),:);
        case 6
          data.(datfield)(:,:,:,:,:,sel) = tmp(:,:,:,:,:,selindx(sel));
        otherwise
          ft_error('unsupported dimension (%d) for making a selection for %s', seldim, datfield);
      end
    end
    
    if avgoverdim
      data.(datfield) = average(data.(datfield), seldim);
    end
end % switch

end % function makeselection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_chan(data, selchan, avgoverchan)
if isempty(selchan)
  %error('no channels were selected');
  data.label = {};
elseif avgoverchan && all(isnan(selchan))
  str = sprintf('%s, ', data.label{:});
  str = str(1:end-2);
  str = sprintf('mean(%s)', str);
  data.label = {str};
elseif avgoverchan && ~any(isnan(selchan))
  str = sprintf('%s, ', data.label{selchan});
  str = str(1:end-2);
  str = sprintf('mean(%s)', str);
  data.label = {str};                 % remove the last '+'
elseif all(isfinite(selchan))
  data.label = data.label(selchan);
  data.label = data.label(:);
elseif numel(selchan)==1 && any(~isfinite(selchan))
  % do nothing
elseif numel(selchan)>1  && any(~isfinite(selchan))
  tmp = cell(numel(selchan),1);
  for k = 1:numel(tmp)
    if isfinite(selchan(k))
      tmp{k} = data.label{selchan(k)};
    end
  end
  data.label = tmp;
else
  % this should never happen
  ft_error('cannot figure out how to select channels');
end
end % function makeselection_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_chancmb(data, selchancmb, avgoverchancmb)
if isempty(selchancmb)
  ft_error('no channel combinations were selected');
elseif avgoverchancmb && all(isnan(selchancmb))
  % naming the channel combinations becomes ambiguous, but should not
  % suggest that the mean was computed prior to combining
  str1 = sprintf('%s, ', data.labelcmb{:,1});
  str1 = str1(1:end-2);
  % str1 = sprintf('mean(%s)', str1);
  str2 = sprintf('%s, ', data.labelcmb{:,2});
  str2 = str2(1:end-2);
  % str2 = sprintf('mean(%s)', str2);
  data.label = {str1, str2};
elseif avgoverchancmb && ~any(isnan(selchancmb))
  % naming the channel combinations becomes ambiguous, but should not
  % suggest that the mean was computed prior to combining
  str1 = sprintf('%s, ', data.labelcmb{selchancmb,1});
  str1 = str1(1:end-2);
  % str1 = sprintf('mean(%s)', str1);
  str2 = sprintf('%s, ', data.labelcmb{selchancmb,2});
  str2 = str2(1:end-2);
  % str2 = sprintf('mean(%s)', str2);
  data.label = {str1, str2};
elseif all(isfinite(selchancmb))
  data.labelcmb = data.labelcmb(selchancmb,:);
elseif numel(selchancmb)==1 && any(~isfinite(selchancmb))
  % do nothing
elseif numel(selchancmb)>1  && any(~isfinite(selchancmb))
  tmp = cell(numel(selchancmb),2);
  for k = 1:size(tmp,1)
    if isfinite(selchan(k))
      tmp{k,1} = data.labelcmb{selchan(k),1};
      tmp{k,2} = data.labelcmb{selchan(k),2};
    end
  end
  data.labelcmb = tmp;
else
  % this should never happen
  ft_error('cannot figure out how to select channelcombinations');
end
end % function makeselection_chancmb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chanindx, cfg] = getselection_chan(cfg, varargin)

selmode  = varargin{end};
ndata    = numel(varargin)-1;
varargin = varargin(1:ndata);

% loop over data once to initialize
chanindx = cell(ndata,1);
label    = cell(1,0);

for k = 1:ndata
  if isfield(varargin{k}, 'grad') && isfield(varargin{k}.grad, 'type')
    % this makes channel selection more robust
    selchannel = ft_channelselection(cfg.channel, varargin{k}.label, varargin{k}.grad.type);
  elseif isfield(varargin{k}, 'elec') && isfield(varargin{k}.elec, 'type')
    % this makes channel selection more robust
    selchannel = ft_channelselection(cfg.channel, varargin{k}.label, varargin{k}.elec.type);
  else
    selchannel = ft_channelselection(cfg.channel, varargin{k}.label);
  end
  label      = union(label, selchannel);
end
label = label(:);   % ensure column array

% this call to match_str ensures that that labels are always in the
% order of the first input argument see bug_2917, but also temporarily keep
% the labels from the other data structures not present in the first one
% (in case selmode = 'union')
[ix, iy] = match_str(varargin{1}.label, label);
label1   = varargin{1}.label(:); % ensure column array
label    = [label1(ix); label(setdiff(1:numel(label),iy))];

indx = nan+zeros(numel(label), ndata);
for k = 1:ndata
  [ix, iy] = match_str(label, varargin{k}.label);
  indx(ix,k) = iy;
end

switch selmode
  case 'intersect'
    sel      = sum(isfinite(indx),2)==ndata;
    indx     = indx(sel,:);
    label    = varargin{1}.label(indx(:,1));
  case 'union'
    % don't do a subselection
  otherwise
    ft_error('invalid value for cfg.select');
end % switch

ok = false(size(indx,1),1);
for k = 1:ndata
  % loop through the columns to preserve the order of the channels, where
  % the order of the input arguments determines the final order
  ix = find(~ok);
  [srt,srtix] = sort(indx(ix,k));
  indx(ix,:)  = indx(ix(srtix),:);
  ok = ok | isfinite(indx(:,k));
end

for k = 1:ndata
  % do a sanity check on double occurrences
  if numel(unique(indx(isfinite(indx(:,k)),k)))<sum(isfinite(indx(:,k)))
    ft_error('the selection of channels across input arguments leads to double occurrences');
  end
  chanindx{k} = indx(:,k);
end

for k = 1:ndata
  if isequal(chanindx{k}, (1:numel(varargin{k}.label))')
    % no actual selection is needed for this data structure
    chanindx{k} = nan;
  end
end

cfg.channel = label;

end % function getselection_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chancmbindx, cfg] = getselection_chancmb(cfg, varargin)

selmode  = varargin{end};
ndata    = numel(varargin)-1;
varargin = varargin(1:ndata);

chancmbindx = cell(ndata,1);

if ~isfield(cfg, 'channelcmb')
  for k=1:ndata
    % the nan return value specifies that no selection was specified
    chancmbindx{k} = nan;
  end
  
else
  
  switch selmode
    case 'intersect'
      for k=1:ndata
        if ~isfield(varargin{k}, 'label')
          cfg.channelcmb = ft_channelcombination(cfg.channelcmb, unique(varargin{k}.labelcmb(:)));
        else
          cfg.channelcmb = ft_channelcombination(cfg.channelcmb, varargin{k}.label);
        end
      end
      
      ncfgcmb = size(cfg.channelcmb,1);
      cfgcmb  = cell(ncfgcmb, 1);
      for i=1:ncfgcmb
        cfgcmb{i} = sprintf('%s&%s', cfg.channelcmb{i,1}, cfg.channelcmb{i,2});
      end
      
      for k=1:ndata
        ndatcmb = size(varargin{k}.labelcmb,1);
        datcmb = cell(ndatcmb, 1);
        for i=1:ndatcmb
          datcmb{i} = sprintf('%s&%s', varargin{k}.labelcmb{i,1}, varargin{k}.labelcmb{i,2});
        end
        
        % return the order according to the (joint) configuration, not according to the (individual) data
        % FIXME this should adhere to the general code guidelines, where
        % the order returned will be according to the first data argument!
        [dum, chancmbindx{k}] = match_str(cfgcmb, datcmb);
      end
      
    case 'union'
      % FIXME this is not yet implemented
      ft_error('union of channel combination is not yet supported');
      
    otherwise
      ft_error('invalid value for cfg.select');
  end % switch
  
  
end

end % function getselection_chancmb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikeindx, cfg] = getselection_spike(cfg, varargin)
% possible specifications are
% cfg.latency = string -> 'all'
% cfg.latency = [beg end]
% cfg.trials  = string -> 'all'
% cfg.trials  = vector with indices

ndata    = numel(varargin);
varargin = varargin(1:ndata);

if isequal(cfg.latency, 'all') && isequal(cfg.trials, 'all')
  spikeindx = cell(1,ndata);
  for i=1:ndata
    spikeindx{i} = num2cell(nan(1, length(varargin{i}.time)));
  end
  return
end

trialbeg = varargin{1}.trialtime(:,1);
trialend = varargin{1}.trialtime(:,2);
for i=2:ndata
  trialbeg = cat(1, trialbeg, varargin{1}.trialtime(:,1));
  trialend = cat(1, trialend, varargin{1}.trialtime(:,2));
end

% convert string into a numeric selection
if ischar(cfg.latency)
  switch cfg.latency
    case 'all'
      cfg.latency = [-inf inf];
    case 'maxperiod'
      cfg.latency = [min(trialbeg) max(trialend)];
    case 'minperiod'
      cfg.latency = [max(trialbeg) min(trialend)];
    case 'prestim'
      cfg.latency = [min(trialbeg) 0];
    case 'poststim'
      cfg.latency = [0 max(trialend)];
    otherwise
      ft_error('incorrect specification of cfg.latency');
  end % switch
end

spikeindx = cell(1,ndata);
for i=1:ndata
  nchan = length(varargin{i}.time);
  spikeindx{i} = cell(1,nchan);
  for j=1:nchan
    selbegtime = varargin{i}.time{j}>=cfg.latency(1);
    selendtime = varargin{i}.time{j}<=cfg.latency(2);
    if isequal(cfg.trials, 'all')
      seltrial = true(size(varargin{i}.trial{j}));
    else
      seltrial = ismember(varargin{i}.trial{j}, cfg.trials);
    end
    spikeindx{i}{j} = find(selbegtime & selendtime & seltrial);
  end
end

end % function getselection_spiketime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeindx, cfg] = getselection_time(cfg, varargin)
% possible specifications are
% cfg.latency = value     -> can be 'all'
% cfg.latency = [beg end]

if ft_datatype(varargin{1}, 'spike')
  ft_error('latency selection in spike data is not supported')
end

selmode  = varargin{end};
tol      = varargin{end-1};
ndata    = numel(varargin)-2;
varargin = varargin(1:ndata);

if isequal(cfg.latency, 'all') && iscell(varargin{1}.time)
  % for raw data this means that all trials should be selected as they are
  % for timelock/freq data it is still needed to make the intersection between data arguments
  timeindx = cell(1,ndata);
  for i=1:ndata
    % the nan return value specifies that no selection was specified
    timeindx{i} = num2cell(nan(1, length(varargin{i}.time)));
  end
  return
end

% if there is a single timelock/freq input, there is one time vector
% if there are multiple timelock/freq inputs, there are multiple time vectors
% if there is a single raw input, there are multiple time vectors
% if there are multiple raw inputs, there are multiple time vectors

% collect all time axes in one large cell-array
alltimecell = {};
if iscell(varargin{1}.time)
  for k = 1:ndata
    alltimecell = [alltimecell varargin{k}.time{:}];
  end
else
  for k = 1:ndata
    alltimecell = [alltimecell {varargin{k}.time}];
  end
end

% the nan return value specifies that no selection was specified
timeindx = repmat({nan}, size(alltimecell));

% loop over data once to determine the union of all time axes
alltimevec = zeros(1,0);
for k = 1:length(alltimecell)
  alltimevec = union(alltimevec, round(alltimecell{k}/tol)*tol);
end

indx = nan(numel(alltimevec), numel(alltimecell));
for k = 1:numel(alltimecell)
  [dum, ix, iy] = intersect(alltimevec, round(alltimecell{k}/tol)*tol);
  indx(ix,k) = iy;
end

if iscell(varargin{1}.time) && ischar(cfg.latency)&& ~strcmp(cfg.latency, 'minperiod')
  % if the input data arguments are of type 'raw', temporarily set the
  % selmode to union, otherwise the potentially different length trials
  % will be truncated to the shorted epoch, prior to latency selection.
  selmode = 'union';
elseif ischar(cfg.latency) && strcmp(cfg.latency, 'minperiod')
  % enforce intersect
  selmode = 'intersect';
end
switch selmode
  case 'intersect'
    sel        = sum(isfinite(indx),2)==numel(alltimecell);
    indx       = indx(sel,:);
    alltimevec = alltimevec(sel);
  case 'union'
    % don't do a subselection
  otherwise
    ft_error('invalid value for cfg.select');
end

% Note that cfg.toilim handling has been removed, as it was renamed to cfg.latency

% convert a string selection into a numeric selection
if ischar(cfg.latency)
  switch cfg.latency
    case {'all' 'maxperlen' 'maxperiod'}
      cfg.latency = [min(alltimevec) max(alltimevec)];
    case 'prestim'
      cfg.latency = [min(alltimevec) 0];
    case 'poststim'
      cfg.latency = [0 max(alltimevec)];
    case 'minperiod'
      % the time vector has been pruned above
      cfg.latency = [min(alltimevec) max(alltimevec)];
    otherwise
      ft_error('incorrect specification of cfg.latency');
  end % switch
end

% deal with numeric selection
if isempty(cfg.latency)
  for k = 1:numel(alltimecell)
    % FIXME I do not understand this
    % this signifies that all time bins are deselected and should be removed
    timeindx{k} = [];
  end
  
elseif numel(cfg.latency)==1
  % this single value should be within the time axis of each input data structure
  if numel(alltimevec)>1
    tbin = nearest(alltimevec, cfg.latency, true, true); % determine the numerical tolerance
  else
    tbin = nearest(alltimevec, cfg.latency, true, false); % don't consider tolerance
  end
  cfg.latency = alltimevec(tbin);
  
  for k = 1:ndata
    timeindx{k} = indx(tbin, k);
  end
  
elseif numel(cfg.latency)==2
  % the [min max] range can be specifed with +inf or -inf, but should
  % at least partially overlap with the time axis of the input data
  mintime = min(alltimevec);
  maxtime = max(alltimevec);
  if all(cfg.latency<mintime) || all(cfg.latency>maxtime)
    ft_error('the selected time range falls outside the time axis in the data');
  end
  tbeg = nearest(alltimevec, cfg.latency(1), false, false);
  tend = nearest(alltimevec, cfg.latency(2), false, false);
  cfg.latency = alltimevec([tbeg tend]);
  
  for k = 1:numel(alltimecell)
    timeindx{k} = indx(tbeg:tend, k);
    % if the input data arguments are of type 'raw', the non-finite values
    % need to be removed from the individual cells to ensure correct
    % behavior
    if iscell(varargin{1}.time)
      timeindx{k} = timeindx{k}(isfinite(timeindx{k}));
    end
  end
  
elseif size(cfg.latency,2)==2
  % this may be used for specification of the computation, not for data selection
  
else
  ft_error('incorrect specification of cfg.latency');
end

for k = 1:numel(alltimecell)
  if ~iscell(varargin{1}.time)
    if isequal(timeindx{k}(:)', 1:length(alltimecell{k}))
      % no actual selection is needed for this data structure
      timeindx{k} = nan;
    end
  else
    % if the input data arguments are of type 'raw', they need to be
    % handled differently, because the individual trials can be of
    % different length
  end
end

if iscell(varargin{1}.time)
  % split all time axes again over the different input raw data structures
  dum = cell(1,ndata);
  for k = 1:ndata
    sel = 1:length(varargin{k}.time);
    dum{k} = timeindx(sel); % get the first selection
    timeindx(sel) = []; % remove the first selection
  end
  timeindx = dum;
else
  % no splitting is needed, each input data structure has one selection
end

end % function getselection_time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freqindx, cfg] = getselection_freq(cfg, varargin)
% possible specifications are
% cfg.frequency = value     -> can be 'all'
% cfg.frequency = [beg end]

selmode  = varargin{end};
tol      = varargin{end-1};
ndata    = numel(varargin)-2;
varargin = varargin(1:ndata);

% loop over data once to initialize
freqindx = cell(ndata,1);
freqaxis = zeros(1,0);
for k = 1:ndata
  % the nan return value specifies that no selection was specified
  freqindx{k} = nan;
  
  % update the axis along which the frequencies are defined
  freqaxis = union(freqaxis, round(varargin{k}.freq(:)/tol)*tol);
end

indx = nan+zeros(numel(freqaxis), ndata);
for k = 1:ndata
  [dum, ix, iy] = intersect(freqaxis, round(varargin{k}.freq(:)/tol)*tol);
  indx(ix,k) = iy;
end

switch selmode
  case 'intersect'
    sel      = sum(isfinite(indx),2)==ndata;
    indx     = indx(sel,:);
    freqaxis = varargin{1}.freq(indx(:,1));
  case 'union'
    % don't do a subselection
  otherwise
    ft_error('invalid value for cfg.select');
end

if isfield(cfg, 'frequency')
  % deal with string selection
  % some of these do not make sense, but are here for consistency with ft_multiplotER
  if ischar(cfg.frequency)
    if strcmp(cfg.frequency, 'all')
      cfg.frequency = [min(freqaxis) max(freqaxis)];
    elseif strcmp(cfg.frequency, 'maxmin')
      cfg.frequency = [min(freqaxis) max(freqaxis)]; % the same as 'all'
    elseif strcmp(cfg.frequency, 'minzero')
      cfg.frequency = [min(freqaxis) 0];
    elseif strcmp(cfg.frequency, 'maxabs')
      cfg.frequency = [-max(abs(freqaxis)) max(abs(freqaxis))];
    elseif strcmp(cfg.frequency, 'zeromax')
      cfg.frequency = [0 max(freqaxis)];
    elseif strcmp(cfg.frequency, 'zeromax')
      cfg.frequency = [0 max(freqaxis)];
    else
      ft_error('incorrect specification of cfg.frequency');
    end
  end
  
  % deal with numeric selection
  if isempty(cfg.frequency)
    for k = 1:ndata
      % FIXME I do not understand this
      % this signifies that all frequency bins are deselected and should be removed
      freqindx{k} = [];
    end
    
  elseif numel(cfg.frequency)==1
    % this single value should be within the frequency axis of each input data structure
    if numel(freqaxis)>1
      fbin = nearest(freqaxis, cfg.frequency, true, true); % determine the numerical tolerance
    else
      fbin = nearest(freqaxis, cfg.frequency, true, false); % don't consider tolerance
    end
    cfg.frequency = freqaxis(fbin);
    
    for k = 1:ndata
      freqindx{k} = indx(fbin,k);
    end
    
  elseif numel(cfg.frequency)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the freq axis of the input data
    minfreq = min(freqaxis);
    maxfreq = max(freqaxis);
    if all(cfg.frequency<minfreq) || all(cfg.frequency>maxfreq)
      ft_error('the selected range falls outside the frequency axis in the data');
    end
    fbeg = nearest(freqaxis, cfg.frequency(1), false, false);
    fend = nearest(freqaxis, cfg.frequency(2), false, false);
    cfg.frequency = freqaxis([fbeg fend]);
    
    for k = 1:ndata
      freqindx{k} = indx(fbeg:fend,k);
    end
    
  elseif size(cfg.frequency,2)==2
    % this may be used for specification of the computation, not for data selection
    
  else
    ft_error('incorrect specification of cfg.frequency');
  end
end % if cfg.frequency

for k = 1:ndata
  if isequal(freqindx{k}, 1:length(varargin{k}.freq))
    % the cfg was updated, but no selection is needed for the data
    freqindx{k} = nan;
  end
end

end % function getselection_freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rptindx, cfg, rptdim, rpttapindx] = getselection_rpt(cfg, varargin)
% this should deal with cfg.trials

dimord   = varargin{end};
ndata    = numel(varargin)-1;
data     = varargin{1:ndata}; % this syntax ensures that it will only work on a single data input

dimtok = tokenize(dimord, '_');
rptdim = find(strcmp(dimtok, '{rpt}') | strcmp(dimtok, 'rpt') | strcmp(dimtok, 'rpttap') | strcmp(dimtok, 'subj'));

if isequal(cfg.trials, 'all')
  rptindx    = nan; % the nan return value specifies that no selection was specified
  rpttapindx = nan; % the nan return value specifies that no selection was specified
  
elseif isempty(rptdim)
  % FIXME should [] not mean that none of the trials is to be selected?
  rptindx    = nan; % the nan return value specifies that no selection was specified
  rpttapindx = nan; % the nan return value specifies that no selection was specified
  
else
  rptindx = ft_getopt(cfg, 'trials');
  
  if islogical(rptindx)
    % convert from booleans to indices
    rptindx = find(rptindx);
  end
  
  rptindx = unique(sort(rptindx));
  
  if strcmp(dimtok{rptdim}, 'rpttap') && isfield(data, 'cumtapcnt')
    % there are tapers in the data
    
    % determine for each taper to which trial it belongs
    nrpt = size(data.cumtapcnt, 1);
    taper = zeros(nrpt, 1);
    sumtapcnt = cumsum([0; data.cumtapcnt(:)]);
    begtapcnt = sumtapcnt(1:end-1)+1;
    endtapcnt = sumtapcnt(2:end);
    for i=1:nrpt
      taper(begtapcnt(i):endtapcnt(i)) = i;
    end
    rpttapindx = find(ismember(taper, rptindx));
    
  else
    % there are no tapers in the data
    rpttapindx = rptindx;
  end
end

end % function getselection_rpt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posindx, cfg] = getselection_pos(cfg, varargin)
% possible specifications are <none>

ndata   = numel(varargin)-2;
tol     = varargin{end-1}; % FIXME this is still ignored
selmode = varargin{end};   % FIXME this is still ignored
data    = varargin(1:ndata);

for i=1:ndata
  if ~isequal(varargin{i}.pos, varargin{1}.pos)
    % FIXME it would be possible here to make a selection based on intersect or union
    ft_error('not yet implemented');
  end
end

if strcmp(cfg.select, 'union')
  % FIXME it would be possible here to make a selection based on intersect or union
  ft_error('not yet implemented');
end

for i=1:ndata
  posindx{i} = nan;    % the nan return value specifies that no selection was specified
end
end % function getselection_pos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = squeezedim(x, dim)
siz = size(x);
for i=(numel(siz)+1):numel(dim)
  % all trailing singleton dimensions have length 1
  siz(i) = 1;
end
if isvector(x) && ~(isrow(x) && dim(1) && numel(x)>1)
  % there is no harm to keep it as it is, unless the data matrix is 1xNx1x1
elseif istable(x)
  % there is no harm to keep it as it is
else
  x = reshape(x, [siz(~dim) 1]);
end
end % function squeezedim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = makeunion(x, field)
old = cellfun(@getfield, x, repmat({field}, size(x)), 'uniformoutput', false);
if iscell(old{1})
  % empty is indicated to represent missing value for a cell-array (label, labelcmb)
  new = old{1};
  for i=2:length(old)
    sel = ~cellfun(@isempty, old{i});
    new(sel) = old{i}(sel);
  end
else
  % nan is indicated to represent missing value for a numeric array (time, freq, pos)
  new = old{1};
  for i=2:length(old)
    sel = ~isnan(old{i});
    new(sel) = old{i}(sel);
  end
end
x = cellfun(@setfield, x, repmat({field}, size(x)), repmat({new}, size(x)), 'uniformoutput', false);
end % function makeunion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to make a selextion in data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = cellmatselect(x, seldim, selindx, maybevector)
if nargin<4
  % some fields are a vector with an unspecified singleton dimension, these can be transposed
  % if the singleton dimension represents something explicit, they should not be transposed
  % they might for example represent a single trial, or a single channel
  maybevector = true;
end
if iscell(x)
  if seldim==1
    x = x(selindx);
  else
    for i=1:numel(x)
      if isempty(x{i})
        continue
      end
      switch seldim
        case 2
          if maybevector && isvector(x{i})
            % sometimes the data is 1xN, whereas the dimord describes only the first dimension
            % in this case a row and column vector can be interpreted as equivalent
            x{i} = x{i}(selindx);
          elseif istable(x)
            % multidimensional indexing is not supported
            x{i} = x{i}(selindx,:);
          else
            x{i} = x{i}(selindx,:,:,:,:);
          end
        case 3
          x{i} = x{i}(:,selindx,:,:,:);
        case 4
          x{i} = x{i}(:,:,selindx,:,:);
        case 5
          x{i} = x{i}(:,:,:,selindx,:);
        case 6
          x{i} = x{i}(:,:,:,:,selindx);
        otherwise
          ft_error('unsupported dimension (%d) for making a selection', seldim);
      end % switch
    end % for
  end
else
  switch seldim
    case 1
      if maybevector && isvector(x)
        % sometimes the data is 1xN, whereas the dimord describes only the first dimension
        % in this case a row and column vector can be interpreted as equivalent
        x = x(selindx);
      elseif istable(x)
        % multidimensional indexing is not supported
        x = x(selindx,:);
      else
        x = x(selindx,:,:,:,:,:);
      end
    case 2
      x = x(:,selindx,:,:,:,:);
    case 3
      x = x(:,:,selindx,:,:,:);
    case 4
      x = x(:,:,:,selindx,:,:);
    case 5
      x = x(:,:,:,:,selindx,:);
    case 6
      x = x(:,:,:,:,:,selindx);
    otherwise
      ft_error('unsupported dimension (%d) for making a selection', seldim);
  end
end
end % function cellmatselect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to take an average in data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = cellmatmean(x, seldim, average)
if iscell(x)
  if seldim==1
    for i=2:numel(x)
      x{1} = x{1} + x{i};
    end
    x = {x{1}/numel(x)};
  else
    for i=1:numel(x)
      x{i} = average(x{i}, seldim-1);
    end % for
  end
else
  x = average(x, seldim);
end
end % function cellmatmean
                                                                                                                                                                                             n = fieldnames(data);
fd = nan(size(fn));
for i=1:numel(fn)
  fd(i) = ndims(data.(fn{i}));
end

if ~isfield(data, 'dim')
  % this part depends on the assumption that the list of positions is describing a full 3D volume in
  % an ordered way which allows for the extraction of a transformation matrix, i.e. slice by slice
  data.dim = pos2dim(data.pos);
  try
    % if the dim is correct, it should be possible to obtain the transform
    ws = warning('off', 'MATLAB:rankDeficientMatrix');
    pos2transform(data.pos, data.dim);
    warning(ws);
  catch
    % remove the incorrect dim
    data = rmfield(data, 'dim');
  end
end

if isfield(data, 'dim')
  data.transform = pos2transform(data.pos, data.dim);
end

% remove the unwanted fields
data = removefields(data, {'pos', 'xgrid', 'ygrid', 'zgrid', 'tri', 'tet', 'hex'});

% make inside a volume
data = fixinside(data, 'logical');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = freq2raw(freq)

if isfield(freq, 'powspctrm')
  param = 'powspctrm';
elseif isfield(freq, 'fourierspctrm')
  param = 'fourierspctrm';
else
  ft_error('not supported for this data representation');
end

if strcmp(freq.dimord, 'rpt_chan_freq_time') || strcmp(freq.dimord, 'rpttap_chan_freq_time')
  dat = freq.(param);
elseif strcmp(freq.dimord, 'chan_freq_time')
  dat = freq.(param);
  dat = reshape(dat, [1 size(dat)]); % add a singleton dimension
else
  ft_error('not supported for dimord %s', freq.dimord);
end

nrpt  = size(dat,1);
nchan = size(dat,2);
nfreq = size(dat,3);
ntime = size(dat,4);
data = [];
% create the channel labels like "MLP11@12Hz""
k = 0;
for i=1:nfreq
  for j=1:nchan
    k = k+1;
    data.label{k} = sprintf('%s@%dHz', freq.label{j}, freq.freq(i));
  end
end
% reshape and copy the data as if it were timecourses only
for i=1:nrpt
  data.time{i}  = freq.time;
  data.trial{i} = reshape(dat(i,:,:,:), nchan*nfreq, ntime);
  if any(sum(isnan(data.trial{i}),1)==size(data.trial{i},1))
    tmp = sum(~isfinite(data.trial{i}),1)==size(data.trial{i},1);
    begsmp = find(~tmp,1, 'first');
    endsmp = find(~tmp,1, 'last' );
    data.trial{i} = data.trial{i}(:, begsmp:endsmp);
    data.time{i}  = data.time{i}(begsmp:endsmp);
  end
end

if isfield(freq, 'trialinfo'), data.trialinfo = freq.trialinfo; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tlck] = raw2timelock(data)

data   = ft_checkdata(data, 'hassampleinfo', 'yes');
ntrial = numel(data.trial);
nchan  = numel(data.label);

if ntrial==1
  tlck.time   = data.time{1};
  tlck.avg    = data.trial{1};
  tlck.label  = data.label;
  tlck.dimord = 'chan_time';
  tlck        = copyfields(data, tlck, {'grad', 'elec', 'opto', 'cfg', 'trialinfo', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord'});
  
else
  % the code below tries to construct a general time-axis where samples of all trials can fall on
  % find the earliest beginning and latest ending
  begtime = min(cellfun(@min, data.time));
  endtime = max(cellfun(@max, data.time));
  % find 'common' sampling rate
  fsample = 1./nanmean(cellfun(@mean, cellfun(@diff,data.time, 'uniformoutput', false)));
  % estimate number of samples
  nsmp = round((endtime-begtime)*fsample) + 1; % numerical round-off issues should be dealt with by this round, as they will/should never cause an extra sample to appear
  % construct general time-axis
  time = linspace(begtime,endtime,nsmp);
  
  % concatenate all trials
  tmptrial = nan(ntrial, nchan, length(time));
  
  begsmp = nan(ntrial, 1);
  endsmp = nan(ntrial, 1);
  for i=1:ntrial
    begsmp(i) = nearest(time, data.time{i}(1));
    endsmp(i) = nearest(time, data.time{i}(end));
    tmptrial(i,:,begsmp(i):endsmp(i)) = data.trial{i};
  end
  
  % construct the output timelocked data
  tlck.trial   = tmptrial;
  tlck.time    = time;
  tlck.dimord  = 'rpt_chan_time';
  tlck.label   = data.label;
  tlck         = copyfields(data, tlck, {'grad', 'elec', 'opto', 'cfg', 'trialinfo', 'sampleinfo', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord'});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = timelock2raw(data)
fn = getdatfield(data);
if any(ismember(fn, {'trial', 'individual', 'avg'}))
  % trial, individual and avg (in that order) should be preferred over all other data fields
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2965#c12
  fn = fn(ismember(fn, {'trial', 'individual', 'avg'}));
end
dimord = cell(size(fn));
for i=1:numel(fn)
  % determine the dimensions of each of the data fields
  dimord{i} = getdimord(data, fn{i});
end
% the fields trial, individual and avg (with their corresponding default dimord) are preferred
if sum(strcmp(dimord, 'rpt_chan_time'))==1
  fn = fn{strcmp(dimord, 'rpt_chan_time')};
  ft_info('constructing trials from "%s"\n', fn);
  dimsiz = getdimsiz(data, fn);
  ntrial = dimsiz(1);
  nchan  = dimsiz(2);
  ntime  = dimsiz(3);
  tmptrial = {};
  tmptime  = {};
  for j=1:ntrial
    tmptrial{j} = reshape(data.(fn)(j,:,:), [nchan, ntime]);
    tmptime{j}  = data.time;
  end
  data       = rmfield(data, fn);
  data.trial = tmptrial;
  data.time  = tmptime;
elseif sum(strcmp(dimord, 'subj_chan_time'))==1
  fn = fn{strcmp(dimord, 'subj_chan_time')};
  ft_info('constructing trials from "%s"\n', fn);
  dimsiz = getdimsiz(data, fn);
  nsubj = dimsiz(1);
  nchan  = dimsiz(2);
  ntime  = dimsiz(3);
  tmptrial = {};
  tmptime  = {};
  for j=1:nsubj
    tmptrial{j} = reshape(data.(fn)(j,:,:), [nchan, ntime]);
    tmptime{j}  = data.time;
  end
  data       = rmfield(data, fn);
  data.trial = tmptrial;
  data.time  = tmptime;
elseif sum(strcmp(dimord, 'chan_time'))==1
  fn = fn{strcmp(dimord, 'chan_time')};
  ft_info('constructing single trial from "%s"\n', fn);
  data.time  = {data.time};
  data.trial = {data.(fn)};
  data = rmfield(data, fn);
else
  ft_error('unsupported data structure');
end
% remove unwanted fields
data = removefields(data, {'avg', 'var', 'cov', 'dimord', 'numsamples' ,'dof'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = chan2freq(data)
data.dimord = [data.dimord '_freq'];
data.freq   = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = chan2timelock(data)
data.dimord = [data.dimord '_time'];
data.time   = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike] = raw2spike(data)
ft_info('converting raw data into spike data\n');
nTrials 	 = length(data.trial);
[spikelabel] = detectspikechan(data);
spikesel     = match_str(data.label, spikelabel);
nUnits       = length(spikesel);
if nUnits==0
  ft_error('cannot convert raw data to spike format since the raw data structure does not contain spike channels');
end

trialTimes  = zeros(nTrials,2);
for iUnit = 1:nUnits
  unitIndx = spikesel(iUnit);
  spikeTimes  = []; % we dont know how large it will be, so use concatenation inside loop
  trialInds   = [];
  for iTrial = 1:nTrials
    
    % read in the spike times
    [spikeTimesTrial]    = getspiketimes(data, iTrial, unitIndx);
    nSpikes              = length(spikeTimesTrial);
    spikeTimes           = [spikeTimes; spikeTimesTrial(:)];
    trialInds            = [trialInds; ones(nSpikes,1)*iTrial];
    
    % get the begs and ends of trials
    hasNum = find(~isnan(data.time{iTrial}));
    if iUnit==1, trialTimes(iTrial,:) = data.time{iTrial}([hasNum(1) hasNum(end)]); end
  end
  
  spike.label{iUnit}     = data.label{unitIndx};
  spike.waveform{iUnit}  = [];
  spike.time{iUnit}      = spikeTimes(:)';
  spike.trial{iUnit}     = trialInds(:)';
  
  if iUnit==1, spike.trialtime             = trialTimes; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = spike2raw(spike, fsample)

if nargin<2 || isempty(fsample)
  timeDiff = abs(diff(sort([spike.time{:}])));
  fsample  = 1/min(timeDiff(timeDiff>0));
  ft_warning('Desired sampling rate for spike data not specified, automatically resampled to %f', fsample);
end

% get some sizes
nUnits  = length(spike.label);
nTrials = size(spike.trialtime,1);

% preallocate
data.trial(1:nTrials) = {[]};
data.time(1:nTrials)  = {[]};
for iTrial = 1:nTrials
  
  % make bins: note that the spike.time is already within spike.trialtime
  x = [spike.trialtime(iTrial,1):(1/fsample):spike.trialtime(iTrial,2)];
  timeBins   = [x x(end)+1/fsample] - (0.5/fsample);
  time       = (spike.trialtime(iTrial,1):(1/fsample):spike.trialtime(iTrial,2));
  
  % convert to continuous
  trialData = zeros(nUnits,length(time));
  for iUnit = 1:nUnits
    
    % get the timestamps and only select those timestamps that are in the trial
    ts       = spike.time{iUnit};
    hasTrial = spike.trial{iUnit}==iTrial;
    ts       = ts(hasTrial);
    
    N = histc(ts,timeBins);
    if isempty(N)
      N = zeros(1,length(timeBins)-1);
    else
      N(end) = [];
    end
    
    % store it in a matrix
    trialData(iUnit,:) = N;
  end
  
  data.trial{iTrial} = trialData;
  data.time{iTrial}  = time;
  
end % for all trials

% create the associated labels and other aspects of data such as the header
data.label = spike.label;
data.fsample = fsample;
if isfield(spike,'hdr'), data.hdr = spike.hdr; end
if isfield(spike,'cfg'), data.cfg = spike.cfg; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = source2raw(source)

fn = fieldnames(source);
fn = setdiff(fn, {'pos', 'dim', 'transform', 'time', 'freq', 'cfg'});
for i=1:length(fn)
  dimord{i} = getdimord(source, fn{i});
end
sel = strcmp(dimord, 'pos_time');
assert(sum(sel)>0, 'the source structure does not contain a suitable field to represent as raw channel-level data');
assert(sum(sel)<2, 'the source structure contains multiple fields that can be represented as raw channel-level data');
fn     = fn{sel};
dimord = dimord{sel};

switch dimord
  case 'pos_time'
    % add fake raw channel data to the original data structure
    data.trial{1} = source.(fn);
    data.time{1}  = source.time;
    % add fake channel labels
    data.label = {};
    for i=1:size(source.pos,1)
      data.label{i} = sprintf('source%d', i);
    end
    data.label = data.label(:);
    data.cfg = source.cfg;
  otherwise
    % FIXME other formats could be implemented as well
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for detection of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikelabel, eeglabel] = detectspikechan(data)

maxRate = 2000; % default on what we still consider a neuronal signal: this firing rate should never be exceeded

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    hasAllInts    = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:) == round(data.trial{i}(j,:)));
    hasAllPosInts = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:)>=0);
    T = nansum(diff(data.time{i}),2); % total time
    fr            = nansum(data.trial{i}(j,:),2) ./ T;
    spikechan(j)  = spikechan(j) + double(hasAllInts & hasAllPosInts & fr<=maxRate);
  end
end
spikechan = (spikechan==ntrial);

spikelabel = data.label(spikechan);
eeglabel   = data.label(~spikechan);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikeTimes] = getspiketimes(data, trial, unit)
spikeIndx       = logical(data.trial{trial}(unit,:));
spikeCount      = data.trial{trial}(unit,spikeIndx);
spikeTimes      = data.time{trial}(spikeIndx);
if isempty(spikeTimes), return; end
multiSpikes     = find(spikeCount>1);
% get the additional samples and spike times, we need only loop through the bins
[addSamples, addTimes]   = deal([]);
for iBin = multiSpikes(:)' % looping over row vector
  addTimes     = [addTimes ones(1,spikeCount(iBin))*spikeTimes(iBin)];
  addSamples   = [addSamples ones(1,spikeCount(iBin))*spikeIndx(iBin)];
end
% before adding these times, first remove the old ones
spikeTimes(multiSpikes) = [];
spikeTimes              = sort([spikeTimes(:); addTimes(:)]);
                                                                                                                    'biosemi128'
      label = {
        'A1'
        'A2'
        'A3'
        'A4'
        'A5'
        'A6'
        'A7'
        'A8'
        'A9'
        'A10'
        'A11'
        'A12'
        'A13'
        'A14'
        'A15'
        'A16'
        'A17'
        'A18'
        'A19'
        'A20'
        'A21'
        'A22'
        'A23'
        'A24'
        'A25'
        'A26'
        'A27'
        'A28'
        'A29'
        'A30'
        'A31'
        'A32'
        'B1'
        'B2'
        'B3'
        'B4'
        'B5'
        'B6'
        'B7'
        'B8'
        'B9'
        'B10'
        'B11'
        'B12'
        'B13'
        'B14'
        'B15'
        'B16'
        'B17'
        'B18'
        'B19'
        'B20'
        'B21'
        'B22'
        'B23'
        'B24'
        'B25'
        'B26'
        'B27'
        'B28'
        'B29'
        'B30'
        'B31'
        'B32'
        'C1'
        'C2'
        'C3'
        'C4'
        'C5'
        'C6'
        'C7'
        'C8'
        'C9'
        'C10'
        'C11'
        'C12'
        'C13'
        'C14'
        'C15'
        'C16'
        'C17'
        'C18'
        'C19'
        'C20'
        'C21'
        'C22'
        'C23'
        'C24'
        'C25'
        'C26'
        'C27'
        'C28'
        'C29'
        'C30'
        'C31'
        'C32'
        'D1'
        'D2'
        'D3'
        'D4'
        'D5'
        'D6'
        'D7'
        'D8'
        'D9'
        'D10'
        'D11'
        'D12'
        'D13'
        'D14'
        'D15'
        'D16'
        'D17'
        'D18'
        'D19'
        'D20'
        'D21'
        'D22'
        'D23'
        'D24'
        'D25'
        'D26'
        'D27'
        'D28'
        'D29'
        'D30'
        'D31'
        'D32'
        };
      
    case 'biosemi256'
      label = {
        'A1'
        'A2'
        'A3'
        'A4'
        'A5'
        'A6'
        'A7'
        'A8'
        'A9'
        'A10'
        'A11'
        'A12'
        'A13'
        'A14'
        'A15'
     