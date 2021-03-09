function [dataout] = ft_denoise_dssp(cfg, datain)

% FT_DENOISE_DSSP implements a dual signal subspace projection algorithm
% to suppress interference outside a predefined source region of
% interest. It is based on: Sekihara et al. J. Neural Eng. 2016 13(3), and
% Sekihara et al. J. Neural Eng. 2018 15(3).
%
% Use as
%   dataout = ft_denoise_dssp(cfg, datain)
% where cfg is a configuration structure that contains
%   cfg.channel          = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.sourcemodel      = structure, source model with precomputed leadfields, see FT_PREPARE_LEADFIELD
%   cfg.dssp             = structure with parameters that determine the behavior of the algorithm
%   cfg.dssp.n_space     = 'all', or scalar. Number of dimensions for the
%                          initial spatial projection.
%   cfg.dssp.n_in        = 'all', or scalar. Number of dimensions of the
%                          subspace describing the field inside the ROI.
%   cfg.dssp.n_out       = 'all', or scalar. Number of dimensions of the
%                          subspace describing the field outside the ROI.
%   cfg.dssp.n_intersect = scalar (default = 0.9). Number of dimensions (if
%                          value is an integer>=1), or threshold for the
%                          included eigenvalues (if value<1), determining
%                          the dimensionality of the intersection.
%
% See also FT_DENOISE_PCA, FT_DENOISE_SYNTHETIC, FT_DENOISE_TSR

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% check the input data
datain = ft_checkdata(datain, 'datatype', {'raw'}); % FIXME how about timelock and freq?

cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'grid',    'sourcemodel'});

% get the options
cfg.trials       = ft_getopt(cfg, 'trials',  'all', 1);
cfg.channel      = ft_getopt(cfg, 'channel', 'all');
cfg.sourcemodel  = ft_getopt(cfg, 'sourcemodel');
cfg.dssp         = ft_getopt(cfg, 'dssp');         % sub-structure to hold the parameters
cfg.dssp.n_space = ft_getopt(cfg.dssp, 'n_space', 'all'); % number of spatial components to retain from the Gram matrix
cfg.dssp.n_in    = ft_getopt(cfg.dssp, 'n_in', 'all');    % dimensionality of the Bin subspace to be used for the computation of the intersection
cfg.dssp.n_out   = ft_getopt(cfg.dssp, 'n_out', 'all');   % dimensionality of the Bout subspace to be used for the computation of the intersection
cfg.dssp.n_intersect = ft_getopt(cfg.dssp, 'n_intersect', 0.9); % dimensionality of the intersection
cfg.output       = ft_getopt(cfg, 'output', 'original');

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});
datain = ft_selectdata(tmpcfg, datain);
% restore the provenance information
[cfg, datain] = rollback_provenance(cfg, datain);

% match the input data's channels with the labels in the leadfield
sourcemodel = cfg.sourcemodel;
if ~isfield(sourcemodel, 'leadfield')
  ft_error('cfg.sourcemodel needs to contain leadfields');
end
[indx1, indx2] = match_str(datain.label, sourcemodel.label);
if ~isequal(indx1(:),(1:numel(datain.label))')
  ft_error('unsupported mismatch between data channels and leadfields');
end
if islogical(sourcemodel.inside)
  inside = find(sourcemodel.inside);
else
  inside = sourcemodel.inside;
end
for k = inside(:)'
  sourcemodel.leadfield{k} = sourcemodel.leadfield{k}(indx2,:);
end

% compute the Gram-matrix of the supplied forward model
lf = cat(2, sourcemodel.leadfield{:});
G  = lf*lf';

dat     = cat(2,datain.trial{:});
[dum, Ae, N, Nspace, Sout, Sin, Sspace, S] = dssp(dat, G, cfg.dssp.n_in, cfg.dssp.n_out, cfg.dssp.n_space, cfg.dssp.n_intersect);
datAe   = dat*Ae; % the projection is a right multiplication
% with a matrix (eye(size(Ae,1))-Ae*Ae'), since Ae*Ae' can become quite
% sizeable, it's computed slightly differently here.

% put some diagnostic information in the output cfg.
cfg.dssp.S_space        = Sspace;
cfg.dssp.n_space        = Nspace;
cfg.dssp.S_out          = Sout;
cfg.dssp.S_in           = Sin;
cfg.dssp.S_intersect    = S;
cfg.dssp.n_intersect    = N;

% compute the cleaned data and put in a cell-array
nsmp  = cellfun(@numel, datain.time);
csmp  = cumsum([0 nsmp]);
trial = cell(size(datain.trial));
switch cfg.output
  case 'original'
    for k = 1:numel(datain.trial)
      trial{k} = datain.trial{k} - datAe*Ae((csmp(k)+1):csmp(k+1),:)';
    end
  case 'complement'
    for k = 1:numel(datain.trial)
      trial{k} = datAe*Ae((csmp(k)+1):csmp(k+1),:)';
    end
  otherwise
    ft_error(sprintf('cfg.output = ''%s'' is not implemented',cfg.output));
end

% create the output argument
dataout       = keepfields(datain, {'label','time','fsample','trialinfo','sampleinfo','grad', 'elec', 'opto'}); % grad can be kept and does not need to be balanced, since the cleaned data is a mixture over time, not space.
dataout.trial = trial;

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions for the computation of the projection matrix
% kindly provided by Kensuke, and adjusted a bit by Jan-Mathijs
function [Bclean, Ae, Nee, Nspace, Sout, Sin, Sspace, S] = dssp(B, G, Nin, Nout, Nspace, Nee)

% Nc: number of sensors
% Nt: number of time points
% inputs
% B(Nc,Nt):  interference overlapped sensor data
% G(Nc,Nc): Gram matrix of voxel lead field
% Nout and Nin: dimensions of the two row spaces
% recom_Nspace: recommended value for the dimension of the pseudo-signal subspace
% outputs
% Bclean(Nc,Nt): cleaned sensor data
% Nee: dimension of the intersection
% Nspace: dimension of the pseudo-signal subspace
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%
% The code below is modified by Jan-Mathijs, no functional changes
% merely cosmetics

% eigen decomposition of the Gram matrix, matrix describing the spatial
% components
[U,S]   = eig(G);
Sspace  = abs(diag(S));

[Sspace, iorder] = sort(-Sspace);
Sspace           = -Sspace;
U(:,:)           = U(:,iorder);

if isempty(Nspace)
  ttext = 'enter the spatial dimension: ';
  Nspace    = input(ttext);
elseif ischar(Nspace) && isequal(Nspace, 'interactive')
  figure, plot(log10(Sspace),'-o');
  Nspace = input('enter spatial dimension of the ROI subspace: ');
elseif ischar(Nspace) && isequal(Nspace, 'all')
  Nspace = find(Sspace./Sspace(1)>1e5*eps, 1, 'last');
end
fprintf('Using %d spatial dimensions\n', Nspace);

% spatial subspace projector
Us   = U(:,1:Nspace);
USUS = Us*Us';

% Bin and Bout creations
Bin  =                  USUS  * B;
Bout = (eye(size(USUS))-USUS) * B;

% create the temporal subspace projector and apply it to the data
%[AeAe, Nee] = CSP01(Bin, Bout, Nout, Nin, Nee);
%Bclean      = B*(eye(size(AeAe))-AeAe);

[Ae, Nee, Sout, Sin, S] = CSP01(Bin, Bout, Nin, Nout, Nee);
Bclean    = B - (B*Ae)*Ae'; % avoid computation of Ae*Ae'


function [Ae, Nee, Sout, Sin, S] = CSP01(Bin, Bout, Nin, Nout, Nee)
%
% interference rejection by removing the common temporal subspace of the two subspaces
% K. Sekihara,  March 28, 2012
% Golub and Van Loan, Matrix computations, The Johns Hopkins University Press, 1996
%
%  Nc: number of channels
%  Nt: number of time points
% inputs
%  Bout(1:Nc,1:Nt): interference data
%  Bin(1:Nc,1:Nt): signal plus interference data
%  ypost(1:Nc,1:Nt): denoised data
%  Nout: dimension of the interference subspace
%  Nin: dimension of the signal plus interference subspace
%  Nee: dimension of the intersection of the two subspaces
% outputs
% Ae = matrix from which the projector onto the intersection can
%      be obtained:
% AeAe: projector onto the intersection, which is equal to the
%       interference subspace.
% Nee: dimension of the intersection
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%

[dum,Sout,Vout] = svd(Bout,'econ');
[dum,Sin, Vin]  = svd(Bin, 'econ');
Sout = diag(Sout);
Sin  = diag(Sin);

if isempty(Nout)
  ttext = 'enter the spatial dimension for the outside field: ';
  Nout  = input(ttext);
elseif ischar(Nout) && isequal(Nout, 'interactive')
  figure, plot(Sout,'-o');
  Nout = input('enter dimension of the outside field: ');
elseif ischar(Nout) && isequal(Nout, 'all')
  Nout = find(Sout./Sout(1)>1e5*eps, 1, 'last');
end
fprintf('Using %d spatial dimensions for the outside field\n', Nout);

if isempty(Nin)
  ttext = 'enter the spatial dimension for the inside field: ';
  Nin  = input(ttext);
elseif ischar(Nin) && isequal(Nin, 'interactive')
  figure, plot(log10(Sin),'-o');
  Nin = input('enter dimension of the inside field: ');
elseif ischar(Nin) && isequal(Nin, 'all')
  Nin = find(Sin./Sin(1)>1e5*eps, 1, 'last');
end
fprintf('Using %d spatial dimensions for the inside field\n', Nin);

Qout = Vout(:,1:Nout);
Qin  = Vin(:, 1:Nin);

C     = Qin'*Qout;
[U,S] = svd(C);
S     = diag(S);
if (ischar(Nee) && strcmp(Nee, 'auto'))
  ft_error('automatic determination of intersection dimension is not supported');
elseif ischar(Nee) && strcmp(Nee, 'interactive')
  figure, plot(S,'-o');
  Nee  = input('enter dimension of the intersection: ');
elseif Nee<1
  % treat a numeric value < 1 as a threshold
  Nee = find(S>Nee,1,'last');
  if isempty(Nee), Nee = 1; end
end
fprintf('Using %d dimensions for the interaction\n', Nee);

Ae   = Qin*U;
Ae   = Ae(:,1:Nee);
%AeAe = Ae*Ae';
                                                                                                                                                                                                                                                                                             exchange', 3, 1); % not required
  end

  try
    % these directories deal with compatibility with older MATLAB versions
    if ft_platform_supports('matlabversion', -inf, '2008a'), ft_hastoolbox('compat/matlablt2008b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2008b'), ft_hastoolbox('compat/matlablt2009a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2009a'), ft_hastoolbox('compat/matlablt2009b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2009b'), ft_hastoolbox('compat/matlablt2010a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2010a'), ft_hastoolbox('compat/matlablt2010b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2010b'), ft_hastoolbox('compat/matlablt2011a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2011a'), ft_hastoolbox('compat/matlablt2011b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2011b'), ft_hastoolbox('compat/matlablt2012a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2012a'), ft_hastoolbox('compat/matlablt2012b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2012b'), ft_hastoolbox('compat/matlablt2013a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2013a'), ft_hastoolbox('compat/matlablt2013b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2013b'), ft_hastoolbox('compat/matlablt2014a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2014a'), ft_hastoolbox('compat/matlablt2014b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2014d'), ft_hastoolbox('compat/matlablt2015a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2015a'), ft_hastoolbox('compat/matlablt2015b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2015b'), ft_hastoolbox('compat/matlablt2016a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2016a'), ft_hastoolbox('compat/matlablt2016b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2016b'), ft_hastoolbox('compat/matlablt2017a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2017a'), ft_hastoolbox('compat/matlablt2017b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2017b'), ft_hastoolbox('compat/matlablt2018a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2018a'), ft_hastoolbox('compat/matlablt2018b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2018b'), ft_hastoolbox('compat/matlablt2019a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2019a'), ft_hastoolbox('compat/matlablt2019b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2019b'), ft_hastoolbox('compat/matlablt2020a', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2020a'), ft_hastoolbox('compat/matlablt2020b', 3, 1); end
    if ft_platform_supports('matlabversion', -inf, '2020b'), ft_hastoolbox('compat/matlablt2021a', 3, 1); end
    % this deals with compatibility with all OCTAVE versions
    if ft_platform_supports('octaveversion', -inf, +inf),    ft_hastoolbox('compat/octave', 3, 1); end
  end

  try
    % these contains template layouts, neighbour structures, MRIs and cortical meshes
    ft_hastoolbox('template/layout',      1, 1);
    ft_hastoolbox('template/anatomy',     1, 1);
    ft_hastoolbox('template/headmodel',   1, 1);
    ft_hastoolbox('template/electrode',   1, 1);
    ft_hastoolbox('template/neighbours',  1, 1);
    ft_hastoolbox('template/sourcemodel', 1, 1);
  end

  try
    % this is used in ft_statistics
    ft_hastoolbox('statfun', 1, 1);
  end

  try
    % this is used in ft_definetrial
    ft_hastoolbox('trialfun', 1, 1);
  end

  try
    % this contains the low-level reading functions
    ft_hastoolbox('fileio', 1, 1);
  end

  try
    % this is for filtering etc. on time-series data
    ft_hastoolbox('preproc', 1, 1);
  end

  try
    % this contains forward models for the EEG and MEG volume conductor
    ft_hastoolbox('forward', 1, 1);
  end

  try
    % this contains inverse source estimation methods
    ft_hastoolbox('inverse', 1, 1);
  end

  try
    % this contains intermediate-level plotting functions, e.g. multiplots and 3-d objects
    ft_hastoolbox('plotting', 1, 1);
  end

  try
    % this contains intermediate-level functions for spectral analysis
    ft_hastoolbox('specest', 1, 1);
  end

  try
    % this contains the functions to compute connectivity metrics
    ft_hastoolbox('connectivity', 1, 1);
  end

  try
    % this contains test scripts
    ft_hastoolbox('test', 1, 1);
  end

  try
    % this contains the functions for spike and spike-field analysis
    ft_hastoolbox('contrib/spike', 1, 1);
  end

  try
    % this contains user contributed functions
    ft_hastoolbox('contrib/misc', 1, 1);
  end

  try
    % this contains specific code and examples for realtime processing
    ft_hastoolbox('realtime/example', 3, 1);    % not required
    ft_hastoolbox('realtime/online_mri', 3, 1); % not required
    ft_hastoolbox('realtime/online_meg', 3, 1); % not required
    ft_hastoolbox('realtime/online_eeg', 3, 1); % not required
  end

end

% the toolboxes added by this function should not be removed by FT_POSTAMBLE_HASTOOLBOX
ft_default.toolbox.cleanup = prevcleanup;

% track the usage of this function, this only happens once at startup
ft_trackusage('startup');

% remember that the function has executed in a persistent variable
initialized = true;

end % function ft_default


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkMultipleToolbox(toolbox, keyfile)

persistent warned
if isempty(warned)
  warned = false;
end

if ~ft_platform_supports('which-all')
  return;
end

list = which(keyfile, '-all');
if length(list)>1
  ft_warning('Multiple versions of %s on your path will confuse FieldTrip', toolbox);
  if ~warned % only throw the following warnings once
    warned = true;
    for i=1:length(list)
      ft_warning('one version of %s is found here: %s', toolbox, list{i});
    end
  end
  ft_warning('You probably used addpath(genpath(''path_to_fieldtrip'')), this can lead to unexpected behavior. See http://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path');
end
end % function checkMultipleToolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkIncorrectPath
p = fileparts(mfilename('fullpath'));
incorrect = fullfile(p, 'compat', 'incorrect');
if ~isempty(strfind(path, incorrect))
  ft_warning('Your path is set up incorrectly. You probably used addpath(genpath(''path_to_fieldtrip'')), this can lead to unexpected behavior. See http://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path');
end
end % function checkIncorrectPath
                                                                                                                             'normalized', 'style', 'text', 'string', 'translate', 'callback', [])
uicontrol('tag', 'tx', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.translate(1)), 'callback', @cb_redraw)
uicontrol('tag', 'ty', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.translate(2)), 'callback', @cb_redraw)
uicontrol('tag', 'tz', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.translate(3)), 'callback', @cb_redraw)
ft_uilayout(figHandle, 'tag', 'translateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'tx',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'ty',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'tz',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);

% somehow the toolbar gets lost in 2012b
set(figHandle, 'toolbar', 'figure');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_close(figHandle, varargin)
% the figure will be closed in the main function after collecting the guidata
uiresume;
                                                                                                                                                                                                                                                                                                                                                                                                   t = ft_read_event(cfg.dataset);
  end
  
  cfg.channel = ft_channelselection(cfg.channel, hdr.label);
  chansel = match_str(hdr.label, cfg.channel);
  Nchans  = length(chansel);
  
  if strcmp(cfg.continuous, 'yes')
    Ntrials = 1;
  else
    Ntrials = hdr.nTrials;
  end
  
  % construct trl-matrix for the data on disk
  if ~isempty(cfg.trl)
    % take the one specified by the user
    trlorg = cfg.trl;
  else
    % make one that is according to the segments/trials in the data on disk
    trlorg = zeros(Ntrials,3);
    if strcmp(cfg.continuous, 'yes')
      trlorg(1, [1 2]) = [1 hdr.nSamples*hdr.nTrials];
    else
      for k = 1:Ntrials
        trlorg(k, [1 2]) = [1 hdr.nSamples] + [hdr.nSamples hdr.nSamples] .* (k-1);
      end
    end
  end
end % if hasdata

if strcmp(cfg.continuous, 'no') && isempty(cfg.blocksize)
  cfg.blocksize = (trlorg(1,2) - trlorg(1,1)+1) ./ hdr.Fs;
elseif strcmp(cfg.continuous, 'yes') && isempty(cfg.blocksize)
  cfg.blocksize = 1;
end

if cfg.blocksize<round(10*1/hdr.Fs)
  ft_warning('the blocksize is very small given the samping rate, increasing blocksize to 10 samples');
  cfg.blocksize = round(10*1/hdr.Fs);
end

% FIXME make a check for the consistency of cfg.continuous, cfg.blocksize, cfg.trl and the data header

if Nchans == 0
  ft_error('no channels to display');
end

if Ntrials == 0
  ft_error('no trials to display');
end

if ischar(cfg.selectfeature)
  % ensure that it is a cell-array
  cfg.selectfeature = {cfg.selectfeature};
end
if ~isempty(cfg.selectfeature)
  for i=1:length(cfg.selectfeature)
    if ~isfield(cfg.artfctdef, cfg.selectfeature{i})
      cfg.artfctdef.(cfg.selectfeature{i})          = [];
      cfg.artfctdef.(cfg.selectfeature{i}).artifact = zeros(0,2);
    end
  end
end

% determine the vertical scaling
if ischar(cfg.ylim)
  if hasdata
    sel = 1;
    while all(isnan(reshape(data.trial{sel}(chansel,:),[],1)))
      sel = sel+1;
    end
    % the first trial is used to determine the vertical scaling
    dat = data.trial{sel}(chansel,:);
  else
    % one second of data is read from file to determine the vertical scaling
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', 1, 'endsample', round(hdr.Fs), 'chanindx', chansel, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat, 'headerformat', cfg.headerformat);
  end % if hasdata
  % convert the data to another numeric precision, i.e. double, single or int32
  if ~isempty(cfg.precision)
    dat = cast(dat, cfg.precision);
  end
  minval = min(dat(:));
  maxval = max(dat(:));
  switch cfg.ylim
    case 'maxabs'
      maxabs   = max(abs([minval maxval]));
      scalefac = 10^(fix(log10(maxabs)));
      if scalefac==0
        % this happens if the data is all zeros
        scalefac=1;
      end
      maxabs   = (round(maxabs / scalefac * 100) / 100) * scalefac;
      cfg.ylim = [-maxabs maxabs];
    case 'maxmin'
      if minval==maxval
        % this happens if the data is constant, e.g. all zero or clipping
        minval = minval - eps;
        maxval = maxval + eps;
      end
      cfg.ylim = [minval maxval];
    otherwise
      ft_error('unsupported value for cfg.ylim');
  end % switch ylim
  % zoom in a bit when viemode is vertical
  if strcmp(cfg.viewmode, 'vertical')
    cfg.ylim = cfg.ylim/10;
  end
else
  if (numel(cfg.ylim) ~= 2) || ~isnumeric(cfg.ylim)
    ft_error('cfg.ylim needs to be a 1x2 vector [ymin ymax], describing the upper and lower limits')
  end
end

% determine the coloring of channels
if hasdata
  linecolor = linecolor_common(cfg, data);
else
  linecolor = linecolor_common(cfg, hdr);
end

% collect the artifacts that have been detected from cfg.artfctdef.xxx.artifact
artlabel = fieldnames(cfg.artfctdef);
sel      = zeros(size(artlabel));
artifact = cell(size(artlabel));
for i=1:length(artlabel)
  sel(i) = isfield(cfg.artfctdef.(artlabel{i}), 'artifact');
  if sel(i)
    artifact{i} = cfg.artfctdef.(artlabel{i}).artifact;
    ft_info('detected %3d %s artifacts\n', size(artifact{i}, 1), artlabel{i});
  end
end

% get the subset of the artfctdef fields that seem to contain artifacts
artifact = artifact(sel==1);
artlabel = artlabel(sel==1);

if length(artlabel) > 9
  ft_error('only up to 9 artifacts groups supported')
end

% make artdata representing all artifacts in a "raw data" format
datendsample = max(trlorg(:,2));

artdata = [];
artdata.trial{1}       = convert_event(artifact, 'boolvec', 'endsample', datendsample); % every artifact is a "channel"
artdata.time{1}        = offset2time(0, hdr.Fs, datendsample);
artdata.label          = artlabel;
artdata.fsample        = hdr.Fs;
artdata.cfg.trl        = [1 datendsample 0];

% determine amount of unique event types (for cfg.ploteventlabels)
if ~isempty(event) && isstruct(event)
  eventtypes = unique({event.type});
else
  eventtypes = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up default functions to be available in the right-click menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cfg.selfun - labels that are presented in rightclick menu, and is appended using ft_getuserfun(..., 'browse') later on to create a function handle
% cfg.selcfg - cfgs for functions to be executed

if ~isempty(cfg.selfun) || ~isempty(cfg.selcfg)
  if ischar(cfg.selfun)
    cfg.selfun = {cfg.selfun};
  end
  if isstruct(cfg.selcfg)
    cfg.selcfg = {cfg.selcfg};
  end
elseif isempty(cfg.selfun) && isempty(cfg.selcfg)
  % simplefft
  cfg.selcfg{1} = [];
  cfg.selcfg{1}.linecolor = linecolor;
  cfg.selfun{1} = 'simpleFFT';
  % multiplotER
  cfg.selcfg{2} = [];
  cfg.selcfg{2}.linecolor = linecolor;
  cfg.selcfg{2}.layout = cfg.layout;
  cfg.selfun{2} = 'multiplotER';
  % topoplotER
  cfg.selcfg{3} = [];
  cfg.selcfg{3}.linecolor = linecolor;
  cfg.selcfg{3}.layout = cfg.layout;
  cfg.selfun{3} = 'topoplotER';
  % topoplotVAR
  cfg.selcfg{4} = [];
  cfg.selcfg{4}.layout = cfg.layout;
  cfg.selfun{4} = 'topoplotVAR';
  % movieplotER
  cfg.selcfg{5} = [];
  cfg.selcfg{5}.layout = cfg.layout;
  cfg.selcfg{5}.interactive = 'yes';
  cfg.selfun{5} = 'movieplotER';
  % audiovideo
  cfg.selcfg{6} = [];
  cfg.selcfg{6}.audiofile = ft_getopt(cfg, 'audiofile');
  cfg.selcfg{6}.videofile = ft_getopt(cfg, 'videofile');
  cfg.selcfg{6}.anonimize = ft_getopt(cfg, 'anonimize');
  cfg.selfun{6} = 'audiovideo';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the data structures used in the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% opt represents the global data/settings, it should contain
% - the original data, epoched or continuous
% - the artifacts represented as continuous data
% - the redraw_cb settings
% - the preproc   settings
% - the select_range_cb settings (also used in keyboard_cb)

% these elements are stored inside the figure, so that the callback routines can modify them
opt = [];
if hasdata
  opt.orgdata   = data;
else
  opt.orgdata   = [];      % this means that it will read from cfg.dataset
end
if strcmp(cfg.continuous, 'yes')
  opt.trialviewtype = 'segment';
else
  opt.trialviewtype = 'trial';
end
opt.artdata     = artdata;
opt.hdr         = hdr;
opt.event       = event;
opt.trlop       = 1;          % the active trial being displayed
opt.ftsel       = find(strcmp(artlabel, cfg.selectfeature)); % current artifact/feature being selected
opt.trlorg      = trlorg;
opt.fsample     = hdr.Fs;
opt.artifactcolors = [0.9686 0.7608 0.7686; 0.7529 0.7098 0.9647; 0.7373 0.9725 0.6824; 0.8118 0.8118 0.8118; 0.9725 0.6745 0.4784; 0.9765 0.9176 0.5686; 0.6863 1 1; 1 0.6863 1; 0 1 0.6000];
opt.linecolor   = linecolor;
opt.cleanup     = false;      % this is needed for a corrent handling if the figure is closed (either in the corner or by "q")
opt.chanindx    = [];         % this is used to check whether the component topographies need to be redrawn
opt.eventtypes  = eventtypes;
opt.eventtypescolors = [0 0 0; 1 0 0; 0 0 1; 0 1 0; 1 0 1; 0.5 0.5 0.5; 0 1 1; 1 1 0];
opt.eventtypecolorlabels = {'black', 'red', 'blue', 'green', 'cyan', 'grey', 'light blue', 'yellow'};
opt.nanpaddata  = []; % this is used to allow horizontal scaling to be constant (when looking at last segment continuous data, or when looking at segmented/zoomed-out non-continuous data)
opt.trllock     = []; % this is used when zooming into trial based data


% save original layout when viewmode = component
if strcmp(cfg.viewmode, 'component')
  opt.layorg = cfg.layout;
end

% determine labelling of channels
if strcmp(cfg.plotlabels, 'yes')
  opt.plotLabelFlag = 1;
elseif strcmp(cfg.plotlabels, 'some')
  opt.plotLabelFlag = 2;
else
  opt.plotLabelFlag = 0;
end

% set changedchanflg as true for initialization
opt.changedchanflg = true; % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)

% create fig
if isfield(cfg, 'position')
  h = figure('Position', cfg.position);
else
  h = figure;
end

% put appdata in figure
setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);
if ~isempty(cfg.renderer)
  set(h, 'renderer', cfg.renderer);
end

% set interruptible to off, see bug 3123
set(h, 'Interruptible', 'off', 'BusyAction', 'queue'); % enforce busyaction to queue to be sure

% enable custom data cursor text
dcm = datacursormode(h);
set(dcm, 'updatefcn', @cb_datacursortext);

% set the figure window title
funcname = mfilename();
if ~hasdata
  if isfield(cfg, 'dataset')
    dataname = cfg.dataset;
  elseif isfield(cfg, 'datafile')
    dataname = cfg.datafile;
  else
    dataname = [];
  end
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  dataname = cfg.inputfile;
else
  dataname = inputname(2);
end
set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), funcname, join_str(', ',dataname)));
set(gcf, 'NumberTitle', 'off');

% set zoom option to on
% zoom(h, 'on')
% set(zoom(h), 'actionPostCallback', @zoom_drawlabels_cb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the figure and callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(h, 'KeyPressFcn',           @keyboard_cb);
set(h, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonDownFcn'});
set(h, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonUpFcn'});
set(h, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonMotionFcn'});
if any(strcmp(cfg.viewmode, {'component', 'vertical'}))
  set(h, 'ReSizeFcn',           @winresize_cb); % resize will now trigger redraw and replotting of labels
end

% make the user interface elements for the data view
uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', opt.trialviewtype, 'userdata', 't')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'leftarrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'rightarrow')

if strcmp(cfg.viewmode, 'component')
  uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'component', 'userdata', 'c')
else
  uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'channel', 'userdata', 'c')
end

uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'uparrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'downarrow')

uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'horizontal', 'userdata', 'h')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+leftarrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+rightarrow')

uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'vertical', 'userdata', 'v')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+downarrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+uparrow')

% legend artifacts/features
for iArt = 1:length(artlabel)
  uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', artlabel{iArt}, 'userdata', num2str(iArt), 'position', [0.91, 0.9 - ((iArt-1)*0.09), 0.08, 0.04], 'backgroundcolor', opt.artifactcolors(iArt,:))
  uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', ['shift+' num2str(iArt)], 'position', [0.91, 0.855 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artifactcolors(iArt,:))
  uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', ['control+' num2str(iArt)], 'position', [0.96, 0.855 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artifactcolors(iArt,:))
end
if length(artlabel)>1 % highlight the first one as active
  arth = findobj(h, 'tag', 'artifactui');
  arth = arth(end:-1:1); % order is reversed so reverse it again
  hsel = [1 2 3] + (opt.ftsel-1) .*3;
  set(arth(hsel), 'fontweight', 'bold')
end

if true % strcmp(cfg.viewmode, 'butterfly')
  % button to find label of nearest channel to datapoint
  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'identify', 'userdata', 'i', 'position', [0.91, 0.1, 0.08, 0.05], 'backgroundcolor', [1 1 1])
end

% 'edit preproc'-button
uicontrol('tag', 'preproccfg', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'preproc cfg', 'position', [0.91, 0.55 - ((iArt-1)*0.09), 0.08, 0.04], 'callback', @preproc_cfg1_cb)

ft_uilayout(h, 'tag', 'labels',  'width', 0.10, 'height', 0.05);
ft_uilayout(h, 'tag', 'buttons', 'width', 0.05, 'height', 0.05);

ft_uilayout(h, 'tag', 'labels',     'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'buttons',    'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'artifactui', 'style', 'pushbutton', 'callback', @keyboard_cb);

ft_uilayout(h, 'tag', 'labels',  'retag', 'viewui');
ft_uilayout(h, 'tag', 'buttons', 'retag', 'viewui');
ft_uilayout(h, 'tag', 'viewui', 'BackgroundColor', [0.8 0.8 0.8], 'hpos', 'auto', 'vpos', 0);

definetrial_cb(h);
redraw_cb(h);

% %% Scrollbar
%
% % set initial scrollbar value
% dx = maxtime;
%
% % set scrollbar position
% fig_pos=get(gca, 'position');
% scroll_pos=[fig_pos(1) fig_pos(2) fig_pos(3) 0.02];
%
% % define callback
% S=['set(gca, ''xlim'',get(gcbo, ''value'')+[ ' num2str(mintime) ', ' num2str(maxtime) '])'];
%
% % Creating Uicontrol
% s=uicontrol('style', 'slider',...
%     'units', 'normalized', 'position',scroll_pos,...
%     'callback',S, 'min',0, 'max',0, ...
%     'visible', 'off'); %'value', xmin

% set initial scrollbar value
% dx = maxtime;
%
% % set scrollbar position
% fig_pos=get(gca, 'position');
% scroll_pos=[fig_pos(1) fig_pos(2) fig_pos(3) 0.02];
%
% % define callback
% S=['set(gca, ''xlim'',get(gcbo, ''value'')+[ ' num2str(mintime) ', ' num2str(maxtime) '])'];
%
% % Creating Uicontrol
% s=uicontrol('style', 'slider',...
%     'units', 'normalized', 'position',scroll_pos,...
%     'callback',S, 'min',0, 'max',0, ...
%     'visible', 'off'); %'value', xmin
%initialize postion of plot
% set(gca, 'xlim', [xmin xmin+dx]);

if nargout
  % wait until the user interface is closed, get the user data with the updated artifact details
  set(h, 'CloseRequestFcn', @cleanup_cb);
  
  while ishandle(h)
    uiwait(h);
    opt = getappdata(h, 'opt');
    if opt.cleanup
      delete(h);
    end
  end
  
  % add the updated artifact definitions to the output cfg
  for i=1:length(opt.artdata.label)
    cfg.artfctdef.(opt.artdata.label{i}).artifact = convert_event(opt.artdata.trial{1}(i,:), 'artifact');
  end
  
  % add the updated preproc to the output
  try
    browsecfg = getappdata(h, 'cfg');
    cfg.preproc = browsecfg.preproc;
  end
  
  % add the update event to the output cfg
  cfg.event = opt.event;
  
end % if nargout

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance

if ~nargout
  clear cfg
end

end % main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(h, eventdata)
opt = getappdata(h, 'opt');
opt.cleanup = true;
setappdata(h, 'opt', opt);
uiresume
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function definetrial_cb(h, eventdata)
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

if strcmp(cfg.continuous, 'no')
  
  % when zooming in, lock the trial! one can only go to the next trial when horizontal scaling doesn't segment the data - from ft-meeting: this might be relaxed later on - roevdmei
  if isempty(opt.trllock)
    opt.trllock = opt.trlop;
  end
  locktrllen = ((opt.trlorg(opt.trllock,2)-opt.trlorg(opt.trllock,1)+1) ./ opt.fsample);
  % if cfg.blocksize is close to the length of the locked trial, set it to that
  if (abs(locktrllen-cfg.blocksize) / locktrllen) < 0.1
    cfg.blocksize = locktrllen;
  end
  
  %%%%%%%%%
  % trial is locked, change subdivision of trial
  if cfg.blocksize < locktrllen
    % lock the trial if it wasn't locked (and thus trlop refers to the actual trial)
    if isempty(opt.trllock)
      opt.trllock = trlop;
    end
    % save current position if already
    if isfield(opt, 'trlvis')
      thissegbeg = opt.trlvis(opt.trlop,1);
    end
    datbegsample = min(opt.trlorg(opt.trllock,1));
    datendsample = max(opt.trlorg(opt.trllock,2));
    smpperseg  = round(opt.fsample * cfg.blocksize);
    begsamples = datbegsample:smpperseg:datendsample;
    endsamples = datbegsample+smpperseg-1:smpperseg:datendsample;
    offset     = (((1:numel(begsamples))-1)*smpperseg) + opt.trlorg(opt.trllock,3);
    if numel(endsamples)<numel(begsamples)
      endsamples(end+1) = datendsample;
    end
    trlvis = [];
    trlvis(:,1) = begsamples';
    trlvis(:,2) = endsamples';
    trlvis(:,3) = offset;
    % determine length of each trial, and determine the offset with the current requested zoom-level
    trllen   = (trlvis(:,2) - trlvis(:,1)+1);
    sizediff = smpperseg - trllen;
    opt.nanpaddata = sizediff;
    
    if isfield(opt, 'trlvis')
      % update the current trial counter and try to keep the current sample the same
      opt.trlop   = nearest(begsamples, thissegbeg);
    end
    % update trialviewtype
    opt.trialviewtype = 'trialsegment';
    % update button
    set(findobj(get(h, 'children'), 'string', 'trial'), 'string', 'segment');
    %%%%%%%%%
    
    
    %%%%%%%%%
    % trial is not locked, go to original trial division and zoom out
  elseif cfg.blocksize >= locktrllen
    trlvis = opt.trlorg;
    % set current trlop to locked trial if it was locked before
    if ~isempty(opt.trllock)
      opt.trlop = opt.trllock;
    end
    smpperseg  = round(opt.fsample * cfg.blocksize);
    % determine length of each trial, and determine the offset with the current requested zoom-level
    trllen   = (trlvis(:,2) - trlvis(:,1)+1);
    sizediff = smpperseg - trllen;
    opt.nanpaddata = sizediff;
    
    % update trialviewtype
    opt.trialviewtype = 'trial';
    % update button
    set(findobj(get(h, 'children'), 'string', 'trialsegment'), 'string',opt.trialviewtype);
    
    % release trial lock
    opt.trllock = [];
    %%%%%%%%%
  end
  
  % save trlvis
  opt.trlvis  = trlvis;
  
else
  % construct a trial definition for visualisation
  if isfield(opt, 'trlvis') % if present, remember where we were
    thistrlbeg = opt.trlvis(opt.trlop,1);
  end
  % look at cfg.blocksize and make opt.trl accordingly
  datbegsample = min(opt.trlorg(:,1));
  datendsample = max(opt.trlorg(:,2));
  smpperseg  = round(opt.fsample * cfg.blocksize);
  begsamples = datbegsample:smpperseg:datendsample;
  endsamples = datbegsample+smpperseg-1:smpperseg:datendsample;
  if numel(endsamples)<numel(begsamples)
    endsamples(end+1) = datendsample;
  end
  trlvis = [];
  trlvis(:,1) = begsamples';
  trlvis(:,2) = endsamples';
  % compute the offset. In case if opt.trlorg has multiple trials, the first sample is t=0, otherwise use the offset in opt.trlorg
  if size(opt.trlorg,1)==1
    offset = begsamples - repmat(begsamples(1), [1 numel(begsamples)]); % offset for all segments compared to the first
    offset = offset + opt.trlorg(1,3);
    trlvis(:,3) = offset;
  else
    offset = begsamples - repmat(begsamples(1), [1 numel(begsamples)]);
    trlvis(:,3) = offset;
  end
  
  if isfield(opt, 'trlvis')
    % update the current trial counter and try to keep the current sample the same
    % opt.trlop   = nearest(round((begsamples+endsamples)/2), thissample);
    opt.trlop   = nearest(begsamples, thistrlbeg);
  end
  opt.trlvis  = trlvis;
  
  % NaN-padding when horizontal scaling is bigger than the data
  % two possible situations, 1) zoomed out so far that all data is one segment, or 2) multiple segments but last segment is smaller than the rest
  sizediff = smpperseg-(endsamples-begsamples+1);
  opt.nanpaddata = sizediff;
end % if continuous
setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function help_cb(h, eventdata)
fprintf('------------------------------------------------------------------------------------\n')
fprintf('You can use the following keyboard buttons in the databrowser\n')
fprintf('1-9                : select artifact type 1-9\n');
fprintf('shift 1-9          : select previous artifact of type 1-9\n');
fprintf('                     (does not work with numpad keys)\n');
fprintf('control 1-9        : select next artifact of type 1-9\n');
fprintf('alt 1-9            : select next artifact of type 1-9\n');
fprintf('arrow-left         : previous trial\n');
fprintf('arrow-right        : next trial\n');
fprintf('shift arrow-up     : increase vertical scaling\n');
fprintf('shift arrow-down   : decrease vertical scaling\n');
fprintf('shift arrow-left   : increase horizontal scaling\n');
fprintf('shift arrow-down   : decrease horizontal scaling\n');
fprintf('s                  : toggles between cfg.selectmode options\n');
fprintf('q                  : quit\n');
fprintf('------------------------------------------------------------------------------------\n')
fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_range_cb(h, range, cmenulab) %range 1X4 in sec relative to current trial
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

% the range should be in the displayed box
range(1) = max(opt.hpos-opt.width/2, range(1));
range(2) = max(opt.hpos-opt.width/2, range(2));
range(1) = min(opt.hpos+opt.width/2, range(1));
range(2) = min(opt.hpos+opt.width/2, range(2));
range = (range-(opt.hpos-opt.width/2)) / opt.width; % left side of the box becomes 0, right side becomes 1
range = range * (opt.hlim(2) - opt.hlim(1)) + opt.hlim(1);   % 0 becomes hlim(1), 1 becomes hlim(2)

begsample = opt.trlvis(opt.trlop,1);
endsample = opt.trlvis(opt.trlop,2);
offset    = opt.trlvis(opt.trlop,3);

% determine the selection
begsel = round(range(1)*opt.fsample+begsample-offset-1);
endsel = round(range(2)*opt.fsample+begsample-offset);
% artifact selection is now always based on begsample/endsample/offset
% -roevdmei

% the selection should always be confined to the current trial
begsel = max(begsample, begsel);
endsel = min(endsample, endsel);

% mark or execute selfun
if isempty(cmenulab)
  % the left button was clicked INSIDE a selected range, update the artifact definition or event
  
  if strcmp(cfg.selectmode, 'markartifact')
    % mark or unmark artifacts
    artval = opt.artdata.trial{1}(opt.ftsel, begsel:endsel);
    artval = any(artval,1);
    if any(artval)
      fprintf('there is overlap with the active artifact (%s), disabling this artifact\n',opt.artdata.label{opt.ftsel});
      opt.artdata.trial{1}(opt.ftsel, begsel:endsel) = 0;
    else
      fprintf('there is no overlap with the active artifact (%s), marking this as a new artifact\n',opt.artdata.label{opt.ftsel});
      opt.artdata.trial{1}(opt.ftsel, begsel:endsel) = 1;
    end
    
    % redraw only when marking (so the focus doesn't go back to the databrowser after calling selfuns
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h);
  elseif strcmp(cfg.selectmode, 'markpeakevent') || strcmp(cfg.selectmode, 'marktroughevent')
    %mark or unmark events, marking at peak/trough of window
    if any(intersect(begsel:endsel, [opt.event.sample]))
      fprintf('there is overlap with one or more event(s), disabling this/these event(s)\n');
      ind_rem = intersect(begsel:endsel, [opt.event.sample]);
      for iRemove = 1:length(ind_rem)
        opt.event([opt.event.sample]==ind_rem(iRemove)) = [];
      end
    else
      fprintf('there is no overlap with any event, adding an event to the peak/trough value\n');
      % check if only 1 chan, other wise not clear max in which channel. %
      % ingnie: would be cool to add the option to select the channel when multiple channels
      if size(dat,1) > 1
        ft_error('cfg.selectmode = ''markpeakevent'' and ''marktroughevent'' only supported with 1 channel in the data')
      end
      if strcmp(cfg.selectmode, 'markpeakevent')
        [dum ind_minmax] = max(dat(begsel-begsample+1:endsel-begsample+1));
        val = 'peak';
      elseif strcmp(cfg.selectmode, 'marktroughevent')
        [dum ind_minmax] = min(dat(begsel-begsample+1:endsel-begsample+1));
        val = 'trough';
      end
      samp_minmax = begsel + ind_minmax - 1;
      event_new.type     = 'ft_databrowser_manual';
      event_new.sample   = samp_minmax;
      event_new.value    = val;
      event_new.duration = 1;
      event_new.offset   = 0;
      % add new event to end opt.event
      % check if events are in order now
      if  min(diff([opt.event.sample]))>0
        % add new event in line with old ones
        nearest_event = nearest([opt.event.sample], samp_minmax);
        if opt.event(nearest_event).sample > samp_minmax
          %place new event before nearest
          ind_event_new = nearest_event;
        else
          %place new event after nearest
          ind_event_new = nearest_event +1;
        end
        event_lastpart = opt.event(ind_event_new:end);
        opt.event(ind_event_new) = event_new;
        opt.event(ind_event_new+1:end+1) = event_lastpart;
      else
        %just add to end
        opt.event(end+1) = event_new;
      end
      clear event_new ind_event_new event_lastpart val dum ind_minmax
    end
    % redraw only when marking (so the focus doesn't go back to the databrowser after calling selfuns
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h);
  end
  
else
  % the right button was used to activate the context menu and the user made a selection from that menu
  % execute the corresponding function
  
  % get index into cfgs
  selfunind = strcmp(cfg.selfun, cmenulab);
  
  % cut out the requested data segment
  switch cfg.seldat
    case 'current'
      seldata             = keepfields(opt.curdata, {'label', 'grad', 'elec', 'opto', 'hdr'});
      seldata.trial{1}    = ft_fetch_data(opt.curdata, 'begsample', begsel, 'endsample', endsel);
    case 'all'
      seldata             = keepfields(opt.org, {'label', 'grad', 'elec', 'opto', 'hdr'});
      seldata.trial{1}    = ft_fetch_data(opt.orgdata, 'begsample', begsel, 'endsample', endsel);
  end
  seldata.time{1}     = offset2time(offset+begsel-begsample, opt.fsample, endsel-begsel+1);
  seldata.fsample     = opt.fsample;
  seldata.sampleinfo  = [begsel endsel];
  
  % prepare input
  funhandle = ft_getuserfun(cmenulab, 'browse');
  funcfg    = cfg.selcfg{selfunind};
  % get windowname and give as input (can be used for the other functions as well, not implemented yet)
  if ~strcmp(opt.trialviewtype, 'trialsegment')
    str = sprintf('%s %d/%d, time from %g to %g s', opt.trialviewtype, opt.trlop, size(opt.trlvis,1), seldata.time{1}(1), seldata.time{1}(end));
  else
    str = sprintf('trial %d/%d: segment: %d/%d , time from %g to %g s', opt.trllock, size(opt.trlorg,1), opt.trlop, size(opt.trlvis,1), seldata.time{1}(1), seldata.time{1}(end));
  end
  funcfg.figurename = [cmenulab ': ' str];
  feval(funhandle, funcfg, seldata);
end

end % function select_range_cb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preproc_cfg1_cb(h,eventdata)
parent = get(h, 'parent');
cfg = getappdata(parent, 'cfg');

editfontsize = cfg.editfontsize;
editfontunits = cfg.editfontunits;

% parse cfg.preproc
if ~isempty(cfg.preproc)
  code = printstruct('cfg.preproc', cfg.preproc);
else
  code = '';
end

% add descriptive lines
code = {
  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  '% Add or change options for on-the-fly preprocessing'
  '% Use as cfg.preproc.option=value'
  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  code
  };

% make figure displaying the edit box
pph = figure;
axis off
% add save button
uicontrol('tag', 'preproccfg_l2', 'parent', pph, 'units', 'normalized', 'style', 'pushbutton', 'string', 'save and close', 'position', [0.81, 0.6 , 0.18, 0.10], 'callback', @preproc_cfg2_cb);

% add edit box
ppeh = uicontrol('style', 'edit');
set(pph, 'toolBar', 'none')
set(pph, 'menuBar', 'none')
set(pph, 'Name', 'cfg.preproc editor')
set(pph, 'NumberTitle', 'off')
set(ppeh, 'Units', 'normalized');
set(ppeh, 'Position', [0 0 .8 1]);
set(ppeh, 'backgroundColor', [1 1 1]);
set(ppeh, 'horizontalAlign', 'left');
set(ppeh, 'max', 2);
set(ppeh, 'min', 0);
set(ppeh, 'FontName', 'Courier', 'FontSize', editfontsize, 'FontUnits', editfontunits);
set(ppeh, 'string', code);

% add handle for the edit style to figure
setappdata(pph, 'superparent', parent); % superparent is the main ft_databrowser window
setappdata(pph, 'ppeh', ppeh);

end % function preproc_cfg1_cb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preproc_cfg2_cb(h,eventdata)
parent = get(h, 'parent');
superparent = getappdata(parent, 'superparent');
ppeh = getappdata(parent, 'ppeh');
code = get(ppeh, 'string');

% get rid of empty lines and white space
skip = [];
for iline = 1:numel(code)
  code{iline} = strtrim(code{iline});
  if isempty(code{iline})
    skip = [skip iline];
    continue
  end
  if code{iline}(1)=='%'
    skip = [skip iline];
    continue
  end
end
code(skip) = [];

if ~isempty(code)
  ispreproccfg = strncmp(code, 'cfg.preproc.',12);
  if ~all(ispreproccfg)
    errordlg('cfg-options must be specified as cfg.preproc.xxx', 'cfg.preproc editor', 'modal')
  end
  % eval the code
  for icomm = 1:numel(code)
    eval([code{icomm} ';']);
  end
  
  % check for cfg and output into the original appdata-window
  if ~exist('cfg', 'var')
    cfg = [];
    cfg.preproc = [];
  end
  maincfg = getappdata(superparent, 'cfg');
  maincfg.preproc = cfg.preproc;
  setappdata(superparent, 'cfg', maincfg)
end

close(parent)
redraw_cb(superparent)
end % function preproc_cfg2_cb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyboard_cb(h, eventdata)

if (isempty(eventdata) && ft_platform_supports('matlabversion',-Inf, '2014a')) || isa(eventdata, 'matlab.ui.eventdata.ActionData')
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

switch key
  case {'1' '2' '3' '4' '5' '6' '7' '8' '9'}
    % switch to another artifact type
    opt.ftsel = str2double(key);
    numart = size(opt.artdata.trial{1}, 1);
    if opt.ftsel > numart
      fprintf('data has no artifact type %i \n', opt.ftsel)
    else
      % bold the active one
      arth = findobj(h, 'tag', 'artifactui');
      arth = arth(end:-1:1); % order is reversed so reverse it again
      hsel = [1 2 3] + (opt.ftsel-1) .*3 ;
      set(arth(hsel), 'fontweight', 'bold')
      % unbold the passive ones
      set(arth(setdiff(1:numel(arth),hsel)), 'fontweight', 'normal')
      % redraw
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      fprintf('switching to the "%s" artifact\n', opt.artdata.label{opt.ftsel});
      redraw_cb(h, eventdata);
    end
  case {'shift+1' 'shift+2' 'shift+3' 'shift+4' 'shift+5' 'shift+6' 'shift+7' 'shift+8' 'shift+9'}
    % go to previous artifact
    opt.ftsel = str2double(key(end));
    numart = size(opt.artdata.trial{1}, 1);
    if opt.ftsel > numart
      fprintf('data has no artifact type %i \n', opt.ftsel)
    else
      % find the previous occuring artifact, keeping in mind that:
      % 1) artifacts can cross trial boundaries
      % 2) artifacts might not occur inside a trial boundary (when data is segmented differently than during artifact detection)
      % fetch trl representation of current artifact type
      arttrl = convert_event(opt.artdata.trial{1}(opt.ftsel,:), 'trl');
      % discard artifacts in the future
      curvisend = opt.trlvis(opt.trlop,2);
      arttrl(arttrl(:,1) > curvisend,:) = [];
      % find nearest artifact by searching in each trl (we have to do this here everytime, because trlvis can change on the fly because of x-zooming)
      newtrlop = [];
      for itrlvis = opt.trlop-1:-1:1
        % is either the start or the end of any artifact present?
        if any(any(opt.trlvis(itrlvis,1)<=arttrl(:,1:2) & opt.trlvis(itrlvis,2)>=arttrl(:,1:2)))
          % if so, we're done
          newtrlop = itrlvis;
          break
        end
      end
      if isempty(newtrlop)
        fprintf('no earlier %s with "%s" artifact found\n', opt.trialviewtype, opt.artdata.label{opt.ftsel});
      else
        fprintf('going to previous %s with "%s" artifact\n', opt.trialviewtype, opt.artdata.label{opt.ftsel});
        opt.trlop = newtrlop;
        % other artifact type potentially selected, bold the active one
        arth = findobj(h, 'tag', 'artifactui');
        arth = arth(end:-1:1); % order is reversed so reverse it again
        hsel = [1 2 3] + (opt.ftsel-1) .*3 ;
        set(arth(hsel), 'fontweight', 'bold')
        % unbold the passive ones
        set(arth(setdiff(1:numel(arth),hsel)), 'fontweight', 'normal')
        % export into fig and continue
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
      end
    end
  case {'control+1' 'control+2' 'control+3' 'control+4' 'control+5' 'control+6' 'control+7' 'control+8' 'control+9' 'alt+1' 'alt+2' 'alt+3' 'alt+4' 'alt+5' 'alt+6' 'alt+7' 'alt+8' 'alt+9'}
    % go to next artifact
    opt.ftsel = str2double(key(end));
    numart = size(opt.artdata.trial{1}, 1);
    if opt.ftsel > numart
      fprintf('data has no artifact type %i \n', opt.ftsel)
    else
      % find the next occuring artifact, keeping in mind that:
      % 1) artifacts can cross trial boundaries
      % 2) artifacts might not occur inside a trial boundary (when data is segmented differently than during artifact detection)
      % fetch trl representation of current artifact type
      arttrl = convert_event(opt.artdata.trial{1}(opt.ftsel,:), 'trl');
      % discard artifacts in the past
      curvisbeg = opt.trlvis(opt.trlop,1);
      arttrl(arttrl(:,2) < curvisbeg,:) = [];
      % find nearest artifact by searching in each trl (we have to do this here everytime, because trlvis can change on the fly because of x-zooming)
      newtrlop = [];
      for itrlvis = opt.trlop+1:size(opt.trlvis,1)
        % is either the start or the end of any artifact present?
        if any(any(opt.trlvis(itrlvis,1)<=arttrl(:,1:2) & opt.trlvis(itrlvis,2)>=arttrl(:,1:2)))
          % if so, we're done
          newtrlop = itrlvis;
          break
        end
      end
      if isempty(newtrlop)
        fprintf('no later %s with "%s" artifact found\n', opt.trialviewtype, opt.artdata.label{opt.ftsel});
      else
        fprintf('going to next %s with "%s" artifact\n', opt.trialviewtype, opt.artdata.label{opt.ftsel});
        opt.trlop = newtrlop;
        % other artifact type potentially selected, bold the active one
        arth = findobj(h, 'tag', 'artifactui');
        arth = arth(end:-1:1); % order is reversed so reverse it again
        hsel = [1 2 3] + (opt.ftsel-1) .*3 ;
        set(arth(hsel), 'fontweight', 'bold')
        % unbold the passive ones
        set(arth(setdiff(1:numel(arth),hsel)), 'fontweight', 'normal')
        % export into fig and continue
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
      end
    end
  case 'leftarrow'
    opt.trlop = max(opt.trlop - 1, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'rightarrow'
    opt.trlop = min(opt.trlop + 1, size(opt.trlvis,1)); % should not be larger than the number of trials
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'uparrow'
    chansel = match_str(opt.hdr.label, cfg.channel);
    minchan = min(chansel);
    numchan = length(chansel);
    chansel = minchan - numchan : minchan - 1;
    if min(chansel)<1
      chansel = chansel - min(chansel) + 1;
    end
    % convert numeric array into cell-array with channel labels
    cfg.channel = opt.hdr.label(chansel);
    opt.changedchanflg = true; % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h, 'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
    redraw_cb(h, eventdata);
  case 'downarrow'
    chansel = match_str(opt.hdr.label, cfg.channel);
    maxchan = max(chansel);
    numchan = length(chansel);
    chansel = maxchan + 1 : maxchan + numchan;
    if max(chansel)>length(opt.hdr.label)
      chansel = chansel - (max(chansel) - length(opt.hdr.label));
    end
    % convert numeric array into cell-array with channel labels
    cfg.channel = opt.hdr.label(chansel);
    opt.changedchanflg = true; % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h, 'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
    redraw_cb(h, eventdata);
  case 'shift+leftarrow'
    cfg.blocksize = cfg.blocksize*sqrt(2);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    definetrial_cb(h, eventdata);
    redraw_cb(h, eventdata);
  case 'shift+rightarrow'
    cfg.blocksize = cfg.blocksize/sqrt(2);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    definetrial_cb(h, eventdata);
    redraw_cb(h, eventdata);
  case 'shift+uparrow'
    cfg.ylim = cfg.ylim/sqrt(2);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'shift+downarrow'
    cfg.ylim = cfg.ylim*sqrt(2);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'q'
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    cleanup_cb(h);
  case 't'
    % select the trial to display
    if ~strcmp(opt.trialviewtype, 'trialsegment')
      str = sprintf('%s to display (current trial = %d/%d)', opt.trialviewtype, opt.trlop, size(opt.trlvis,1));
    else
      str = sprintf('segment to display (current segment = %d/%d)', opt.trlop, size(opt.trlvis,1));
    end
    response = inputdlg(str, 'specify', 1, {num2str(opt.trlop)});
    if ~isempty(response)
      opt.trlop = str2double(response);
      opt.trlop = min(opt.trlop, size(opt.trlvis,1)); % should not be larger than the number of trials
      opt.trlop = max(opt.trlop, 1); % should not be smaller than 1
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      redraw_cb(h, eventdata);
    end
  case 'h'
    % select the horizontal scaling
    response = inputdlg('horizontal scale', 'specify', 1, {num2str(cfg.blocksize)});
    if ~isempty(response)
      cfg.blocksize = str2double(response);
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      definetrial_cb(h, eventdata);
      redraw_cb(h, eventdata);
    end
  case 'v'
    % select the vertical scaling
    response = inputdlg('vertical scale, [ymin ymax], ''maxabs'' or ''maxmin''', 'specify', 1, {['[ ' num2str(cfg.ylim) ' ]']});
    if ~isempty(response)
      if strcmp(response, 'maxmin')
        minval = min(opt.curdata.trial{1}(:));
        maxval = max(opt.curdata.trial{1}(:));
        cfg.ylim = [minval maxval];
      elseif strcmp(response, 'maxabs')
        minval = min(opt.curdata.trial{1}(:));
        maxval = max(opt.curdata.trial{1}(:));
        cfg.ylim = [-max(abs([minval maxval])) max(abs([minval maxval]))];
      else
        % convert to string and add brackets, just to ensure that str2num will work
        tmp = str2num(['[' response{1} ']']);
        if numel(tmp)==2
          cfg.ylim = tmp;
        else
          ft_warning('incorrect specification of cfg.ylim, not changing the limits for the vertical axes')
        end
      end
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      redraw_cb(h, eventdata);
    end
  case 'c'
    % select channels
    select = match_str(opt.hdr.label, cfg.channel);
    select = select_channel_list(opt.hdr.label, select);
    cfg.channel = opt.hdr.label(select);
    opt.changedchanflg = true; % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h, 'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
    redraw_cb(h, eventdata);
  case 'i'
    delete(findobj(h, 'tag', 'identify'));
    % click in data and get name of nearest channel
    fprintf('click in the figure to identify the name of the closest channel\n');
    val = ginput(1);
    pos = val(1);
    if strcmp(cfg.viewmode, 'butterfly') || strcmp(cfg.viewmode, 'vertical')
      switch cfg.viewmode
        case 'butterfly'
          % transform 'val' to match data
          val(1) = val(1) * range(opt.hlim) + opt.hlim(1);
          val(2) = val(2) * range(opt.vlim) + opt.vlim(1);
          channame = val2nearestchan(opt.curdata,val);
          channb   = match_str(opt.curdata.label,channame);
          % set chanposind
        