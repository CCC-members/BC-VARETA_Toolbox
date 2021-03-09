function [pnt, ori, lab] = channelposition(sens)

% CHANNELPOSITION computes the channel positions and orientations from the
% coils or electrodes
%
% Use as
%   [pos, ori, lab] = channelposition(sens)
% where sens is an electrode,  gradiometer or optode array.
%
% See also FT_DATATYPE_SENS

% Copyright (C) 2009-2018, Robert Oostenveld & Vladimir Litvak
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

% remove the balancing from the sensor definition, e.g. planar gradients, 3rd-order gradients, PCA-cleaned data or ICA projections
if isfield(sens, 'balance') && ~strcmp(sens.balance.current, 'none')
  sens = undobalancing(sens);
end

% keep it backward compatible with sensor definitions prior to 2011v1 (see ft_datatype_sens), which have pnt/ori instead of coilpos/coilori.
if isfield(sens, 'ori')
  sens.coilori = sens.ori;
  sens = rmfield(sens, 'ori');
end
if isfield(sens, 'pnt')
  sens.coilpos = sens.pnt;
  sens = rmfield(sens, 'pnt');
end

% treat all sensor arrays similar, i.e. as gradiometer systems
if     ~isfield(sens, 'coilori') && isfield(sens, 'coilpos')
  sens.coilori = nan(size(sens.coilpos));
elseif ~isfield(sens, 'coilori') && isfield(sens, 'elecpos')
  sens.coilori = nan(size(sens.elecpos));
elseif ~isfield(sens, 'coilori') && isfield(sens, 'optopos')
  sens.coilori = nan(size(sens.optopos));
end

switch ft_senstype(sens)
  case {'ctf64', 'ctf151', 'ctf275' 'bti148', 'bti248', 'bti248grad', 'itab28', 'itab153', 'yokogawa64', 'yokogawa160', 'babysquid74'}
    % the following code applies to systems with only axial gradiometers or magnetometers

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % do the MEG sensors first
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    sensorig   = sens;
    sel        = ft_chantype(sens, 'meg');
    sens.label = sens.label(sel);
    sens.tra   = sens.tra(sel,:);

    % subsequently remove the unused coils
    used = any(abs(sens.tra)>0.0001, 1);  % allow a little bit of rounding-off error
    sens.coilpos = sens.coilpos(used,:);
    sens.coilori = sens.coilori(used,:);
    sens.tra     = sens.tra(:,used);

    % compute distances from the center of the helmet
    center = mean(sens.coilpos);
    dist   = sqrt(sum((sens.coilpos - repmat(center, size(sens.coilpos, 1), 1)).^2, 2));

    % put the corresponding distances instead of non-zero tra entries
    maxval = repmat(max(abs(sens.tra),[],2), [1 size(sens.tra,2)]);
    maxval = min(maxval, ones(size(maxval))); %a value > 1 sometimes leads to problems; this is an empirical fix
    dist = (abs(sens.tra)>0.7.*maxval).*repmat(dist', size(sens.tra, 1), 1);

    % for the occasional case where there are nans: -> 0's will be
    % converted to inf anyhow
    dist(isnan(dist)) = 0;

    % put infs instead of the zero entries
    dist(~dist) = inf;

    % use the matrix to find coils with minimal distance to the center,
    % i.e. the bottom coil in the case of axial gradiometers
    % this only works for a full-rank unbalanced tra-matrix

    numcoils = sum(isfinite(dist),2);

    if all(numcoils==numcoils(1))
      % add the additional constraint that coils cannot be used twice,
      % i.e. for the position of 2 channels. A row of the dist matrix can end
      % up with more than 1 (magnetometer array) or 2 (axial gradiometer array)
      % non-zero entries when the input grad structure is rank-reduced
      % FIXME: I don't know whether this works for a vector-gradiometer
      % system. It also does not work when the system has mixed gradiometers
      % and magnetometers

      % use the magic that Jan-Mathijs implemented
      tmp      = mode(numcoils);
      niter    = 0;
      while ~all(numcoils==tmp)
        niter    = niter + 1;
        selmode  = find(numcoils==tmp);
        selrest  = setdiff((1:size(dist,1))', selmode);
        dist(selrest,sum(~isinf(dist(selmode,:)))>0) = inf;
        numcoils = sum(isfinite(dist),2);
        if niter>500
          ft_error('Failed to extract the positions of the channels. This is most likely due to the balancing matrix being rank deficient. Please replace data.grad with the original grad-structure obtained after reading the header.');
        end
      end
    else
      % assume that the solution is not so hard and just determine the bottom coil
    end

    [junk, ind] = min(dist, [], 2);

    lab(sel)   = sens.label;
    pnt(sel,:) = sens.coilpos(ind, :);
    ori(sel,:) = sens.coilori(ind, :);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % then do the references
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sens = sensorig;
    sel  = ft_chantype(sens, 'ref');

    sens.label = sens.label(sel);
    sens.tra   = sens.tra(sel,:);

    % subsequently remove the unused coils
    used = any(abs(sens.tra)>0.0001, 1);  % allow a little bit of rounding-off error
    sens.coilpos = sens.coilpos(used,:);
    sens.coilori = sens.coilori(used,:);
    sens.tra = sens.tra(:,used);

    [nchan, ncoil] = size(sens.tra);
    refpnt = zeros(nchan,3);
    refori = zeros(nchan,3); % FIXME not sure whether this will work
    for i=1:nchan
      weight = abs(sens.tra(i,:));
      weight = weight ./ sum(weight);
      refpnt(i,:) = weight * sens.coilpos;
      refori(i,:) = weight * sens.coilori;
    end
    reflab = sens.label;

    lab(sel)   = reflab;
    pnt(sel,:) = refpnt;
    ori(sel,:) = refori;

    sens = sensorig;

  case {'ctf64_planar', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar', 'bti248grad_planar', 'itab28_planar', 'itab153_planar', 'yokogawa64_planar', 'yokogawa160_planar'}
    % create a list with planar channel names
    chan = {};
    for i=1:length(sens.label)
      if endsWith(sens.label{i}, '_dH') || endsWith(sens.label{i}, '_dV')
        chan{i} = sens.label{i}(1:(end-3));
      end
    end
    chan = unique(chan);
    % find the matching channel-duplets
    ind = false(size(chan));
    lab = cell(length(chan),2);
    pnt = nan(length(chan),3);
    ori = nan(length(chan),3);
    for i=1:length(chan)
      ch1 =  [chan{i} '_dH'];
      ch2 =  [chan{i} '_dV'];
      sel = match_str(sens.label, {ch1, ch2});
      if length(sel)==2
        ind(i)   = true;
        lab(i,:) = {ch1, ch2};
        meanpnt1 = mean(sens.coilpos(abs(sens.tra(sel(1),:))>0.5, :), 1);
        meanpnt2 = mean(sens.coilpos(abs(sens.tra(sel(2),:))>0.5, :), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2], 1);
        meanori1 = mean(sens.coilori(abs(sens.tra(sel(1),:))>0.5, :), 1);
        meanori2 = mean(sens.coilori(abs(sens.tra(sel(2),:))>0.5, :), 1);
        ori(i,:) = mean([meanori1; meanori2], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    ori = ori(ind,:);

  case 'neuromag122'
    % find the matching channel-duplets
    ind = [];
    lab = {};
    for i=1:2:140
      % first try MEG channel labels with a space
      ch1 = sprintf('MEG %03d', i);
      ch2 = sprintf('MEG %03d', i+1);
      sel = match_str(sens.label, {ch1, ch2});
      % then try MEG channel labels without a space
      if (length(sel)~=2)
        ch1 = sprintf('MEG%03d', i);
        ch2 = sprintf('MEG%03d', i+1);
        sel = match_str(sens.label, {ch1, ch2});
      end
      % then try to determine the channel locations
      if (length(sel)==2)
        ind = [ind; i];
        lab(i,:) = {ch1, ch2};
        meanpnt1 = mean(sens.coilpos(abs(sens.tra(sel(1),:))>0.5,:), 1);
        meanpnt2 = mean(sens.coilpos(abs(sens.tra(sel(2),:))>0.5,:), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2], 1);
        meanori1 = mean(sens.coilori(abs(sens.tra(sel(1),:))>0.5,:), 1);
        meanori2 = mean(sens.coilori(abs(sens.tra(sel(2),:))>0.5,:), 1);
        ori(i,:) = mean([meanori1; meanori2], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    ori = ori(ind,:);

  case 'neuromag306'
    % find the matching channel-triplets
    ind = [];
    lab = {};
    for i=1:300
      % first try MEG channel labels with a space
      ch1 = sprintf('MEG %03d1', i);
      ch2 = sprintf('MEG %03d2', i);
      ch3 = sprintf('MEG %03d3', i);
      [sel1, sel2] = match_str(sens.label, {ch1, ch2, ch3});
      % the try MEG channels without a space
      if isempty(sel1)
        ch1 = sprintf('MEG%03d1', i);
        ch2 = sprintf('MEG%03d2', i);
        ch3 = sprintf('MEG%03d3', i);
        [sel1, sel2] = match_str(sens.label, {ch1, ch2, ch3});
      end
      % then try to determine the channel locations
      if (~isempty(sel1) && length(sel1)<=3)
        ind = [ind; i];
        lab(i,sel2) = sens.label(sel1)';
        meanpnt  = [];
        meanori  = [];
        for j = 1:length(sel1)
          meanpnt  = [meanpnt; mean(sens.coilpos(abs(sens.tra(sel1(j),:))>0.5,:), 1)];
          meanori  = [meanori; mean(sens.coilori(abs(sens.tra(sel1(j),:))>0.5,:), 1)];
        end
        pnt(i,:) = mean(meanpnt, 1);
        ori(i,:) = mean(meanori, 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    ori = ori(ind,:);

  otherwise
    % compute the position for each gradiometer or electrode
    nchan = length(sens.label);
    if isfield(sens, 'elecpos')
      nelec = size(sens.elecpos,1); % these are the electrodes
    elseif isfield(sens, 'coilpos')
      ncoil = size(sens.coilpos,1); % these are the coils
    elseif isfield(sens, 'optopos')
      nopto = size(sens.optopos,1); % these are the optodes
    end

    if ~isfield(sens, 'tra') && isfield(sens, 'elecpos') && nchan==nelec
      % there is one electrode per channel, which means that the channel position is identical to the electrode position
      pnt = sens.elecpos;
      if isfield(sens, 'elecori')
        ori = sens.elecori;
      else
        ori = nan(size(pnt));
      end
      lab = sens.label;

    elseif isfield(sens, 'tra') && isfield(sens, 'elecpos') && isequal(sens.tra, eye(nelec))
      % there is one electrode per channel, which means that the channel position is identical to the electrode position
      pnt = sens.elecpos;
      if isfield(sens, 'elecori')
        ori = sens.elecori;
      else
        ori = nan(size(pnt));
      end
      lab = sens.label;

    elseif isfield(sens, 'tra') && isfield(sens, 'elecpos') && isequal(sens.tra, eye(nelec)-1/nelec)
      % there is one electrode per channel, channels are average referenced
      pnt = sens.elecpos;
      if isfield(sens, 'elecori')
        ori = sens.elecori;
      else
        ori = nan(size(pnt));
      end
      lab = sens.label;

    elseif ~isfield(sens, 'tra') && isfield(sens, 'coilpos') && nchan==ncoil
      % there is one coil per channel, which means that the channel position is identical to the coil position
      pnt = sens.coilpos;
      ori = sens.coilori;
      lab = sens.label;

    elseif ~isfield(sens, 'tra') && isfield(sens, 'optopos') && nchan==nopto
      % there is one optode per channel, which means that the channel position is identical to the optode position
      pnt = sens.optopos;
      if isfield(sens, 'optoori')
        ori = sens.optoori;
      else
        ori = nan(size(pnt));
      end
      lab = sens.label;

    elseif isfield(sens, 'tra')
      % each channel depends on multiple sensors (electrodes or coils), compute a weighted position for the channel
      % for MEG gradiometer channels this means that the position is in between the two coils
      % for bipolar EEG channels this means that the position is in between the two electrodes
      % for NIRS channels this means that the position is in between the transmit and receive optode
      pnt = nan(nchan,3);
      ori = nan(nchan,3);
      if isfield(sens, 'coilpos')
        for i=1:nchan
          weight = abs(sens.tra(i,:));
          weight = weight ./ sum(weight);
          pnt(i,:) = weight * sens.coilpos;
          ori(i,:) = weight * sens.coilori;
        end
      elseif isfield(sens, 'elecpos')
        for i=1:nchan
          weight = abs(sens.tra(i,:));
          weight = weight ./ sum(weight);
          pnt(i,:) = weight * sens.elecpos;
        end
      elseif isfield(sens, 'optopos')
        for i=1:nchan
          weight = abs(sens.tra(i,:));
          weight = weight ./ sum(weight);
          pnt(i,:) = weight * sens.optopos;
        end
      end
      lab = sens.label;
    end
    
end % switch senstype

n = size(lab,2);
% this is to fix the planar layouts, which cannot be plotted anyway
if n>1 && size(lab, 1)>1 % this is to prevent confusion when lab happens to be a row array
  pnt = repmat(pnt, n, 1);
  ori = repmat(ori, n, 1);
end

% ensure that the channel order is the same as in sens
[sel1, sel2] = match_str(sens.label, lab);
lab = lab(sel2);
pnt = pnt(sel2, :);
ori = ori(sel2, :);

% ensure that it is a row vector
lab = lab(:);

% do a sanity check on the number of positions
nchan = numel(sens.label);
if length(lab)~=nchan || size(pnt,1)~=nchan || size(ori,1)~=nchan
  ft_warning('the positions were not determined for all channels');
end
                                                                                                                                                                                                                                                                                         axis                  - from patient anterior to posterior*
    % Slices in 'z' axis                  - from patient inferior to superior

    if verbose,
        fprintf('...writing axial flipped (+Y from Anterior to Posterior)\n');
    end

    avw.hdr.hist.orient = uint8(3);

    SliceDim = double(avw.hdr.dime.dim(4)); % z
    RowDim   = double(avw.hdr.dime.dim(3)); % y
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(4));
    RowSz    = double(avw.hdr.dime.pixdim(3));
    PixelSz  = double(avw.hdr.dime.pixdim(2));

    x = 1:PixelDim;
    for z = 1:SliceDim,
        for y = RowDim:-1:1, % flipped in Y
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end

case 4, % coronal flipped

    % For the 'coronal flipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'z' axis                  - from patient inferior to superior
    % Slices in 'y' axis                  - from patient anterior to posterior

    if verbose,
        fprintf('...writing coronal flipped (+Z from Superior to Inferior)\n');
    end

    avw.hdr.hist.orient = uint8(4);

    SliceDim = double(avw.hdr.dime.dim(3)); % y
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(3));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(2));

    x = 1:PixelDim;
    for y = 1:SliceDim,
        for z = RowDim:-1:1,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end

case 5, % sagittal flipped

    % For the 'sagittal flipped' type, the voxels are stored with
    % Pixels in 'y' axis (varies fastest) - from patient posterior to anterior
    % Rows in   'z' axis                  - from patient superior to inferior
    % Slices in 'x' axis                  - from patient right to left

    if verbose,
        fprintf('...writing sagittal flipped (+Z from Superior to Inferior)\n');
    end

    avw.hdr.hist.orient = uint8(5);

    SliceDim = double(avw.hdr.dime.dim(2)); % x
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(3)); % y
    SliceSz  = double(avw.hdr.dime.pixdim(2));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(3));

    y = 1:PixelDim;
    for x = 1:SliceDim,
        for z = RowDim:-1:1, % superior to inferior
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end

otherwise, % transverse/axial unflipped

    % For the 'transverse unflipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'y' axis                  - from patient posterior to anterior
    % Slices in 'z' axis                  - from patient inferior to superior

    if verbose,
        fprintf('...unknown orientation specified, assuming default axial unflipped\n');
    end

    avw.hdr.hist.orient = uint8(0);

    SliceDim = double(avw.hdr.dime.dim(4)); % z
    RowDim   = double(avw.hdr.dime.dim(3)); % y
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(4));
    RowSz    = double(avw.hdr.dime.pixdim(3));
    PixelSz  = double(avw.hdr.dime.pixdim(2));

    x = 1:PixelDim;
    for z = 1:SliceDim,
        for y = 1:RowDim,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end

end

fclose(fid);

% Update the header
avw.hdr.dime.dim(2:4) = int16([PixelDim,RowDim,SliceDim]);
avw.hdr.dime.pixdim(2:4) = single([PixelSz,RowSz,SliceSz]);

return







%----------------------------------------------------------------------------

function write_header(fid,avw,verbose)

    header_key(fid,avw.hdr.hk);
    image_dimension(fid,avw.hdr.dime);
    data_history(fid,avw.hdr.hist);

    % check the file size is 348 bytes
    fbytes = ftell(fid);
    fclose(fid);
    if ~isequal(fbytes,348),
        msg = sprintf('...file size is not 348 bytes!\n');
        ft_warning(msg);
    end

return

%----------------------------------------------------------------------------

function header_key(fid,hk)

	% Original header structures - ANALYZE 7.5
	% struct header_key                      /* header key      */
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char hkey_un0;                   /* 39 +  1         */
	%       };                               /* total=40 bytes  */

    fseek(fid,0,'bof');

    fwrite(fid, hk.sizeof_hdr(1),   'int32');  % must be 348!

    data_type = sprintf('%-10s',hk.data_type); % ensure it is 10 chars
    fwrite(fid, hk.data_type(1:10), 'uchar');

    db_name   = sprintf('%-18s',hk.db_name);   % ensure it is 18 chars
    fwrite(fid, db_name(1:18),      'uchar');

    fwrite(fid, hk.extents(1),      'int32');
    fwrite(fid, hk.session_error(1),'int16');

    regular   = sprintf('%1s',hk.regular);     % ensure it is 1 char
    fwrite(fid, regular(1),         'uchar');  % might be uint8

    %hkey_un0  = sprintf('%1s',hk.hkey_un0);    % ensure it is 1 char
    %fwrite(fid, hkey_un0(1),        'uchar');
    fwrite(fid, hk.hkey_un0(1),     'uint8');

    %    >Would you set hkey_un0 as char or uint8?
    %   Really doesn't make any difference.  As far as anyone here can remember,
    %   this was just to pad to an even byte boundary for that structure.  I guess
    %   I'd suggest setting it to a uint8 value of 0 (i.e, truly zero-valued) so
    %   that it doesn't look like anything important!
    %   Denny <hanson.dennis2@mayo.edu>

return

%----------------------------------------------------------------------------

function image_dimension(fid,dime)

	%struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
	%       char vox_units[4];               /* 16 + 4          */
	%       char cal_units[8];               /* 20 + 8          */
	%       short int unused1;               /* 28 + 2          */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int dim_un0;               /* 34 + 2          */
	%       float pixdim[8];                 /* 36 + 32         */
	%			/*
	%				pixdim[] specifies the voxel dimensions:
	%				pixdim[1] - voxel width
	%				pixdim[2] - voxel height
	%				pixdim[3] - interslice distance
	%					..etc
	%			*/
	%       float vox_offset;                /* 68 + 4          */
	%       float roi_scale;                 /* 72 + 4          */
	%       float funused1;                  /* 76 + 4          */
	%       float funused2;                  /* 80 + 4          */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       int compressed;                  /* 92 + 4          */
	%       int verified;                    /* 96 + 4          */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */

	fwrite(fid, dime.dim(1:8),      'int16');
	fwrite(fid, dime.vox_units(1:4),'uchar');
	fwrite(fid, dime.cal_units(1:8),'uchar');
	fwrite(fid, dime.unused1(1),    'int16');
	fwrite(fid, dime.datatype(1),   'int16');
	fwrite(fid, dime.bitpix(1),     'int16');
	fwrite(fid, dime.dim_un0(1),    'int16');
	fwrite(fid, dime.pixdim(1:8),   'float32');
	fwrite(fid, dime.vox_offset(1), 'float32');

    % Ensure compatibility with SPM (according to MRIcro)
    if dime.roi_scale == 0, dime.roi_scale = 0.00392157; end
	fwrite(fid, dime.roi_scale(1),  'float32');

	fwrite(fid, dime.funused1(1),   'float32');
	fwrite(fid, dime.funused2(1),   'float32');
	fwrite(fid, dime.cal_max(1),    'float32');
	fwrite(fid, dime.cal_min(1),    'float32');
	fwrite(fid, dime.compressed(1), 'int32');
	fwrite(fid, dime.verified(1),   'int32');
	fwrite(fid, dime.glmax(1),      'int32');
	fwrite(fid, dime.glmin(1),      'int32');

return

%----------------------------------------------------------------------------

function data_history(fid,hist)

	% Original header structures - ANALYZE 7.5
	%struct data_history
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       char orient;                     /* 104 + 1         */
	%       char originator[10];             /* 105 + 10        */
	%       char generated[10];              /* 115 + 10        */
	%       char scannum[10];                /* 125 + 10        */
	%       char patient_id[10];             /* 135 + 10        */
	%       char exp_date[10];               /* 145 + 10        */
	%       char exp_time[10];               /* 155 + 10        */
	%       char hist_un0[3];                /* 165 + 3         */
	%       int views                        /* 168 + 4         */
	%       int vols_added;                  /* 172 + 4         */
	%       int start_field;                 /* 176 + 4         */
	%       int field_skip;                  /* 180 + 4         */
	%       int omax;                        /* 184 + 4         */
	%       int omin;                        /* 188 + 4         */
	%       int smax;                        /* 192 + 4         */
	%       int smin;                        /* 196 + 4         */
	%       };                               /* total=200 bytes */

    descrip     = sprintf('%-80s', hist.descrip);       % 80 chars
    aux_file    = sprintf('%-24s', hist.aux_file);      % 24 chars
    originator  = sprintf('%-10s', hist.originator);    % 10 chars
    generated   = sprintf('%-10s', hist.generated);     % 10 chars
    scannum     = sprintf('%-10s', hist.scannum);       % 10 chars
    patient_id  = sprintf('%-10s', hist.patient_id);    % 10 chars
    exp_date    = sprintf('%-10s', hist.exp_date);      % 10 chars
    exp_time    = sprintf('%-10s', hist.exp_time);      % 10 chars
    hist_un0    = sprintf( '%-3s', hist.hist_un0);      %  3 chars

    % ---
    % The following should not be necessary, but I actually
    % found one instance where it was, so this totally anal
    % retentive approach became necessary, despite the
    % apparently elegant solution above to ensuring that variables
    % are the right length.

    if length(descrip) < 80,
      paddingN = 80-length(descrip);
      padding = char(repmat(double(' '),1,paddingN));
      descrip = [descrip,padding];
    end
    if length(aux_file) < 24,
      paddingN = 24-length(aux_file);
      padding = char(repmat(double(' '),1,paddingN));
      aux_file = [aux_file,padding];
    end
    if length(originator) < 10,
      paddingN = 10-length(originator);
      padding = char(repmat(double(' '),1,paddingN));
      originator = [originator, padding];
    end
    if length(generated) < 10,
      paddingN = 10-length(generated);
      padding = char(repmat(double(' '),1,paddingN));
      generated = [generated, padding];
    end
    if length(scannum) < 10,
      paddingN = 10-length(scannum);
      padding = char(repmat(double(' '),1,paddingN));
      scannum = [scannum, padding];
    end
    if length(patient_id) < 10,
      paddingN = 10-length(patient_id);
      padding = char(repmat(double(' '),1,paddingN));
      patient_id = [patient_id, padding];
    end
    if length(exp_date) < 10,
      paddingN = 10-length(exp_date);
      padding = char(repmat(double(' '),1,paddingN));
      exp_date = [exp_date, padding];
    end
    if length(exp_time) < 10,
      paddingN = 10-length(exp_time);
      padding = char(repmat(double(' '),1,paddingN));
      exp_time = [exp_time, padding];
    end
    if length(hist_un0) < 10,
      paddingN = 10-length(hist_un0);
      padding = char(repmat(double(' '),1,paddingN));
      hist_un0 = [hist_un0, padding];
    end

    % -- if you thought that was anal, try this;
    % -- lets check for unusual ASCII char values!

    if find(double(descrip)>128),
      indexStrangeChar = find(double(descrip)>128);
      descrip(indexStrangeChar) = ' ';
    end
    if find(double(aux_file)>128),
      indexStrangeChar = find(double(aux_file)>128);
      aux_file(indexStrangeChar) = ' ';
    end
    if find(double(originator)>128),
      indexStrangeChar = find(double(originator)>128);
      originator(indexStrangeChar) = ' ';
    end
    if find(double(generated)>128),
      indexStrangeChar = find(double(generated)>128);
      generated(indexStrangeChar) = ' ';
    end
    if find(double(scannum)>128),
      indexStrangeChar = find(double(scannum)>128);
      scannum(indexStrangeChar) = ' ';
    end
    if find(double(patient_id)>128),
      indexStrangeChar = find(double(patient_id)>128);
      patient_id(indexStrangeChar) = ' ';
    end
    if find(double(exp_date)>128),
      indexStrangeChar = find(double(exp_date)>128);
      exp_date(indexStrangeChar) = ' ';
    end
    if find(double(exp_time)>128),
      indexStrangeChar = find(double(exp_time)>128);
      exp_time(indexStrangeChar) = ' ';
    end
    if find(double(hist_un0)>128),
      indexStrangeChar = find(double(hist_un0)>128);
      hist_un0(indexStrangeChar) = ' ';
    end


    % --- finally, we write the fields

    fwrite(fid, descrip(1:80),    'uchar');
    fwrite(fid, aux_file(1:24),   'uchar');


    %orient      = sprintf(  '%1s', hist.orient);        %  1 char
    %fwrite(fid, orient(1),        'uchar');
    fwrite(fid, hist.orient(1),   'uint8');     % see note below on char

    fwrite(fid, originator(1:10), 'uchar');
    fwrite(fid, generated(1:10),  'uchar');
    fwrite(fid, scannum(1:10),    'uchar');
    fwrite(fid, patient_id(1:10), 'uchar');
    fwrite(fid, exp_date(1:10),   'uchar');
    fwrite(fid, exp_time(1:10),   'uchar');
    fwrite(fid, hist_un0(1:3),    'uchar');

    fwrite(fid, hist.views(1),      'int32');
    fwrite(fid, hist.vols_added(1), 'int32');
    fwrite(fid, hist.start_field(1),'int32');
    fwrite(fid, hist.field_skip(1), 'int32');
    fwrite(fid, hist.omax(1),       'int32');
    fwrite(fid, hist.omin(1),       'int32');
    fwrite(fid, hist.smax(1),       'int32');
    fwrite(fid, hist.smin(1),       'int32');

return



% Note on using char:
% The 'char orient' field in the header is intended to
% hold simply an 8-bit unsigned integer value, not the ASCII representation
% of the character for that value.  A single 'char' byte is often used to
% represent an integer value in Analyze if the known value range doesn't
% go beyond 0-255 - saves a byte over a short int, which may not mean
% much in today's computing environments, but given that this format
% has been around since the early 1980's, saving bytes here and there on
% older systems was important!  In this case, 'char' simply provides the
% byte of storage - not an indicator of the format for what is stored in
% this byte.  Generally speaking, anytime a single 'char' is used, it is
% probably meant to hold an 8-bit integer value, whereas if this has
% been dimensioned as an array, then it is intended to hold an ASCII
% character string, even if that was only a single character.
% Denny  <hanson.dennis2@mayo.edu>

% See other notes in avw_hdr_read
                                                                                                                                                                                 
        'MRP21_dH'  'MRP21_dV'  'MRP21'
        'MRP22_dH'  'MRP22_dV'  'MRP22'
        'MRP23_dH'  'MRP23_dV'  'MRP23'
        'MRP31_dH'  'MRP31_dV'  'MRP31'
        'MRP32_dH'  'MRP32_dV'  'MRP32'
        'MRP33_dH'  'MRP33_dV'  'MRP33'
        'MRP34_dH'  'MRP34_dV'  'MRP34'
        'MRP35_dH'  'MRP35_dV'  'MRP35'
        'MRP41_dH'  'MRP41_dV'  'MRP41'
        'MRP42_dH'  'MRP42_dV'  'MRP42'
        'MRP43_dH'  'MRP43_dV'  'MRP43'
        'MRP44_dH'  'MRP44_dV'  'MRP44'
        'MRP45_dH'  'MRP45_dV'  'MRP45'
        'MRP51_dH'  'MRP51_dV'  'MRP51'
        'MRP52_dH'  'MRP52_dV'  'MRP52'
        'MRP53_dH'  'MRP53_dV'  'MRP53'
        'MRP54_dH'  'MRP54_dV'  'MRP54'
        'MRP55_dH'  'MRP55_dV'  'MRP55'
        'MRP56_dH'  'MRP56_dV'  'MRP56'
        'MRP57_dH'  'MRP57_dV'  'MRP57'
        'MRT11_dH'  'MRT11_dV'  'MRT11'
        'MRT12_dH'  'MRT12_dV'  'MRT12'
        'MRT13_dH'  'MRT13_dV'  'MRT13'
        'MRT14_dH'  'MRT14_dV'  'MRT14'
        'MRT15_dH'  'MRT15_dV'  'MRT15'
        'MRT16_dH'  'MRT16_dV'  'MRT16'
        'MRT21_dH'  'MRT21_dV'  'MRT21'
        'MRT22_dH'  'MRT22_dV'  'MRT22'
        'MRT23_dH'  'MRT23_dV'  'MRT23'
        'MRT24_dH'  'MRT24_dV'  'MRT24'
        'MRT25_dH'  'MRT25_dV'  'MRT25'
        'MRT26_dH'  'MRT26_dV'  'MRT26'
        'MRT27_dH'  'MRT27_dV'  'MRT27'
        'MRT31_dH'  'MRT31_dV'  'MRT31'
        'MRT32_dH'  'MRT32_dV'  'MRT32'
        'MRT33_dH'  'MRT33_dV'  'MRT33'
        'MRT34_dH'  'MRT34_dV'  'MRT34'
        'MRT35_dH'  'MRT35_dV'  'MRT35'
        'MRT36_dH'  'MRT36_dV'  'MRT36'
        'MRT37_dH'  'MRT37_dV'  'MRT37'
        'MRT41_dH'  'MRT41_dV'  'MRT41'
        'MRT42_dH'  'MRT42_dV'  'MRT42'
        'MRT43_dH'  'MRT43_dV'  'MRT43'
        'MRT44_dH'  'MRT44_dV'  'MRT44'
        'MRT45_dH'  'MRT45_dV'  'MRT45'
        'MRT46_dH'  'MRT46_dV'  'MRT46'
        'MRT47_dH'  'MRT47_dV'  'MRT47'
        'MRT51_dH'  'MRT51_dV'  'MRT51'
        'MRT52_dH'  'MRT52_dV'  'MRT52'
        'MRT53_dH'  'MRT53_dV'  'MRT53'
        'MRT54_dH'  'MRT54_dV'  'MRT54'
        'MRT55_dH'  'MRT55_dV'  'MRT55'
        'MRT56_dH'  'MRT56_dV'  'MRT56'
        'MRT57_dH'  'MRT57_dV'  'MRT57'
        'MZC01_dH'  'MZC01_dV'  'MZC01'
        'MZC02_dH'  'MZC02_dV'  'MZC02'
        'MZC03_dH'  'MZC03_dV'  'MZC03'
        'MZC04_dH'  'MZC04_dV'  'MZC04'
        'MZF01_dH'  'MZF01_dV'  'MZF01'
        'MZF02_dH'  'MZF02_dV'  'MZF02'
        'MZF03_dH'  'MZF03_dV'  'MZF03'
        'MZO01_dH'  'MZO01_dV'  'MZO01'
        'MZO02_dH'  'MZO02_dV'  'MZO02'
        'MZO03_dH'  'MZO03_dV'  'MZO03'
        'MZP01_dH'  'MZP01_dV'  'MZP01'
        };
      ctf275_planar_combined = label(:,3);
      label = label(:,1:2);
      
    case {'neuromag122' 'neuromag122alt'}
      % this is the combination of the two versions (with and without space)
      label = {
        'MEG 001'  'MEG 002'  'MEG 001+002'
        'MEG 003'  'MEG 004'  'MEG 003+004'
        'MEG 005'  'MEG 006'  'MEG 005+006'
        'MEG 007'  'MEG 008'  'MEG 007+008'
        'MEG 009'  'MEG 010'  'MEG 009+010'
        'MEG 011'  'MEG 012'  'MEG 011+012'
        'MEG 013'  'MEG 014'  'MEG 013+014'
        'MEG 015'  'MEG 016'  'MEG 015+016'
        'MEG 017'  'MEG 018'  'MEG 017+018'
        'MEG 019'  'MEG 020'  'MEG 019+020'
        'MEG 021'  'MEG 022'  'MEG 021+022'
        'MEG 023'  'MEG 024'  'MEG 023+024'
        'MEG 025'  'MEG 026'  'MEG 025+026'
        'MEG 027'  'MEG 028'  'MEG 027+028'
        'MEG 029'  'MEG 030'  'MEG 029+030'
        'MEG 031'  'MEG 032'  'MEG 031+032'
        'MEG 033'  'MEG 034'  'MEG 033+034'
        'MEG 035'  'MEG 036'  'MEG 035+036'
        'MEG 037'  'MEG 038'  'MEG 037+038'
        'MEG 039'  'MEG 040'  'MEG 039+040'
        'MEG 041'  'MEG 042'  'MEG 041+042'
        'MEG 043'  'MEG 044'  'MEG 043+044'
        'MEG 045'  'MEG 046'  'MEG 045+046'
        'MEG 047'  'MEG 048'  'MEG 047+048'
        'MEG 049'  'MEG 050'  'MEG 049+050'
        'MEG 051'  'MEG 052'  'MEG 051+052'
        'MEG 053'  'MEG 054'  'MEG 053+054'
        'MEG 055'  'MEG 056'  'MEG 055+056'
        'MEG 057'  'MEG 058'  'MEG 057+058'
        'MEG 059'  'MEG 060'  'MEG 059+060'
        'MEG 061'  'MEG 062'  'MEG 061+062'
        'MEG 063'  'MEG 064'  'MEG 063+064'
        'MEG 065'  'MEG 066'  'MEG 065+066'
        'MEG 067'  'MEG 068'  'MEG 067+068'
        'MEG 069'  'MEG 070'  'MEG 069+070'
        'MEG 071'  'MEG 072'  'MEG 071+072'
        'MEG 073'  'MEG 074'  'MEG 073+074'
        'MEG 075'  'MEG 076'  'MEG 075+076'
        'MEG 077'  'MEG 078'  'MEG 077+078'
        'MEG 079'  'MEG 080'  'MEG 079+080'
        'MEG 081'  'MEG 082'  'MEG 081+082'
        'MEG 083'  'MEG 084'  'MEG 083+084'
        'MEG 085'  'MEG 086'  'MEG 085+086'
        'MEG 087'  'MEG 088'  'MEG 087+088'
        'MEG 089'  'MEG 090'  'MEG 089+090'
        'MEG 091'  'MEG 092'  'MEG 091+092'
        'MEG 093'  'MEG 094'  'MEG 093+094'
        'MEG 095'  'MEG 096'  'MEG 095+096'
        'MEG 097'  'MEG 098'  'MEG 097+098'
        'MEG 099'  'MEG 100'  'MEG 099+100'
        'MEG 101'  'MEG 102'  'MEG 101+102'
        'MEG 103'  'MEG 104'  'MEG 103+104'
        'MEG 105'  'MEG 106'  'MEG 105+106'
        'MEG 107'  'MEG 108'  'MEG 107+108'
        'MEG 109'  'MEG 110'  'MEG 109+110'
        'MEG 111'  'MEG 112'  'MEG 111+112'
        'MEG 113'  'MEG 114'  'MEG 113+114'
        'MEG 115'  'MEG 116'  'MEG 115+116'
        'MEG 117'  'MEG 118'  'MEG 117+118'
        'MEG 119'  'MEG 120'  'MEG 119+120'
        'MEG 121'  'MEG 122'  'MEG 121+122'
        % this is an alternative set of labels without a space in them
        'MEG001'  'MEG002'  'MEG001+002'
        'MEG003'  'MEG004'  'MEG003+004'
        'MEG005'  'MEG006'  'MEG005+006'
        'MEG007'  'MEG008'  'MEG007+008'
        'MEG009'  'MEG010'  'MEG009+010'
        'MEG011'  'MEG012'  'MEG011+012'
        'MEG013'  'MEG014'  'MEG013+014'
        'MEG015'  'MEG016'  'MEG015+016'
        'MEG017'  'MEG018'  'MEG017+018'
        'MEG019'  'MEG020'  'MEG019+020'
        'MEG021'  'MEG022'  'MEG021+022'
        'MEG023'  'MEG024'  'MEG023+024'
        'MEG025'  'MEG026'  'MEG025+026'
        'MEG027'  'MEG028'  'MEG027+028'
        'MEG029'  'MEG030'  'MEG029+030'
        'MEG031'  'MEG032'  'MEG031+032'
        'MEG033'  'MEG034'  'MEG033+034'
        'MEG035'  'MEG036'  'MEG035+036'
        'MEG037'  'MEG038'  'MEG037+038'
        'MEG039'  'MEG040'  'MEG039+040'
        'MEG041'  'MEG042'  'MEG041+042'
        'MEG043'  'MEG044'  'MEG043+044'
        'MEG045'  'MEG046'  'MEG045+046'
        'MEG047'  'MEG048'  'MEG047+048'
        'MEG049'  'MEG050'  'MEG049+050'
        'MEG051'  'MEG052'  'MEG051+052'
        'MEG053'  'MEG054'  'MEG053+054'
        'MEG055'  'MEG056'  'MEG055+056'
        'MEG057'  'MEG058'  'MEG057+058'
        'MEG059'  'MEG060'  'MEG059+060'
        'MEG061'  'MEG062'  'MEG061+062'
        'MEG063'  'MEG064'  'MEG063+064'
        'MEG065'  'MEG066'  'MEG065+066'
        'MEG067'  'MEG068'  'MEG067+068'
        'MEG069'  'MEG070'  'MEG069+070'
        'MEG071'  'MEG072'  'MEG071+072'
        'MEG073'  'MEG074'  'MEG073+074'
        'MEG075'  'MEG076'  'MEG075+076'
        'MEG077'  'MEG078'  'MEG077+078'
        'MEG079'  'MEG080'  'MEG079+080'
        'MEG081'  'MEG082'  'MEG081+082'
        'MEG083'  'MEG084'  'MEG083+084'
        'MEG085'  'MEG086'  'MEG085+086'
        'MEG087'  'MEG088'  'MEG087+088'
        'MEG089'  'MEG090'  'MEG089+090'
        'MEG091'  'MEG092'  'MEG091+092'
        'MEG093'  'MEG094'  'MEG093+094'
        'MEG095'  'MEG096'  'MEG095+096'
        'MEG097'  'MEG098'  'MEG097+098'
        'MEG099'  'MEG100'  'MEG099+100'
        'MEG101'  'MEG102'  'MEG101+102'
        'MEG103'  'MEG104'  'MEG103+104'
        'MEG105'  'MEG106'  'MEG105+106'
        'MEG107'  'MEG108'  'MEG107+108'
        'MEG109'  'MEG110'  'MEG109+110'
        'MEG111'  'MEG112'  'MEG111+112'
        'MEG113'  'MEG114'  'MEG113+114'
        'MEG115'  'MEG116'  'MEG115+116'
        'MEG117'  'MEG118'  'MEG117+118'
        'MEG119'  'MEG120'  'MEG119+120'
        'MEG121'  'MEG122'  'MEG121+122'
        };
      neuromag122_combined = label(:,3);
      neuromag122alt_combined = label(:,3);
      label = label(:,1:2);
      
    case {'neuromag306' 'neuromag306alt'}
      % this is the combination of the two versions (with and without space)
      label = {
        'MEG 0112'  'MEG 0113'  'MEG 0111'  'MEG 0112+0113'
        'MEG 0122'  'MEG 0123'  'MEG 0121'  'MEG 0122+0123'
        'MEG 0132'  'MEG 0133'  'MEG 0131'  'MEG 0132+0133'
        'MEG 0142'  'MEG 0143'  'MEG 0141'  'MEG 0142+0143'
        'MEG 0212'  'MEG 0213'  'MEG 0211'  'MEG 0212+0213'
        'MEG 0222'  'MEG 0223'  'MEG 0221'  'MEG 0222+0223'
        'MEG 0232'  'MEG 0233'  'MEG 0231'  'MEG 0232+0233'
        'MEG 0242'  'MEG 0243'  'MEG 0241'  'MEG 0242+0243'
        'MEG 0312'  'MEG 0313'  'MEG 0311'  'MEG 0312+0313'
        'MEG 0322'  'MEG 0323'  'MEG 0321'  'MEG 0322+0323'
        'MEG 0332'  'MEG 0333'  'MEG 0331'  'MEG 0332+0333'
        'MEG 0342'  'MEG 0343'  'MEG 0341'  'MEG 0342+0343'
        'MEG 0412'  'MEG 0413'  'MEG 0411'  'MEG 0412+0413'
        'MEG 0422'  'MEG 0423'  'MEG 0421'  'MEG 0422+0423'
        'MEG 0432'  'MEG 0433'  'MEG 0431'  'MEG 0432+0433'
        'MEG 0442'  'MEG 0443'  'MEG 0441'  'MEG 0442+0443'
        'MEG 0512'  'MEG 0513'  'MEG 0511'  'MEG 0512+0513'
        'MEG 0522'  'MEG 0523'  'MEG 0521'  'MEG 0522+0523'
        'MEG 0532'  'MEG 0533'  'MEG 0531'  'MEG 0532+0533'
        'MEG 0542'  'MEG 0543'  'MEG 0541'  'MEG 0542+0543'
        'MEG 0612'  'MEG 0613'  'MEG 0611'  'MEG 0612+0613'
        'MEG 0622'  'MEG 0623'  'MEG 0621'  'MEG 0622+0623'
        'MEG 0632'  'MEG 0633'  'MEG 0631'  'MEG 0632+0633'
        'MEG 0642'  'MEG 0643'  'MEG 0641'  'MEG 0642+0643'
        'MEG 0712'  'MEG 0713'  'MEG 0711'  'MEG 0712+0713'
        'MEG 0722'  'MEG 0723'  'MEG 0721'  'MEG 0722+0723'
        'MEG 0732'  'MEG 0733'  'MEG 0731'  'MEG 0732+0733'
        'MEG 0742'  'MEG 0743'  'MEG 0741'  'MEG 0742+0743'
        'MEG 0812'  'MEG 0813'  'MEG 0811'  'MEG 0812+0813'
        'MEG 0822'  'MEG 0823'  'MEG 0821'  'MEG 0822+0823'
        'MEG 0912'  'MEG 0913'  'MEG 0911'  'MEG 0912+0913'
        'MEG 0922'  'MEG 0923'  'MEG 0921'  'MEG 0922+0923'
        'MEG 0932'  'MEG 0933'  'MEG 0931'  'MEG 0932+0933'
        'MEG 0942'  'MEG 0943'  'MEG 0941'  'MEG 0942+0943'
        'MEG 1012'  'MEG 1013'  'MEG 1011'  'MEG 1012+1013'
        'MEG 1022'  'MEG 1023'  'MEG 1021'  'MEG 1022+1023'
        'MEG 1032'  'MEG 1033'  'MEG 1031'  'MEG 1032+1033'
        'MEG 1042'  'MEG 1043'  'MEG 1041'  'MEG 1042+1043'
        'MEG 1112'  'MEG 1113'  'MEG 1111'  'MEG 1112+1113'
        'MEG 1122'  'MEG 1123'  'MEG 1121'  'MEG 1122+1123'
        'MEG 1132'  'MEG 1133'  'MEG 1131'  'MEG 1132+1133'
        'MEG 1142'  'MEG 1143'  'MEG 1141'  'MEG 1142+1143'
        'MEG 1212'  'MEG 1213'  'MEG 1211'  'MEG 1212+1213'
        'MEG 1222'  'MEG 1223'  'MEG 1221'  'MEG 1222+1223'
        'MEG 1232'  'MEG 1233'  'MEG 1231'  'MEG 1232+1233'
        'MEG 1242'  'MEG 1243'  'MEG 1241'  'MEG 1242+1243'
        'MEG 1312'  'MEG 1313'  'MEG 1311'  'MEG 1312+1313'
        'MEG 1322'  'MEG 1323'  'MEG 1321'  'MEG 1322+1323'
        'MEG 1332'  'MEG 1333'  'MEG 1331'  'MEG 1332+1333'
        'MEG 1342'  'MEG 1343'  'MEG 1341'  'MEG 1342+1343'
        'MEG 1412'  'MEG 1413'  'MEG 1411'  'MEG 1412+1413'
        'MEG 1422'  'MEG 1423'  'MEG 1421'  'MEG 1422+1423'
        'MEG 1432'  'MEG 1433'  'MEG 1431'  'MEG 1432+1433'
        'MEG 1442'  'MEG 1443'  'MEG 1441'  'MEG 1442+1443'
        'MEG 1512'  'MEG 1513'  'MEG 1511'  'MEG 1512+1513'
        'MEG 1522'  'MEG 1523'  'MEG 1521'  'MEG 1522+1523'
        'MEG 1532'  'MEG 1533'  'MEG 1531'  'MEG 1532+1533'
        'MEG 1542'  'MEG 1543'  'MEG 1541'  'MEG 1542+1543'
        'MEG 1612'  'MEG 1613'  'MEG 1611'  'MEG 1612+1613'
        'MEG 1622'  'MEG 1623'  'MEG 1621'  'MEG 1622+1623'
        'MEG 1632'  'MEG 1633'  'MEG 1631'  'MEG 1632+1633'
        'MEG 1642'  'MEG 1643'  'MEG 1641'  'MEG 1642+1643'
        'MEG 1712'  'MEG 1713'  'MEG 1711'  'MEG 1712+1713'
        'MEG 1722'  'MEG 1723'  'MEG 1721'  'MEG 1722+1723'
        'MEG 1732'  'MEG 1733'  'MEG 1731'  'MEG 1732+1733'
        'MEG 1742'  'MEG 1743'  'MEG 1741'  'MEG 1742+1743'
        'MEG 1812'  'MEG 1813'  'MEG 1811'  'MEG 1812+1813'
        'MEG 1822'  'MEG 1823'  'MEG 1821'  'MEG 1822+1823'
        'MEG 1832'  'MEG 1833'  'MEG 1831'  'MEG 1832+1833'
        'MEG 1842'  'MEG 1843'  'MEG 1841'  'MEG 1842+1843'
        'MEG 1912'  'MEG 1913'  'MEG 1911'  'MEG 1912+1913'
        'MEG 1922'  'MEG 1923'  'MEG 1921'  'MEG 1922+1923'
        'MEG 1932'  'MEG 1933'  'MEG 1931'  'MEG 1932+1933'
        'MEG 1942'  'MEG 1943'  'MEG 1941'  'MEG 1942+1943'
        'MEG 2012'  'MEG 2013'  'MEG 2011'  'MEG 2012+2013'
        'MEG 2022'  'MEG 2023'  'MEG 2021'  'MEG 2022+2023'
        'MEG 2032'  'MEG 2033'  'MEG 2031'  'MEG 2032+2033'
        'MEG 2042'  'MEG 2043'  'MEG 2041'  'MEG 2042+2043'
        'MEG 2112'  'MEG 2113'  'MEG 2111'  'MEG 2112+2113'
        'MEG 2122'  'MEG 2123'  'MEG 2121'  'MEG 2122+2123'
        'MEG 2132'  'MEG 2133'  'MEG 2131'  'MEG 2132+2133'
        'MEG 2142'  'MEG 2143'  'MEG 2141'  'MEG 2142+2143'
        'MEG 2212'  'MEG 2213'  'MEG 2211'  'MEG 2212+2213'
        'MEG 2222'  'MEG 2223'  'MEG 2221'  'MEG 2222+2223'
        'MEG 2232'  'MEG 2233'  'MEG 2231'  'MEG 2232+2233'
        'MEG 2242'  'MEG 2243'  'MEG 2241'  'MEG 2242+2243'
        'MEG 2312'  'MEG 2313'  'MEG 2311'  'MEG 2312+2313'
        'MEG 2322'  'MEG 2323'  'MEG 2321'  'MEG 2322+2323'
        'MEG 2332'  'MEG 2333'  'MEG 2331'  'MEG 2332+2333'
        'MEG 2342'  'MEG 2343'  'MEG 2341'  'MEG 2342+2343'
        'MEG 2412'  'MEG 2413'  'MEG 2411'  'MEG 2412+2413'
        'MEG 2422'  'MEG 2423'  'MEG 2421'  'MEG 2422+2423'
        'MEG 2432'  'MEG 2433'  'MEG 2431'  'MEG 2432+2433'
        'MEG 2442'  'MEG 2443'  'MEG 2441'  'MEG 2442+2443'
        'MEG 2512'  'MEG 2513'  'MEG 2511'  'MEG 2512+2513'
        'MEG 2522'  'MEG 2523'  'MEG 2521'  'MEG 2522+2523'
        'MEG 2532'  'MEG 2533'  'MEG 2531'  'MEG 2532+2533'
        'MEG 2542'  'MEG 2543'  'MEG 2541'  'MEG 2542+2543'
        'MEG 2612'  'MEG 2613'  'MEG 2611'  'MEG 2612+2613'
        'MEG 2622'  'MEG 2623'  'MEG 2621'  'MEG 2622+2623'
        'MEG 2632'  'MEG 2633'  'MEG 2631'  'MEG 2632+2633'
        'MEG 2642'  'MEG 2643'  'MEG 2641'  'MEG 2642+2643'
        % this is an alternative set of labels without a space in them
        'MEG0112'  'MEG0113'  'MEG0111'  'MEG0112+0113'
        'MEG0122'  'MEG0123'  'MEG0121'  'MEG0122+0123'
        'MEG0132'  'MEG0133'  'MEG0131'  'MEG0132+0133'
        'MEG0142'  'MEG0143'  'MEG0141'  'MEG0142+0143'
        'MEG0212'  'MEG0213'  'MEG0211'  'MEG0212+0213'
        'MEG0222'  'MEG0223'  'MEG0221'  'MEG0222+0223'
        'MEG0232'  'MEG0233'  'MEG0231'  'MEG0232+0233'
        'MEG0242'  'MEG0243'  'MEG0241'  'MEG0242+0243'
        'MEG0312'  'MEG0313'  'MEG0311'  'MEG0312+0313'
        'MEG0322'  'MEG0323'  'MEG0321'  'MEG0322+0323'
        'MEG0332'  'MEG0333'  'MEG0331'  'MEG0332+0333'
        'MEG0342'  'MEG0343'  'MEG0341'  'MEG0342+0343'
        'MEG0412'  'MEG0413'  'MEG0411'  'MEG0412+0413'
        'MEG0422'  'MEG0423'  'MEG0421'  'MEG0422+0423'
        'MEG0432'  'MEG0433'  'MEG0431'  'MEG0432+0433'
        'MEG0442'  'MEG0443'  'MEG0441'  'MEG0442+0443'
        'MEG0512'  'MEG0513'  'MEG0511'  'MEG0512+0513'
        'MEG0522'  'MEG0523'  'MEG0521'  'MEG0522+0523'
        'MEG0532'  'MEG0533'  'MEG0531'  'MEG0532+0533'
        'MEG0542'  'MEG0543'  'MEG0541'  'MEG0542+0543'
        'MEG0612'  'MEG0613'  'MEG0611'  'MEG0612+0613'
        'MEG0622'  'MEG0623'  'MEG0621'  'MEG0622+0623'
        'MEG0632'  'MEG0633'  'MEG0631'  'MEG0632+0633'
        'MEG0642'  'MEG0643'  'MEG0641'  'MEG0642+0643'
        'MEG0712'  'MEG0713'  'MEG0711'  'MEG0712+0713'
        'MEG0722'  'MEG0723'  'MEG0721'  'MEG0722+0723'
        'MEG0732'  'MEG0733'  'MEG0731'  'MEG0732+0733'
        'MEG0742'  'MEG0743'  'MEG0741'  'MEG0742+0743'
        'MEG0812'  'MEG0813'  'MEG0811'  'MEG0812+0813'
        'MEG0822'  'MEG0823'  'MEG0821'  'MEG0822+0823'
        'MEG0912'  'MEG0913'  'MEG0911'  'MEG0912+0913'
        'MEG0922'  'MEG0923'  'MEG0921'  'MEG0922+0923'
        'MEG0932'  'MEG0933'  'MEG0931'  'MEG0932+0933'
        'MEG0942'  'MEG0943'  'MEG0941'  'MEG0942+0943'
        'MEG1012'  'MEG1013'  'MEG1011'  'MEG1012+1013'
        'MEG1022'  'MEG1023'  'MEG1021'  'MEG1022+1023'
        'MEG1032'  'MEG1033'  'MEG1031'  'MEG1032+1033'
        'MEG1042'  'MEG1043'  'MEG1041'  'MEG1042+1043'
        'MEG1112'  'MEG1113'  'MEG1111'  'MEG1112+1113'
        'MEG1122'  'MEG1123'  'MEG1121'  'MEG1122+1123'
        'MEG1132'  'MEG1133'  'MEG1131'  'MEG1132+1133'
        'MEG1142'  'MEG1143'  'MEG1141'  'MEG1142+1143'
        'MEG1212'  'MEG1213'  'MEG1211'  'MEG1212+1213'
        'MEG1222'  'MEG1223'  'MEG1221'  'MEG1222+1223'
        'MEG1232'  'MEG1233'  'MEG1231'  'MEG1232+1233'
        'MEG1242'  'MEG1243'  'MEG1241'  'MEG1242+1243'
        'MEG1312'  'MEG1313'  'MEG1311'  'MEG1312+1313'
        'MEG1322'  'MEG1323'  'MEG1321'  'MEG1322+1323'
        'MEG1332'  'MEG1333'  'MEG1331'  'MEG1332+1333'
        'MEG1342'  'MEG1343'  'MEG1341'  'MEG1342+1343'
        'MEG1412'  'MEG1413'  'MEG1411'  'MEG1412+1413'
        'MEG1422'  'MEG1423'  'MEG1421'  'MEG1422+1423'
        'MEG1432'  'MEG1433'  'MEG1431'  'MEG1432+1433'
        'MEG1442'  'MEG1443'  'MEG1441'  'MEG1442+1443'
        'MEG1512'  'MEG1513'  'MEG1511'  'MEG1512+1513'
        'MEG1522'  'MEG1523'  'MEG1521'  'MEG1522+1523'
        'MEG1532'  'MEG1533'  'MEG1531'  'MEG1532+1533'
        'MEG1542'  'MEG1543'  'MEG1541'  'MEG1542+1543'
        'MEG1612'  'MEG1613'  'MEG1611'  'MEG1612+1613'
        'MEG1622'  'MEG1623'  'MEG1621'  'MEG1622+1623'
        'MEG1632'  'MEG1633'  'MEG1631'  'MEG1632+1633'
        'MEG1642'  'MEG1643'  'MEG1641'  'MEG1642+1643'
        'MEG1712'  'MEG1713'  'MEG1711'  'MEG1712+1713'
        'MEG1722'  'MEG1723'  'MEG1721'  'MEG1722+1723'
        'MEG1732'  'MEG1733'  'MEG1731'  'MEG1732+1733'
        'MEG1742'  'MEG1743'  'MEG1741'  'MEG1742+1743'
        'MEG1812'  'MEG1813'  'MEG1811'  'MEG1812+1813'
        'MEG1822'  'MEG1823'  'MEG1821'  'MEG1822+1823'
        'MEG1832'  'MEG1833'  'MEG1831'  'MEG1832+1833'
        'MEG1842'  'MEG1843'  'MEG1841'  'MEG1842+1843'
        'MEG1912'  'MEG1913'  'MEG1911'  'MEG1912+1913'
        'MEG1922'  'MEG1923'  'MEG1921'  'MEG1922+1923'
        'MEG1932'  'MEG1933'  'MEG1931'  'MEG1932+1933'
        'MEG1942'  'MEG1943'  'MEG1941'  'MEG1942+1943'
        'MEG2012'  'MEG2013'  'MEG2011'  'MEG2012+2013'
        'MEG2022'  'MEG2023'  'MEG2021'  'MEG2022+2023'
        'MEG2032'  'MEG2033'  'MEG2031'  'MEG2032+2033'
        'MEG2042'  'MEG2043'  'MEG2041'  'MEG2042+2043'
        'MEG2112'  'MEG2113'  'MEG2111'  'MEG2112+2113'
        'MEG2122'  'MEG2123'  'MEG2121'  'MEG2122+2123'
        'MEG2132'  'MEG2133'  'MEG2131'  'MEG2132+2133'
        'MEG2142'  'MEG2143'  'MEG2141'  'MEG2142+2143'
        'MEG2212'  'MEG2213'  'MEG2211'  'MEG2212+2213'
        'MEG2222'  'MEG2223'  'MEG2221'  'MEG2222+2223'
        'MEG2232'  'MEG2233'  'MEG2231'  'MEG2232+2233'
        'MEG2242'  'MEG2243'  'MEG2241'  'MEG2242+2243'
        'MEG2312'  'MEG2313'  'MEG2311'  'MEG2312+2313'
        'MEG2322'  'MEG2323'  'MEG2321'  'MEG2322+2323'
        'MEG2332'  'MEG2333'  'MEG2331'  'MEG2332+2333'
        'MEG2342'  'MEG2343'  'MEG2341'  'MEG2342+2343'
        'MEG2412'  'MEG2413'  'MEG2411'  'MEG2412+2413'
        'MEG2422'  'MEG2423'  'MEG2421'  'MEG2422+2423'
        'MEG2432'  'MEG2433'  'MEG2431'  'MEG2432+2433'
        'MEG2442'  'MEG2443'  'MEG2441'  'MEG2442+2443'
        'MEG2512'  'MEG2513'  'MEG2511'  'MEG2512+2513'
        'MEG2522'  'MEG2523'  'MEG2521'  'MEG2522+2523'
        'MEG2532'  'MEG2533'  'MEG2531'  'MEG2532+2533'
        'MEG2542'  'MEG2543'  'MEG2541'  'MEG2542+2543'
        'MEG2612'  'MEG2613'  'MEG2611'  'MEG2612+2613'
        'MEG2622'  'MEG2623'  'MEG2621'  'MEG2622+2623'
        'MEG2632'  'MEG2633'  'MEG2631'  'MEG2632+2633'
        'MEG2642'  'MEG2643'  'MEG2641'  'MEG2642+2643'
        };
      neuromag306_mag      = label(:,3);
      neuromag306_planar   = label(:,[1 2]);
      neuromag306_combined = label(:,[3 4]); % magnetometers and combined channels
      label                = label(:,1:3);
      
    case 'eeg1020'
      label = {
        'Fp1'
        'Fpz'
        'Fp2'
        'F7'
        'F3'
        'Fz'
        'F4'
        'F8'
        'T7'
        'C3'
        'Cz'
        'C4'
        'T8'
        'P7'
        'P3'
        'Pz'
        'P4'
        'P8'
        'O1'
        'Oz'
        'O2'};
      
      % Add also reference and some alternative labels that might be used
      label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');
      
    case 'eeg1010'
      label = {
        'Fp1'
        'Fpz'
        'Fp2'
        'AF9'
        'AF7'
        'AF5'
        'AF3'
        'AF1'
        'AFz'
        'AF2'
        'AF4'
        'AF6'
        'AF8'
        'AF10'
        'F9'
        'F7'
        'F5'
        'F3'
        'F1'
        'Fz'
        'F2'
        'F4'
        'F6'
        'F8'
        'F10'
        'FT9'
        'FT7'
        'FC5'
        'FC3'
        'FC1'
        'FCz'
        'FC2'
        'FC4'
        'FC6'
        'FT8'
        'FT10'
        'T9'
        'T7'
        'C5'
        'C3'
        'C1'
        'Cz'
        'C2'
        'C4'
        'C6'
        'T8'
        'T10'
        'TP9'
        'TP7'
        'CP5'
        'CP3'
        'CP1'
        'CPz'
        'CP2'
        'CP4'
        'CP6'
        'TP8'
        'TP10'
        'P9'
        'P7'
        'P5'
        'P3'
        'P1'
        'Pz'
        'P2'
        'P4'
        'P6'
        'P8'
        'P10'
        'PO9'
        'PO7'
        'PO5'
        'PO3'
        'PO1'
        'POz'
        'PO2'
        'PO4'
        'PO6'
        'PO8'
        'PO10'
        'O1'
        'Oz'
        'O2'
        'I1'
        'Iz'
        'I2'
        };

      % Add also reference and some alternative labels that might be used
      label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

    case 'eeg1005'
      label = {
        'Fp1'
        'Fpz'
        'Fp2'
        'AF9'
        'AF7'
        'AF5'
        'AF3'
        'AF1'
        'AFz'
        'AF2'
        'AF4'
        'AF6'
        'AF8'
        'AF10'
        'F9'
        'F7'
        'F5'
        'F3'
        'F1'
        'Fz'
        'F2'
        'F4'
        'F6'
        'F8'
        'F10'
        'FT9'
        'FT7'
        'FC5'
        'FC3'
        'FC1'
        'FCz'
        'FC2'
        'FC4'
        'FC6'
        'FT8'
        'FT10'
        'T9'
        'T7'
        'C5'
        'C3'
        'C1'
        'Cz'
        'C2'
        'C4'
        'C6'
        'T8'
        'T10'
        'TP9'
        'TP7'
        'CP5'
        'CP3'
        'CP1'
        'CPz'
        'CP2'
        'CP4'
        'CP6'
        'TP8'
        'TP10'
        'P9'
        'P7'
        'P5'
        'P3'
        'P1'
        'Pz'
        'P2'
        'P4'
        'P6'
        'P8'
        'P10'
        'PO9'
        'PO7'
        'PO5'
        'PO3'
        'PO1'
        'POz'
        'PO2'
        'PO4'
        'PO6'
        'PO8'
        'PO10'
        'O1'
        'Oz'
        'O2'
        'I1'
        'Iz'
        'I2'
        'AFp9h'
        'AFp7h'
        'AFp5h'
        'AFp3h'
        'AFp1h'
        'AFp2h'
        'AFp4h'
        'AFp6h'
        'AFp8h'
        'AFp10h'
        'AFF9h'
        'AFF7h'
        'AFF5h'
        'AFF3h'
        'AFF1h'
        'AFF2h'
        'AFF4h'
        'AFF6h'
        'AFF8h'
        'AFF10h'
        'FFT9h'
        'FFT7h'
        'FFC5h'
        'FFC3h'
        'FFC1h'
        'FFC2h'
        'FFC4h'
        'FFC6h'
        'FFT8h'
        'FFT10h'
        'FTT9h'
        'FTT7h'
        'FCC5h'
        'FCC3h'
        'FCC1h'
        'FCC2h'
        'FCC4h'
        'FCC6h'
        'FTT8h'
        'FTT10h'
        'TTP9h'
        'TTP7h'
        'CCP5h'
        'CCP3h'
        'CCP1h'
        'CCP2h'
        'CCP4h'
        'CCP6h'
        'TTP8h'
        'TTP10h'
        'TPP9h'
        'TPP7h'
        'CPP5h'
        'CPP3h'
        'CPP1h'
        'CPP2h'
        'CPP4h'
        'CPP6h'
        'TPP8h'
        'TPP10h'
        'PPO9h'
        'PPO7h'
        'PPO5h'
        'PPO3h'
        'PPO1h'
        'PPO2h'
        'PPO4h'
        'PPO6h'
        'PPO8h'
        'PPO10h'
        'POO9h'
        'POO7h'
        'POO5h'
        'POO3h'
        'POO1h'
        'POO2h'
        'POO4h'
        'POO6h'
        'POO8h'
        'POO10h'
        'OI1h'
        'OI2h'
        'Fp1h'
        'Fp2h'
        'AF9h'
        'AF7h'
        'AF5h'
        'AF3h'
        'AF1h'
        'AF2h'
        'AF4h'
        'AF6h'
        'AF8h'
        'AF10h'
        'F9h'
        'F7h'
        'F5h'
        'F3h'
        'F1h'
        'F2h'
        'F4h'
        'F6h'
        'F8h'
        'F10h'
        'FT9h'
        'FT7h'
        'FC5h'
        'FC3h'
        'FC1h'
        'FC2h'
        'FC4h'
        'FC6h'
        'FT8h'
        'FT10h'
        'T9h'
        'T7h'
        'C5h'
        'C3h'
        'C1h'
        'C2h'
        'C4h'
        'C6h'
        'T8h'
        'T10h'
        'TP9h'
        'TP7h'
        'CP5h'
        'CP3h'
        'CP1h'
        'CP2h'
        'CP4h'
        'CP6h'
        'TP8h'
        'TP10h'
        'P9h'
        'P7h'
        'P5h'
        'P3h'
        'P1h'
        'P2h'
        'P4h'
        'P6h'
        'P8h'
        'P10h'
        'PO9h'
        'PO7h'
        'PO5h'
        'PO3h'
        'PO1h'
        'PO2h'
        'PO4h'
        'PO6h'
        'PO8h'
        'PO10h'
        'O1h'
        'O2h'
        'I1h'
        'I2h'
        'AFp9'
        'AFp7'
        'AFp5'
        'AFp3'
        'AFp1'
        'AFpz'
        'AFp2'
        'AFp4'
        'AFp6'
        'AFp8'
        'AFp10'
        'AFF9'
        'AFF7'
        'AFF5'
        'AFF3'
        'AFF1'
        'AFFz'
        'AFF2'
        'AFF4'
        'AFF6'
        'AFF8'
        'AFF10'
        'FFT9'
        'FFT7'
        'FFC5'
        'FFC3'
        'FFC1'
        'FFCz'
        'FFC2'
        'FFC4'
        'FFC6'
        'FFT8'
        'FFT10'
        'FTT9'
        'FTT7'
        'FCC5'
        'FCC3'
        'FCC1'
        'FCCz'
        'FCC2'
        'FCC4'
        'FCC6'
        'FTT8'
        'FTT10'
        'TTP9'
        'TTP7'
        'CCP5'
        'CCP3'
        'CCP1'
        'CCPz'
        'CCP2'
        'CCP4'
        'CCP6'
        'TTP8'
        'TTP10'
        'TPP9'
        'TPP7'
        'CPP5'
        'CPP3'
        'CPP1'
        'CPPz'
        'CPP2'
        'CPP4'
        'CPP6'
        'TPP8'
        'TPP10'
        'PPO9'
        'PPO7'
        'PPO5'
        'PPO3'
        'PPO1'
        'PPOz'
        'PPO2'
        'PPO4'
        'PPO6'
        'PPO8'
        'PPO10'
        'POO9'
        'POO7'
        'POO5'
        'POO3'
        'POO1'
        'POOz'
        'POO2'
        'POO4'
        'POO6'
        'POO8'
        'POO10'
        'OI1'
        'OIz'
        'OI2'
        };

      % Add also reference and some alternative labels that might be used
      label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');
      
    case 'ext1020'
      % start with the eeg1005 list
      label = {
        'Fp1'
        'Fpz'
        'Fp2'
        'AF9'
        'AF7'
        'AF5'
        'AF3'
        'AF1'
        'AFz'
        'AF2'
        'AF4'
        'AF6'
        'AF8'
        'AF10'
        'F9'
        'F7'
        'F5'
        'F3'
        'F1'
        'Fz'
        'F2'
        'F4'
        'F6'
        'F8'
        'F10'
        'FT9'
        'FT7'
        'FC5'
        'FC3'
        'FC1'
        'FCz'
        'FC2'
        'FC4'
        'FC6'
        'FT8'
        'FT10'
        'T9'
        'T7'
        'C5'
        'C3'
        'C1'
        'Cz'
        'C2'
        'C4'
        'C6'
        'T8'
        'T10'
        'TP9'
        'TP7'
        'CP5'
        'CP3'
        'CP1'
        'CPz'
        'CP2'
        'CP4'
        'CP6'
        'TP8'
        'TP10'
        'P9'
        'P7'
        'P5'
        'P3'
        'P1'
        'Pz'
        'P2'
        'P4'
        'P6'
        'P8'
        'P10'
        'PO9'
        'PO7'
        'PO5'
        'PO3'
        'PO1'
        'POz'
        'PO2'
        'PO4'
        'PO6'
        'PO8'
        'PO10'
        'O1'
        'Oz'
        'O2'
        'I1'
        'Iz'
        'I2'
        'AFp9h'
        'AFp7h'
        'AFp5h'
        'AFp3h'
        'AFp1h'
        'AFp2h'
        'AFp4h'
        'AFp6h'
        'AFp8h'
        'AFp10h'
        'AFF9h'
        'AFF7h'
        'AFF5h'
        'AFF3h'
        'AFF1h'
        'AFF2h'
        'AFF4h'
        'AFF6h'
        'AFF8h'
        'AFF10h'
        'FFT9h'
        'FFT7h'
        'FFC5h'
        'FFC3h'
        'FFC1h'
        'FFC2h'
        'FFC4h'
        'FFC6h'
        'FFT8h'
        'FFT10h'
        'FTT9h'
        'FTT7h'
        'FCC5h'
        'FCC3h'
        'FCC1h'
        'FCC2h'
        'FCC4h'
        'FCC6h'
        'FTT8h'
        'FTT10h'
        'TTP9h'
        'TTP7h'
        'CCP5h'
        'CCP3h'
        'CCP1h'
        'CCP2h'
        'CCP4h'
        'CCP6h'
        'TTP8h'
        'TTP10h'
        'TPP9h'
        'TPP7h'
        'CPP5h'
        'CPP3h'
        'CPP1h'
        'CPP2h'
        'CPP4h'
        'CPP6h'
        'TPP8h'
        'TPP10h'
        'PPO9h'
        'PPO7h'
        'PPO5h'
        'PPO3h'
        'PPO1h'
        'PPO2h'
        'PPO4h'
        'PPO6h'
        'PPO8h'
        'PPO10h'
        'POO9h'
        'POO7h'
        'POO5h'
        'POO3h'
        'POO1h'
        'POO2h'
        'POO4h'
        'POO6h'
        'POO8h'
        'POO10h'
        'OI1h'
        'OI2h'
        'Fp1h'
        'Fp2h'
        'AF9h'
        'AF7h'
        'AF5h'
        'AF3h'
        'AF1h'
        'AF2h'
        'AF4h'
        'AF6h'
        'AF8h'
        'AF10h'
        'F9h'
        'F7h'
        'F5h'
        'F3h'
        'F1h'
        'F2h'
        'F4h'
        'F6h'
        'F8h'
        'F10h'
        'FT9h'
        'FT7h'
        'FC5h'
        'FC3h'
        'FC1h'
        'FC2h'
        'FC4h'
        'FC6h'
        'FT8h'
        'FT10h'
        'T9h'
        'T7h'
        'C5h'
        'C3h'
        'C1h'
        'C2h'
        'C4h'
        'C6h'
        'T8h'
        'T10h'
        'TP9h'
        'TP7h'
        'CP5h'
        'CP3h'
        'CP1h'
        'CP2h'
        'CP4h'
        'CP6h'
        'TP8h'
        'TP10h'
        'P9h'
        'P7h'
        'P5h'
        'P3h'
        'P1h'
        'P2h'
        'P4h'
        'P6h'
        'P8h'
        'P10h'
        'PO9h'
        'PO7h'
        'PO5h'
        'PO3h'
        'PO1h'
        'PO2h'
        'PO4h'
        'PO6h'
        'PO8h'
        'PO10h'
        'O1h'
        'O2h'
        'I1h'
        'I2h'
        'AFp9'
        'AFp7'
        'AFp5'
        'AFp3'
        'AFp1'
        'AFpz'
        'AFp2'
        'AFp4'
        'AFp6'
        'AFp8'
        'AFp10'
        'AFF9'
        'AFF7'
        'AFF5'
        'AFF3'
        'AFF1'
        'AFFz'
        'AFF2'
        'AFF4'
        'AFF6'
        'AFF8'
        'AFF10'
        'FFT9'
        'FFT7'
        'FFC5'
        'FFC3'
        'FFC1'
        'FFCz'
        'FFC2'
        'FFC4'
        'FFC6'
        'FFT8'
        'FFT10'
        'FTT9'
        'FTT7'
        'FCC5'
        'FCC3'
        'FCC1'
        'FCCz'
        'FCC2'
        'FCC4'
        'FCC6'
        'FTT8'
        'FTT10'
        'TTP9'
        'TTP7'
        'CCP5'
        'CCP3'
        'CCP1'
        'CCPz'
        'CCP2'
        'CCP4'
        'CCP6'
        'TTP8'
        'TTP10'
        'TPP9'
        'TPP7'
        'CPP5'
        'CPP3'
        'CPP1'
        'CPPz'
        'CPP2'
        'CPP4'
        'CPP6'
        'TPP8'
        'TPP10'
        'PPO9'
        'PPO7'
        'PPO5'
        'PPO3'
        'PPO1'
        'PPOz'
        'PPO2'
        'PPO4'
        'PPO6'
        'PPO8'
        'PPO10'
        'POO9'
        'POO7'
        'POO5'
        'POO3'
        'POO1'
        'POOz'
        'POO2'
        'POO4'
        'POO6'
        'POO8'
        'POO10'
        'OI1'
        'OIz'
        'OI2'
        };
      
      % Add also reference and some alternative labels that might be used
      label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');
      
      % This is to account for all variants of case in 1020 systems
      label = unique(cat(1, label, upper(label), lower(label)));
      
    case 'biosemi64'
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
        };
      
    case 'biosemi128'
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
     