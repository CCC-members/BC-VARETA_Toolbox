function chantype = ft_chantype(input, desired)

% FT_CHANTYPE determines for each individual channel what chantype of data it
% represents, e.g. a planar gradiometer, axial gradiometer, magnetometer,
% trigger channel, etc. If you want to know what the acquisition system is
% (e.g. ctf151 or neuromag306), you should not use this function but
% FT_SENSTYPE instead.
%
% Use as
%   type = ft_chantype(hdr)
%   type = ft_chantype(sens)
%   type = ft_chantype(label)
% or as
%   type = ft_chantype(hdr,   desired)
%   type = ft_chantype(sens,  desired)
%   type = ft_chantype(label, desired)
%
% If the desired unit is not specified as second input argument, this
% function returns a Nchan*1 cell-array with a string describing the type
% of each channel.
%
% If the desired unit is specified as second input argument, this function
% returns a Nchan*1 boolean vector with "true" for the channels of the
% desired type and "false" for the ones that do not match.
%
% The specification of the channel types depends on the acquisition system,
% for example the ctf275 system includes the following type of channels:
% meggrad, refmag, refgrad, adc, trigger, eeg, headloc, headloc_gof.
%
% See also FT_READ_HEADER, FT_SENSTYPE, FT_CHANUNIT

% Copyright (C) 2008-2015, Robert Oostenveld
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

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

% this is to avoid a recursion loop
persistent recursion
if isempty(recursion)
  recursion = false;
end

if nargin<2
  desired = [];
end

% determine the type of input, this is handled similarly as in FT_CHANUNIT
isheader = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'Fs');
isdata   = isa(input, 'struct')  && ~isheader && (isfield(input, 'hdr') || isfield(input, 'grad') || isfield(input, 'elec') || isfield(input, 'opto'));
isgrad   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  &&  isfield(input, 'ori'); % old style
iselec   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  && ~isfield(input, 'ori'); % old style
isgrad   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'coilpos')) || isgrad;             % new style
iselec   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'elecpos')) || iselec;             % new style
isopto   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'transceiver');
islabel  = isa(input, 'cell')    && ~isempty(input) && isa(input{1}, 'char');

if isheader
  % this speeds up the caching in real-time applications
  input.nSamples = 0;
end

current_argin = {input, desired};
if isequal(current_argin, previous_argin)
  % don't do the chantype detection again, but return the previous output from cache
  chantype = previous_argout{1};
  return
end

if isdata
  % the hdr, grad, elec or opto structure might have a different set of channels
  origlabel = input.label;

  if isfield(input, 'hdr')
    input = input.hdr;
    isheader = true;
  elseif isfield(input, 'grad')
    input = input.grad;
    isgrad = true;
  elseif isfield(input, 'elec')
    input = input.elec;
    iselec = true;
  elseif isfield(input, 'opto')
    input = input.opto;
    isopto = true;
  else
    % at least it contains channel labels
    islabel = true;
  end
end

if isheader
  label = input.label;
  numchan = length(label);
elseif isgrad
  label   = input.label;
  numchan = length(label);
elseif iselec
  label   = input.label;
  numchan = length(label);
elseif islabel
  label   = input;
  numchan = length(label);
elseif isfield(input, 'label')
  % this is a last resort: I don't know what it is, but perhaps the labels are informative
  label   = input.label;
  numchan = length(label);
else
  ft_error('the input that was provided to this function cannot be deciphered');
end

if isfield(input, 'chantype')
  % start with the provided channel types
  chantype = input.chantype(:);
else
  % start with unknown chantype for all channels
  chantype = repmat({'unknown'}, numchan, 1);
end

if ~any(strcmp(chantype, 'unknown'))
  % all channels are known, don't bother doing any further heuristics
  
elseif ft_senstype(input, 'unknown')
  % don't bother doing subsequent checks to determine the chantype
  
elseif isheader && (ft_senstype(input, 'neuromag') || ft_senstype(input, 'babysquid74'))
  % channames-KI is the channel kind, 1=meg, 202=eog, 2=eeg, 3=trigger (I am not sure, but have inferred this from a single test file)
  % chaninfo-TY is the Coil chantype (0=magnetometer, 1=planar gradiometer)
  if isfield(input, 'orig') && isfield(input.orig, 'channames')
    for sel=find(input.orig.channames.KI(:)==202)'
      chantype{sel} = 'eog';
    end
    for sel=find(input.orig.channames.KI(:)==2)'
      chantype{sel} = 'eeg';
    end
    for sel=find(input.orig.channames.KI(:)==3)'
      chantype{sel} = 'digital trigger';
    end
    % determine the MEG channel subtype
    selmeg=find(input.orig.channames.KI(:)==1)';
    for i=1:length(selmeg)
      if input.orig.chaninfo.TY(i)==0
        chantype{selmeg(i)} = 'megmag';
      elseif input.orig.chaninfo.TY(i)==1
        % FIXME this might also be a axial gradiometer in case the BabySQUID data is read with the old reading routines
        chantype{selmeg(i)} = 'megplanar';
      end
    end

  elseif isfield(input, 'orig') && isfield(input.orig, 'chs') && isfield(input.orig.chs, 'coil_type')
    % all the chs.kinds and chs.coil_types are obtained from the MNE manual, p.210-211
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==2)' % planar gradiometers
      chantype(sel) = {'megplanar'}; % Neuromag-122 planar gradiometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3012)' %planar gradiometers
      chantype(sel) = {'megplanar'}; % Type T1 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3013)' %planar gradiometers
      chantype(sel) = {'megplanar'}; % Type T2 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3014)' %planar gradiometers
      chantype(sel) = {'megplanar'}; % Type T3 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3022)' %magnetometers
      chantype(sel) = {'megmag'};    % Type T1 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3023)' %magnetometers
      chantype(sel) = {'megmag'};    % Type T2 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3024)' % magnetometers
      chantype(sel) = {'megmag'};    % Type T3 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==7001)' % axial gradiometer
      chantype(sel) = {'megaxial'};
    end
    for sel=find([input.orig.chs.kind]==301)' % MEG reference channel, located far from head
      chantype(sel) = {'ref'};
    end
    for sel=find([input.orig.chs.kind]==2)'   % EEG channels
      chantype(sel) = {'eeg'};
    end
    for sel=find([input.orig.chs.kind]==201)' % MCG channels
      chantype(sel) = {'mcg'};
    end
    for sel=find([input.orig.chs.kind]==3)' % Stim channels
      if any([input.orig.chs(sel).logno] == 101) % new systems: 101 (and 102, if enabled) are digital; low numbers are 'pseudo-analog' (if enabled)
        chantype(sel([input.orig.chs(sel).logno] == 101)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] == 102)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] <= 32))  = {'analog trigger'};
        others = [input.orig.chs(sel).logno] > 32 & [input.orig.chs(sel).logno] ~= 101 & ...
          [input.orig.chs(sel).logno] ~= 102;
        chantype(sel(others)) = {'other trigger'};
      elseif any(ismember([input.orig.chs(sel).logno], [14 15 16])) % older systems: STI 014/015/016 are digital; lower numbers 'pseudo-analog'(if enabled)
        chantype(sel([input.orig.chs(sel).logno] == 14)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] == 15)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] == 16)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] <= 13)) = {'analog trigger'};
        others = [input.orig.chs(sel).logno] > 16;
        chantype(sel(others)) = {'other trigger'};
      else
        ft_warning('There does not seem to be a suitable trigger channel.');
        chantype(sel) = {'other trigger'};
      end
    end
    for sel=find([input.orig.chs.kind]==202)' % EOG
      chantype(sel) = {'eog'};
    end
    for sel=find([input.orig.chs.kind]==302)' % EMG
      chantype(sel) = {'emg'};
    end
    for sel=find([input.orig.chs.kind]==402)' % ECG
      chantype(sel) = {'ecg'};
    end
    for sel=find([input.orig.chs.kind]==502)' % MISC
      chantype(sel) = {'misc'};
    end
    for sel=find([input.orig.chs.kind]==602)' % Resp
      chantype(sel) = {'respiration'};
    end
  end

elseif ft_senstype(input, 'babysquid74')
  % the name can be something like "MEG 001" or "MEG001" or "MEG 0113" or "MEG0113"
  % i.e. with two or three digits and with or without a space
  sel = myregexp('^MEG', label);
  chantype(sel) = {'megaxial'};

elseif ft_senstype(input, 'neuromag122')
  % the name can be something like "MEG 001" or "MEG001" or "MEG 0113" or "MEG0113"
  % i.e. with two or three digits and with or without a space
  sel = myregexp('^MEG', label);
  chantype(sel) = {'megplanar'};

elseif ft_senstype(input, 'neuromag306') && isgrad
  % there should be 204 planar gradiometers and 102 axial magnetometers
  if isfield(input, 'tra')
    tmp = sum(abs(input.tra)>0,2);
    sel = (tmp==1);
    chantype(sel) = {'megmag'};
    sel = (tmp==2);
    chantype(sel) = {'megplanar'};
  end

elseif ft_senstype(input, 'neuromag306') && islabel
  sel = myregexp('^MEG.*1$', label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^MEG.*2$', label);
  chantype(sel) = {'megplanar'};
  sel = myregexp('^MEG.*3$', label);
  chantype(sel) = {'megplanar'};

elseif ft_senstype(input, 'neuromag306_combined') && islabel
  % the magnetometers are detected, the combined channels remain unknown
  sel = myregexp('^MEG.*1$', label);
  chantype(sel) = {'megmag'};

elseif ft_senstype(input, 'ctf') && isheader
  % The following is according to "CTF MEG(TM) File Formats" pdf, Release 5.2.1
  %
  % eMEGReference      0 Reference magnetometer channel
  % eMEGReference1     1 Reference 1st-order gradiometer channel
  % eMEGReference2     2 Reference 2nd-order gradiometer channel
  % eMEGReference3     3 Reference 3rd-order gradiometer channel
  % eMEGSensor         4 Sensor magnetometer channel located in head shell
  % eMEGSensor1        5 Sensor 1st-order gradiometer channel located in head shell
  % eMEGSensor2        6 Sensor 2nd-order gradiometer channel located in head shell
  % eMEGSensor3        7 Sensor 3rd-order gradiometer channel located in head shell
  % eEEGRef            8 EEG unipolar sensors not on the scalp
  % eEEGSensor         9 EEG unipolar sensors on the scalp
  % eADCRef           10 (see eADCAmpRef below)
  % eADCAmpRef        10 ADC amp channels from HLU or PIU (old electronics)
  % eStimRef          11 Stimulus channel for MEG41
  % eTimeRef          12 Time reference coming from video channel
  % ePositionRef      13 Not used
  % eDACRef           14 DAC channel from ECC or HLU
  % eSAMSensor        15 SAM channel derived through data analysis
  % eVirtualSensor    16 Virtual channel derived by combining two or more physical channels
  % eSystemTimeRef    17 System time showing elapsed time since trial started
  % eADCVoltRef       18 ADC volt channels from ECC
  % eStimAnalog       19 Analog trigger channels
  % eStimDigital      20 Digital trigger channels
  % eEEGBipolar       21 EEG bipolar sensor not on the scalp
  % eEEGAflg          22 EEG ADC over range flags
  % eMEGReset         23 MEG resets (counts sensor jumps for crosstalk purposes)
  % eDipSrc           24 Dipole source
  % eSAMSensorNorm    25 Normalized SAM channel derived through data analy- sis
  % eAngleRef         26 Orientation of head localization field
  % eExtractionRef    27 Extracted signal from each sensor of field generated by each localization coil
  % eFitErr           28 Fit error from each head localization coil
  % eOtherRef         29 Any other type of sensor not mentioned but still valid
  % eInvalidType      30 An invalid sensor

  % start with an empty one
  origSensType = [];
  if isfield(input, 'orig')
    if isfield(input.orig, 'sensType') && isfield(input.orig, 'Chan')
      % the header was read using the open-source MATLAB code that originates from CTF and that was modified by the FCDC
      origSensType = input.orig.sensType;
    elseif isfield(input.orig, 'res4') && isfield(input.orig.res4, 'senres')
      % the header was read using the CTF p-files, i.e. readCTFds
      origSensType =  [input.orig.res4.senres.sensorTypeIndex];
    elseif isfield(input.orig, 'sensor') && isfield(input.orig.sensor, 'info')
      % the header was read using the CTF importer from the NIH and Daren Weber
      origSensType = [input.orig.sensor.info.index];
    end
  end

  if isempty(origSensType)
    ft_warning('could not determine channel chantype from the CTF header');
  end

  for sel=find(origSensType(:)==0)'
    chantype{sel} = 'refmag';
  end
  for sel=find(origSensType(:)==1)'
    chantype{sel} = 'refgrad';
  end
  for sel=find(origSensType(:)==5)'
    chantype{sel} = 'meggrad';
  end
  for sel=find(origSensType(:)==9)'
    chantype{sel} = 'eeg';
  end
  for sel=find(origSensType(:)==11)'
    % Stimulus channel for MEG41
    chantype{sel} = 'trigger';
  end
  for sel=find(origSensType(:)==13)'
    chantype{sel} = 'headloc'; % these represent the x, y, z position of the head coils
  end
  for sel=find(origSensType(:)==17)'
    chantype{sel} = 'clock';
  end
  for sel=find(origSensType(:)==18)'
    chantype{sel} = 'adc';
  end
  for sel=find(origSensType(:)==20)'
    % Digital trigger channels
    chantype{sel} = 'trigger';
  end
  for sel=find(origSensType(:)==28)'
    chantype{sel} = 'headloc_gof'; % these represent the goodness of fit for the head coils
  end
  for sel=find(origSensType(:)==29)'
    chantype{sel} = 'reserved'; % these are "reserved for future use", but relate to head localization
  end

elseif ft_senstype(input, 'ctf') && isgrad
  % in principle it is possible to look at the number of coils, but here the channels are identified based on their name
  sel = myregexp('^M[ZLR][A-Z][0-9][0-9]$', input.label);
  chantype(sel) = {'meggrad'};            % normal gradiometer channels
  sel = myregexp('^S[LR][0-9][0-9]$', input.label);
  chantype(sel) = {'meggrad'};            % normal gradiometer channels in the 64 channel CTF system
  sel = myregexp('^B[GPQR][0-9]$', input.label);
  chantype(sel) = {'refmag'};             % reference magnetometers
  sel = myregexp('^[GPQR][0-9][0-9]$', input.label);
  chantype(sel) = {'refgrad'};            % reference gradiometers

elseif ft_senstype(input, 'ctf') && islabel
  % the channels have to be identified based on their name alone
  sel = myregexp('^M[ZLR][A-Z][0-9][0-9]$', label);
  chantype(sel) = {'meggrad'};            % normal gradiometer channels
  sel = myregexp('^S[LR][0-9][0-9]$', label);
  chantype(sel) = {'meggrad'};            % normal gradiometer channels in the 64 channel CTF system
  sel = myregexp('^B[GPR][0-9]$', label);
  chantype(sel) = {'refmag'};             % reference magnetometers
  sel = myregexp('^[GPQR][0-9][0-9]$', label);
  chantype(sel) = {'refgrad'};            % reference gradiometers

elseif ft_senstype(input, 'bti')
  if isfield(input, 'orig') && isfield(input.orig, 'config')
    configname = {input.orig.config.channel_data.name};
    configtype = [input.orig.config.channel_data.type];

    if ~isequal(configname(:), input.label(:))
      % reorder the channels according to the order in input.label
      [sel1, sel2] = match_str(input.label, configname);
      configname = configname(sel2);
      configtype = configtype(sel2);
      configdata = input.orig.config.channel_data(sel2);
    end
    numloops = zeros(size(configdata));
    for i=1:length(configdata)
      if isfield(configdata(i).device_data, 'total_loops')
        numloops(i) = configdata(i).device_data.total_loops;
      end
    end

    % these are taken from bti2grad
    chantype(configtype==1 & numloops==1) = {'megmag'};
    chantype(configtype==1 & numloops==2) = {'meggrad'};
    chantype(configtype==2) = {'eeg'};
    chantype(configtype==3) = {'ref'}; % not known if mag or grad
    chantype(configtype==4) = {'aux'};
    chantype(configtype==5) = {'trigger'};

    % refine the distinction between refmag and refgrad to make the types
    % in grad and header consistent
    sel = myregexp('^M[CLR][xyz][aA]*$', label);
    chantype(sel) = {'refmag'};
    sel = myregexp('^G[xyz][xyz]A$', label);
    chantype(sel) = {'refgrad'};
  else
    % determine the chantype on the basis of the channel labels
    % all 4D-BTi MEG channels start with "A" followed by a number
    % all 4D-BTi reference channels start with M or G
    % all 4D-BTi EEG channels start with E, except for the 248-MEG/32-EEG system in Warsaw where they end with -1
    sel = myregexp('^A[0-9]+$', label);
    chantype(sel) = {'meg'};
    sel = myregexp('^M[CLR][xyz][aA]*$', label);
    chantype(sel) = {'refmag'};
    sel = myregexp('^G[xyz][xyz]A$', label);
    chantype(sel) = {'refgrad'};

    if isgrad && isfield(input, 'tra')
      gradtype = repmat({'unknown'}, size(input.label));
      gradtype(strncmp('A', input.label, 1)) = {'meg'};
      gradtype(strncmp('M', input.label, 1)) = {'refmag'};
      gradtype(strncmp('G', input.label, 1)) = {'refgrad'};
      % look at the number of coils of the meg channels
      selchan = find(strcmp('meg', gradtype));
      for k = 1:length(selchan)
        ncoils = length(find(input.tra(selchan(k),:)==1));
        if ncoils==1
          gradtype{selchan(k)} = 'megmag';
        elseif ncoils==2
          gradtype{selchan(k)} = 'meggrad';
        end
      end
      [selchan, selgrad] = match_str(label, input.label);
      chantype(selchan) = gradtype(selgrad);
    end

    % deal with additional channel types based on the names
    if isheader && issubfield(input, 'orig.channel_data.chan_label')
      tmplabel = {input.orig.channel_data.chan_label};
      tmplabel = tmplabel(:);
    else
      tmplabel = label; % might work
    end
    sel = find(strcmp('unknown', chantype));
    if ~isempty(sel)
      chantype(sel) = ft_chantype(tmplabel(sel));
      sel       = find(strcmp('unknown', chantype));
      if ~isempty(sel)
        % channels that start with E are assumed to be EEG
        % channels that end with -1 are also assumed to be EEG, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2389
        chantype(sel(cellfun(@(x) strcmp(x(end-1:end),'-1') || strcmp(x(1),'E'), label(sel)))) = {'eeg'};
      end
    end
  end

elseif ft_senstype(input, 'itab') && isheader
  origtype = [input.orig.ch.type];
  chantype(origtype==0) = {'unknown'};
  chantype(origtype==1) = {'ele'};
  chantype(origtype==2) = {'mag'}; % might be magnetometer or gradiometer, look at the number of coils
  chantype(origtype==4) = {'ele ref'};
  chantype(origtype==8) = {'mag ref'};
  chantype(origtype==16) = {'aux'};
  chantype(origtype==32) = {'param'};
  chantype(origtype==64) = {'digit'};
  chantype(origtype==128) = {'flag'};
  % these are the channels that are visible to FieldTrip
  chansel = 1:input.orig.nchan;
  chantype = chantype(chansel);

elseif ft_senstype(input, 'yokogawa') && isheader
  % This is to recognize Yokogawa channel types from the original header
  % This is from the original documentation
  NullChannel                   = 0;
  MagnetoMeter                  = 1;
  AxialGradioMeter              = 2;
  PlannerGradioMeter            = 3;
  RefferenceChannelMark         = hex2dec('0100');
  RefferenceMagnetoMeter        = bitor( RefferenceChannelMark, MagnetoMeter      );
  RefferenceAxialGradioMeter    = bitor( RefferenceChannelMark, AxialGradioMeter  );
  RefferencePlannerGradioMeter  = bitor( RefferenceChannelMark, PlannerGradioMeter);
  TriggerChannel                = -1;
  EegChannel                    = -2;
  EcgChannel                    = -3;
  EtcChannel                    = -4;
  if ft_hastoolbox('yokogawa_meg_reader')
    label = input.label;
    sel = myregexp('[0-9][0-9][0-9]$', label);
    chantype(sel) = {'null'};
    sel = myregexp('^M[0-9][0-9][0-9]$', label);
    chantype(sel) = {'megmag'};
    sel = myregexp('^AG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'meggrad'};
    sel = myregexp('^PG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'megplanar'};
    sel = myregexp('^RM[0-9][0-9][0-9]$', label);
    chantype(sel) = {'refmag'};
    sel = myregexp('^RAG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'refgrad'};
    sel = myregexp('^RPG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'refplanar'};
    sel = myregexp('^TRIG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'trigger'};
    %% Possible labels categorized in "eeg"
    sel_A = myregexp('^A[^G]*[0-9hzZ]$', label);
    sel_P = myregexp('^P[^G]*[0-9hzZ]$', label);
    sel_T = myregexp('^T[^R]*[0-9hzZ]$', label);
    sel_E = myregexp('^E$', label);
    sel_Z = myregexp('^[zZ]$', label);
    sel_M = myregexp('^M[0-9]$', label);
    sel_O = myregexp('^[BCFION]\w*[0-9hzZ]$', label);
    sel_EEG = myregexp('^EEG[0-9][0-9][0-9]$', label);
    sel = logical( sel_A + sel_P + sel_T + sel_E + sel_Z + sel_M + sel_O + sel_EEG );
    clear sel_A sel_P sel_T sel_E sel_Z sel_M sel_O sel_EEG
    chantype(sel) = {'eeg'};
    %% Additional EOG, ECG labels
    sel = myregexp('^EO[0-9]$', label); % EO
    chantype(sel) = {'eog'};
%    sel = myregexp('^ECG[0-9][0-9][0-9]$', label);
    sel_X = myregexp('^X[0-9]$', label); % X
    sel_ECG = myregexp('^ECG[0-9][0-9][0-9]$', label);
    sel = logical( sel_X + sel_ECG );
    clear sel_X sel_ECG
    chantype(sel) = {'ecg'};
    sel = myregexp('^ETC[0-9][0-9][0-9]$', label);
    chantype(sel) = {'etc'};

%   % shorten names
%    ch_info = input.orig.channel_info.channel;
%    type_orig = [ch_info.type];
%    sel = (type_orig == NullChannel);
%    chantype(sel) = {'null'};
%    sel = (type_orig == MagnetoMeter);
%    chantype(sel) = {'megmag'};
%    sel = (type_orig == AxialGradioMeter);
%    chantype(sel) = {'meggrad'};
%    sel = (type_orig == PlannerGradioMeter);
%    chantype(sel) = {'megplanar'};
%    sel = (type_orig == RefferenceMagnetoMeter);
%    chantype(sel) = {'refmag'};
%    sel = (type_orig == RefferenceAxialGradioMeter);
%    chantype(sel) = {'refgrad'};
%    sel = (type_orig == RefferencePlannerGradioMeter);
%    chantype(sel) = {'refplanar'};
%    sel = (type_orig == TriggerChannel);
%    chantype(sel) = {'trigger'};
%    sel = (type_orig == EegChannel);
%    chantype(sel) = {'eeg'};
%    sel = (type_orig == EcgChannel);
%    chantype(sel) = {'ecg'};
%    sel = (type_orig == EtcChannel);
%    chantype(sel) = {'etc'};

  elseif ft_hastoolbox('yokogawa')
    sel = (input.orig.channel_info(:, 2) == NullChannel);
    chantype(sel) = {'null'};
    sel = (input.orig.channel_info(:, 2) == MagnetoMeter);
    chantype(sel) = {'megmag'};
    sel = (input.orig.channel_info(:, 2) == AxialGradioMeter);
    chantype(sel) = {'meggrad'};
    sel = (input.orig.channel_info(:, 2) == PlannerGradioMeter);
    chantype(sel) = {'megplanar'};
    sel = (input.orig.channel_info(:, 2) == RefferenceMagnetoMeter);
    chantype(sel) = {'refmag'};
    sel = (input.orig.channel_info(:, 2) == RefferenceAxialGradioMeter);
    chantype(sel) = {'refgrad'};
    sel = (input.orig.channel_info(:, 2) == RefferencePlannerGradioMeter);
    chantype(sel) = {'refplanar'};
    sel = (input.orig.channel_info(:, 2) == TriggerChannel);
    chantype(sel) = {'trigger'};
    sel = (input.orig.channel_info(:, 2) == EegChannel);
    chantype(sel) = {'eeg'};
    sel = (input.orig.channel_info(:, 2) == EcgChannel);
    chantype(sel) = {'ecg'};
    sel = (input.orig.channel_info(:, 2) == EtcChannel);
    chantype(sel) = {'etc'};
  end

elseif ft_senstype(input, 'yokogawa') && isgrad
  % all channels in the gradiometer definition are meg
  % chantype(1:end) = {'meg'};
  % channels are identified based on their name: only magnetic as isgrad==1
  sel = myregexp('^M[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^AG[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'meggrad'};
  sel = myregexp('^PG[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'megplanar'};
  sel = myregexp('^RM[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'refmag'};
  sel = myregexp('^RAG[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'refgrad'};
  sel = myregexp('^RPG[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'refplanar'};

elseif ft_senstype(input, 'yokogawa') && islabel
  % the yokogawa channel labels are a mess, so autodetection is not possible
  % chantype(1:end) = {'meg'};
  sel = myregexp('[0-9][0-9][0-9]$', label);
  chantype(sel) = {'null'};
  sel = myregexp('^M[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^AG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'meggrad'};
  sel = myregexp('^PG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megplanar'};
  sel = myregexp('^RM[0-9][0-9][0-9]$', label);
  chantype(sel) = {'refmag'};
  sel = myregexp('^RAG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'refgrad'};
  sel = myregexp('^RPG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'refplanar'};
  sel = myregexp('^TRIG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'trigger'};
  %% Possible labels categorized in "eeg"
  sel_A = myregexp('^A[^G]*[0-9hzZ]$', label);
  sel_P = myregexp('^P[^G]*[0-9hzZ]$', label);
  sel_T = myregexp('^T[^R]*[0-9hzZ]$', label);
  sel_E = myregexp('^E$', label);
  sel_Z = myregexp('^[zZ]$', label);
  sel_M = myregexp('^M[0-9]$', label);
  sel_O = myregexp('^[BCFION]\w*[0-9hzZ]$', label);
  sel_EEG = myregexp('^EEG[0-9][0-9][0-9]$', label);
  sel = logical( sel_A + sel_P + sel_T + sel_E + sel_Z + sel_M + sel_O + sel_EEG );
  clear sel_A sel_P sel_T sel_E sel_Z sel_M sel_O sel_EEG
  chantype(sel) = {'eeg'};
  %% Additional EOG, ECG labels
  sel = myregexp('^EO[0-9]$', label); % EO
  chantype(sel) = {'eog'};
% sel = myregexp('^ECG[0-9][0-9][0-9]$', label);
  sel_X = myregexp('^X[0-9]$', label); % X
  sel_ECG = myregexp('^ECG[0-9][0-9][0-9]$', label);
  sel = logical( sel_X + sel_ECG );
  clear sel_X sel_ECG
  chantype(sel) = {'ecg'};
  sel = myregexp('^ETC[0-9][0-9][0-9]$', label);
  chantype(sel) = {'etc'};

elseif ft_senstype(input, 'itab') && isheader
  sel = ([input.orig.ch.type]==0);
  chantype(sel) = {'unknown'};
  sel = ([input.orig.ch.type]==1);
  chantype(sel) = {'unknown'};
  sel = ([input.orig.ch.type]==2);
  chantype(sel) = {'megmag'};
  sel = ([input.orig.ch.type]==8);
  chantype(sel) = {'megref'};
  sel = ([input.orig.ch.type]==16);
  chantype(sel) = {'aux'};
  sel = ([input.orig.ch.type]==64);
  chantype(sel) = {'digital'};
  % not all channels are actually processed by FieldTrip, so only return
  % the types fopr the ones that read_header and read_data return
  chantype = chantype(input.orig.chansel);

elseif ft_senstype(input, 'itab') && isgrad
  % the channels have to be identified based on their name alone
  sel = myregexp('^MAG_[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^MAG_[0-9][0-9]$', label); % for the itab28 system
  chantype(sel) = {'megmag'};
  sel = myregexp('^MAG_[0-9]$', label); % for the itab28 system
  chantype(sel) = {'megmag'};
  sel = myregexp('^REF_[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megref'};
  sel = myregexp('^AUX.*$', label);
  chantype(sel) = {'aux'};

elseif ft_senstype(input, 'itab') && islabel
  % the channels have to be identified based on their name alone
  sel = myregexp('^MAG_[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^REF_[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megref'};
  sel = myregexp('^AUX.*$', label);
  chantype(sel) = {'aux'};

elseif ft_senstype(input, 'eeg') && islabel
  % use an external helper function to define the list with EEG channel names
  chantype(match_str(label, ft_senslabel('eeg1005'))) = {'eeg'};          % this includes all channels from the 1010 and 1020 arrangement
  chantype(match_str(label, ft_senslabel(ft_senstype(input)))) = {'eeg'}; % this will work for biosemi, egi and other detected channel arrangements

elseif ft_senstype(input, 'eeg') && iselec
  % all channels in an electrode definition must be eeg channels
  chantype(:) = {'eeg'};

elseif ft_senstype(input, 'eeg') && isheader
  % use an external helper function to define the list with EEG channel names
  chantype(match_str(input.label, ft_senslabel(ft_senstype(input)))) = {'eeg'};

elseif ft_senstype(input, 'plexon') && isheader
  % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
  for i=1:numchan
    switch input.orig.VarHeader(i).Type
      case 0
        chantype{i} = 'spike';
      case 1
        chantype{i} = 'event';
      case 2
        chantype{i} = 'interval';  % Interval variables?
      case 3
        chantype{i} = 'waveform';
      case 4
        chantype{i} = 'population'; % Population variables ?
      case 5
        chantype{i} = 'analog';
      otherwise
        % keep the default 'unknown' chantype
    end
  end

end % ft_senstype

% if possible, set additional types based on channel labels
label2type = {
  {'ecg', 'ekg'};
  {'emg'};
  {'eog', 'heog', 'veog'};
  {'lfp'};
  {'eeg'};
  {'trigger', 'trig', 'dtrig'};
  };
for i = 1:numel(label2type)
  for j = 1:numel(label2type{i})
    chantype(intersect(strmatch(label2type{i}{j}, lower(label)), find(strcmp(chantype, 'unknown')))) = label2type{i}(1);
  end
end

if isdata
  % the input was replaced by one of hdr, grad, elec, opto
  [sel1, sel2] = match_str(origlabel, input.label);
  origtype = repmat({'unknown'}, size(origlabel));
  origtype(sel1) = chantype(sel2);
  % the hdr, grad, elec or opto structure might have a different set of channels
  chantype = origtype;
end

if all(strcmp(chantype, 'unknown')) && ~recursion
  % try whether only lowercase channel labels makes a difference
  if islabel
    recursion = true;
    chantype = ft_chantype(lower(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = lower(input.label);
    recursion = true;
    chantype = ft_chantype(input);
    recursion = false;
  end
end

if all(strcmp(chantype, 'unknown')) && ~recursion
  % try whether only uppercase channel labels makes a difference
  if islabel
    recursion = true;
    chantype = ft_chantype(upper(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = upper(input.label);
    recursion = true;
    chantype = ft_chantype(input);
    recursion = false;
  end
end

if nargin>1
  % return a boolean vector
  if isequal(desired, 'meg') || isequal(desired, 'ref')
    % only compare the first three characters, i.e. meggrad or megmag should match
    chantype = strncmp(desired, chantype, 3);
  else
    chantype = strcmp(desired, chantype);
  end
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {chantype};
previous_argin  = current_argin;
previous_argout = current_argout;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function match = myregexp(pat, list)
match = false(size(list));
for i=1:numel(list)
  match(i) = ~isempty(regexp(list{i}, pat, 'once'));
end
                                                                                                                                  te cross-arrow directions
basecross  = crossdir + basepoint;
tipcross   = crossdir + tippoint;
sbasecross = crossdir + sbasepoint;
stipcross  = crossdir + stippoint;
ii = find(all(crossdir==0)|any(isnan(crossdir)));
if ~isempty(ii),
    numii = length(ii);
    %   transform start points
        tmp1 = [basepoint(:,ii) tippoint(:,ii) sbasepoint(:,ii) stippoint(:,ii)];
        tmp1 = (tmp1-axm(:,[ii ii ii ii])) ./ axr(:,[ii ii ii ii]);
        tmp1 = [tmp1; ones(1,4*numii)];
        if (oneax), X0=T*tmp1;
        else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T(:,[ii ii ii ii]).*tmp1;
              tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
              X0=zeros(4,4*numii); X0(:)=sum(tmp2)'; end;
        X0=X0./(ones(4,1)*X0(4,:));
    %   transform stop points
        tmp1 = [(2*stop(:,ii)-start(:,ii)-axm(:,ii))./axr(:,ii);ones(1,numii)];
        tmp1 = [tmp1 tmp1 tmp1 tmp1];
        if (oneax), Xf=T*tmp1;
        else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T(:,[ii ii ii ii]).*tmp1;
              tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
              Xf=zeros(4,4*numii); Xf(:)=sum(tmp2)'; end;
        Xf=Xf./(ones(4,1)*Xf(4,:));
    %   compute perpendicular directions
        pixfact = ((limrange(1,ii)./limrange(2,ii)).*(ap(2,ii)./ap(1,ii))).^2;
        pixfact = [pixfact pixfact pixfact pixfact];
        pixfact = [pixfact;1./pixfact];
        [dummyval,jj] = max(abs(Xf(1:2,:)-X0(1:2,:)));
        jj1 = ((1:4)'*ones(1,length(jj))==ones(4,1)*jj);
        jj2 = ((1:4)'*ones(1,length(jj))==ones(4,1)*(3-jj));
        jj3 = jj1(1:2,:);
        Xf(jj1)=Xf(jj1)+(Xf(jj1)-X0(jj1)==0); %eaj new 2/24/98
        Xp = X0;
        Xp(jj2) = X0(jj2) + ones(sum(jj2(:)),1);
        Xp(jj1) = X0(jj1) - (Xf(jj2)-X0(jj2))./(Xf(jj1)-X0(jj1)) .* pixfact(jj3);
    %   inverse transform the cross points
        if (oneax), Xp=invT*Xp;
        else, tmp1=[Xp;Xp;Xp;Xp]; tmp1=invT(:,[ii ii ii ii]).*tmp1;
              tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
              Xp=zeros(4,4*numii); Xp(:)=sum(tmp2)'; end;
        Xp=(Xp(1:3,:)./(ones(3,1)*Xp(4,:))).*axr(:,[ii ii ii ii])+axm(:,[ii ii ii ii]);
        basecross(:,ii)  = Xp(:,0*numii+(1:numii));
        tipcross(:,ii)   = Xp(:,1*numii+(1:numii));
        sbasecross(:,ii) = Xp(:,2*numii+(1:numii));
        stipcross(:,ii)  = Xp(:,3*numii+(1:numii));
end;

% compute all points
%   compute start points
    axm11 = [axm axm axm axm axm axm axm axm axm axm axm];
    axr11 = [axr axr axr axr axr axr axr axr axr axr axr];
    st = [stoppoint tippoint basepoint sbasepoint stippoint startpoint stippoint sbasepoint basepoint tippoint stoppoint];
    tmp1 = (st - axm11) ./ axr11;
    tmp1 = [tmp1; ones(1,size(tmp1,2))];
    if (oneax), X0=T*tmp1;
    else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[T T T T T T T T T T T].*tmp1;
          tmp2=zeros(4,44*narrows); tmp2(:)=tmp1(:);
          X0=zeros(4,11*narrows); X0(:)=sum(tmp2)'; end;
    X0=X0./(ones(4,1)*X0(4,:));
%   compute stop points
    tmp1 = ([start tipcross basecross sbasecross stipcross stop stipcross sbasecross basecross tipcross start] ...
         - axm11) ./ axr11;
    tmp1 = [tmp1; ones(1,size(tmp1,2))];
    if (oneax), Xf=T*tmp1;
    else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[T T T T T T T T T T T].*tmp1;
          tmp2=zeros(4,44*narrows); tmp2(:)=tmp1(:);
          Xf=zeros(4,11*narrows); Xf(:)=sum(tmp2)'; end;
    Xf=Xf./(ones(4,1)*Xf(4,:));
%   compute lengths
    len0  = len.*((ends==1)|(ends==3)).*tan(tipangle/180*pi);
    slen0 = len.*((ends==2)|(ends==3)).*tan(tipangle/180*pi);
    le = [zeros(1,narrows) len0 wid/2 wid/2 slen0 zeros(1,narrows) -slen0 -wid/2 -wid/2 -len0 zeros(1,narrows)];
    aprange = ap./limrange;
    aprange = [aprange aprange aprange aprange aprange aprange aprange aprange aprange aprange aprange];
    D = sqrt(sum(((Xf(1:2,:)-X0(1:2,:)).*aprange).^2));
    Dii=find(D==0); if ~isempty(Dii), D=D+(D==0); le(Dii)=zeros(1,length(Dii)); end; %should fix DivideByZero warnings
    tmp1 = X0.*(ones(4,1)*(1-le./D)) + Xf.*(ones(4,1)*(le./D));
%   inverse transform
    if (oneax), tmp3=invT*tmp1;
    else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[invT invT invT invT invT invT invT invT invT invT invT].*tmp1;
          tmp2=zeros(4,44*narrows); tmp2(:)=tmp1(:);
          tmp3=zeros(4,11*narrows); tmp3(:)=sum(tmp2)'; end;
    pts = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)) .* axr11 + axm11;

% correct for ones where the crossdir was specified
ii = find(~(all(crossdir==0)|any(isnan(crossdir))));
if ~isempty(ii),
    D1 = [pts(:,1*narrows+ii)-pts(:,9*narrows+ii) ...
          pts(:,2*narrows+ii)-pts(:,8*narrows+ii) ...
          pts(:,3*narrows+ii)-pts(:,7*narrows+ii) ...
          pts(:,4*narrows+ii)-pts(:,6*narrows+ii) ...
          pts(:,6*narrows+ii)-pts(:,4*narrows+ii) ...
          pts(:,7*narrows+ii)-pts(:,3*narrows+ii) ...
          pts(:,8*narrows+ii)-pts(:,2*narrows+ii) ...
          pts(:,9*narrows+ii)-pts(:,1*narrows+ii)]/2;
    ii = ii'*ones(1,8) + ones(length(ii),1)*[1:4 6:9]*narrows;
    ii = ii(:)';
    pts(:,ii) = st(:,ii) + D1;
end;


% readjust for reverse directions
iicols=(1:narrows)'; iicols=iicols(:,ones(1,11)); iicols=iicols(:).';
tmp1=axrev(:,iicols);
ii = find(tmp1(:)); if ~isempty(ii), pts(ii)=-pts(ii); end;

% readjust for log scale on axes
tmp1=xyzlog(:,iicols);
ii = find(tmp1(:)); if ~isempty(ii), pts(ii)=10.^pts(ii); end;

% compute the x,y,z coordinates of the patches;
ii = narrows*(0:10)'*ones(1,narrows) + ones(11,1)*(1:narrows);
ii = ii(:)';
x = zeros(11,narrows);
y = zeros(11,narrows);
z = zeros(11,narrows);
x(:) = pts(1,ii)';
y(:) = pts(2,ii)';
z(:) = pts(3,ii)';

% do the output
if (nargout<=1),
%   % create or modify the patches
    newpatch = trueornan(ispatch) & (isempty(oldh)|~strcmp(get(oldh,'Type'),'patch'));
    newline = ~trueornan(ispatch) & (isempty(oldh)|~strcmp(get(oldh,'Type'),'line'));
    if isempty(oldh), H=zeros(narrows,1); else, H=oldh; end;
%   % make or modify the arrows
    for k=1:narrows,
        if all(isnan(ud(k,[3 6])))&arrow_is2DXY(ax(k)), zz=[]; else, zz=z(:,k); end;
        % work around a MATLAB 6.x OpenGL bug -- 7/28/02
          xx=x(:,k); yy=y(:,k); 
          mask=any([ones(1,2+size(zz,2));diff([xx yy zz],[],1)],2);
          xx=xx(mask); yy=yy(mask); if ~isempty(zz), zz=zz(mask); end;
        % plot the patch or line
        xyz = {'XData',xx,'YData',yy,'ZData',zz,'Tag',ArrowTag};
        if newpatch(k)|newline(k),
            if newpatch(k),
                H(k) = patch(xyz{:});
            else,
                H(k) = line(xyz{:});
            end;
            if ~isempty(oldh), arrow_copyprops(oldh(k),H(k)); end;
        else,
            if ispatch(k), xyz={xyz{:},'CData',[]}; end;
            set(H(k),xyz{:});
        end;
    end;
    if ~isempty(oldh), delete(oldh(oldh~=H)); end;
%   % additional properties
    set(H,'Clipping','off');
    set(H,{'UserData'},num2cell(ud,2));
    if (length(extraprops)>0), set(H,extraprops{:}); end;
    % handle choosing arrow Start and/or Stop locations if unspecified
    [H,oldaxlims,errstr] = arrow_clicks(H,ud,x,y,z,ax,oldaxlims);
    if ~isempty(errstr), error([upper(mfilename) ' got ' errstr]); end;
    % set the output
    if (nargout>0), h=H; end;
    % make sure the axis limits did not change
    if isempty(oldaxlims),
        ARROW_AXLIMITS = [];
    else,
        lims = get(oldaxlims(:,1),{'XLim','YLim','ZLim'})';
        lims = reshape(cat(2,lims{:}),6,size(lims,2));
        mask = arrow_is2DXY(oldaxlims(:,1));
        oldaxlims(mask,6:7) = lims(5:6,mask)';
        ARROW_AXLIMITS = oldaxlims(find(any(oldaxlims(:,2:7)'~=lims)),:);
        if ~isempty(ARROW_AXLIMITS),
            warning(arrow_warnlimits(ARROW_AXLIMITS,narrows));
        end;
    end;
else,
    % don't create the patch, just return the data
    h=x;
    yy=y;
    zz=z;
end;



function out = arrow_defcheck(in,def,prop)
% check if we got 'default' values
    out = in;
    if ~ischar(in), return; end;
    if size(in,1)==1 & strncmp(lower(in),'def',3),
        out = def;
    elseif ~isempty(prop),
        error([upper(mfilename) ' does not recognize ''' in(:)' ''' as a valid ''' prop ''' string.']);
    end;



function [H,oldaxlims,errstr] = arrow_clicks(H,ud,x,y,z,ax,oldaxlims)
% handle choosing arrow Start and/or Stop locations if necessary
    errstr = '';
    if isempty(H)|isempty(ud)|isempty(x), return; end;
    % determine which (if any) need Start and/or Stop
    needStart = all(isnan(ud(:,1:3)'))';
    needStop  = all(isnan(ud(:,4:6)'))';
    mask = any(needStart|needStop);
    if ~any(mask), return; end;
    ud(~mask,:)=[]; ax(:,~mask)=[];
    x(:,~mask)=[]; y(:,~mask)=[]; z(:,~mask)=[];
    % make them invisible for the time being
    set(H,'Visible','off');
    % save the current axes and limits modes; set to manual for the time being
    oldAx  = gca;
    limModes=get(ax(:),{'XLimMode','YLimMode','ZLimMode'});
    set(ax(:),{'XLimMode','YLimMode','ZLimMode'},{'manual','manual','manual'});
    % loop over each arrow that requires attention
    jj = find(mask);
    for ii=1:length(jj),
        h = H(jj(ii));
        axes(ax(ii));
        % figure out correct call
        if needStart(ii), prop='Start'; else, prop='Stop'; end;
        [wasInterrupted,errstr] = arrow_click(needStart(ii)&needStop(ii),h,prop,ax(ii));
        % handle errors and control-C
        if wasInterrupted,
            delete(H(jj(ii:end)));
            H(jj(ii:end))=[];
            oldaxlims(jj(ii:end),:)=[];
            break;
        end;
    end;
    % restore the axes and limit modes
    axes(oldAx);
    set(ax(:),{'XLimMode','YLimMode','ZLimMode'},limModes);

function [wasInterrupted,errstr] = arrow_click(lockStart,H,prop,ax)
% handle the clicks for one arrow
    fig = get(ax,'Parent');
    % save some things
    oldFigProps = {'Pointer','WindowButtonMotionFcn','WindowButtonUpFcn'};
    oldFigValue = get(fig,oldFigProps);
    oldArrowProps = {'EraseMode'};
    oldArrowValue = get(H,oldArrowProps);
    set(H,'EraseMode','background'); %because 'xor' makes shaft invisible unless Width>1
    global ARROW_CLICK_H ARROW_CLICK_PROP ARROW_CLICK_AX ARROW_CLICK_USE_Z
    ARROW_CLICK_H=H; ARROW_CLICK_PROP=prop; ARROW_CLICK_AX=ax;
    ARROW_CLICK_USE_Z=~arrow_is2DXY(ax)|~arrow_planarkids(ax);
    set(fig,'Pointer','crosshair');
    % set up the WindowButtonMotion so we can see the arrow while moving around
    set(fig,'WindowButtonUpFcn','set(gcf,''WindowButtonUpFcn'','''')', ...
            'WindowButtonMotionFcn','');
    if ~lockStart,
        set(H,'Visible','on');
        set(fig,'WindowButtonMotionFcn',[mfilename '(''callback'',''motion'');']);
    end;
    % wait for the button to be pressed
    [wasKeyPress,wasInterrupted,errstr] = arrow_wfbdown(fig);
    % if we wanted to click-drag, set the Start point
    if lockStart & ~wasInterrupted,
        pt = arrow_point(ARROW_CLICK_AX,ARROW_CLICK_USE_Z);
        feval(mfilename,H,'Start',pt,'Stop',pt);
        set(H,'Visible','on');
        ARROW_CLICK_PROP='Stop';
        set(fig,'WindowButtonMotionFcn',[mfilename '(''callback'',''motion'');']);
        % wait for the mouse button to be released
        eval('waitfor(fig,''WindowButtonUpFcn'','''');','wasInterrupted=1;');
        if wasInterrupted, errstr=lasterr; end;
    end;
    if ~wasInterrupted, feval(mfilename,'callback','motion'); end;
    % restore some things
    set(gcf,oldFigProps,oldFigValue);
    set(H,oldArrowProps,oldArrowValue);

function arrow_callback(varargin)
% handle redrawing callbacks
    if nargin==0, return; end;
    str = varargin{1};
    if ~ischar(str), error([upper(mfilename) ' got an invalid Callback command.']); end;
    s = lower(str);
    if strcmp(s,'motion'),
        % motion callback
        global ARROW_CLICK_H ARROW_CLICK_PROP ARROW_CLICK_AX ARROW_CLICK_USE_Z
        feval(mfilename,ARROW_CLICK_H,ARROW_CLICK_PROP,arrow_point(ARROW_CLICK_AX,ARROW_CLICK_USE_Z));
        drawnow;
    else,
        error([upper(mfilename) ' does not recognize ''' str(:).' ''' as a valid Callback option.']);
    end;

function out = arrow_point(ax,use_z)
% return the point on the given axes
    if nargin==0, ax=gca; end;
    if nargin<2, use_z=~arrow_is2DXY(ax)|~arrow_planarkids(ax); end;
    out = get(ax,'CurrentPoint');
    out = out(1,:);
    if ~use_z, out=out(1:2); end;

function [wasKeyPress,wasInterrupted,errstr] = arrow_wfbdown(fig)
% wait for button down ignoring object ButtonDownFcn's
    if nargin==0, fig=gcf; end;
    errstr = '';
    % save ButtonDownFcn values
    objs = findobj(fig);
    buttonDownFcns = get(objs,'ButtonDownFcn');
    mask=~strcmp(buttonDownFcns,''); objs=objs(mask); buttonDownFcns=buttonDownFcns(mask);
    set(objs,'ButtonDownFcn','');
    % save other figure values
    figProps = {'KeyPressFcn','WindowButtonDownFcn'};
    figValue = get(fig,figProps);
    % do the real work
    set(fig,'KeyPressFcn','set(gcf,''KeyPressFcn'','''',''WindowButtonDownFcn'','''');', ...
            'WindowButtonDownFcn','set(gcf,''WindowButtonDownFcn'','''')');
    lasterr('');
    wasInterrupted=0; eval('waitfor(fig,''WindowButtonDownFcn'','''');','wasInterrupted=1;');
    wasKeyPress = ~wasInterrupted & strcmp(get(fig,'KeyPressFcn'),'');
    if wasInterrupted, errstr=lasterr; end;
    % restore ButtonDownFcn and other figure values
    set(objs,'ButtonDownFcn',buttonDownFcns);
    set(fig,figProps,figValue);



function [out,is2D] = arrow_is2DXY(ax)
% check if axes are 2-D X-Y plots
    % may not work for modified camera angles, etc.
    out = logical(zeros(size(ax))); % 2-D X-Y plots
    is2D = out;                     % any 2-D plots
    views = get(ax(:),{'View'});
    views = cat(1,views{:});
    out(:) = abs(views(:,2))==90;
    is2D(:) = out(:) | all(rem(views',90)==0)';

function out = arrow_planarkids(ax)
% check if axes descendents all have empty ZData (lines,patches,surfaces)
    out = logical(ones(size(ax)));
    allkids = get(ax(:),{'Children'});
    for k=1:length(allkids),
        kids = get([findobj(allkids{k},'flat','Type','line')
                    findobj(allkids{k},'flat','Type','patch')
                    findobj(allkids{k},'flat','Type','surface')],{'ZData'});
        for j=1:length(kids),
            if ~isempty(kids{j}), out(k)=logical(0); break; end;
        end;
    end;



function arrow_fixlimits(axlimits)
% reset the axis limits as necessary
    if isempty(axlimits), disp([upper(mfilename) ' does not remember any axis limits to reset.']); end;
    for k=1:size(axlimits,1),
        if any(get(axlimits(k,1),'XLim')~=axlimits(k,2:3)), set(axlimits(k,1),'XLim',axlimits(k,2:3)); end;
        if any(get(axlimits(k,1),'YLim')~=axlimits(k,4:5)), set(axlimits(k,1),'YLim',axlimits(k,4:5)); end;
        if any(get(axlimits(k,1),'ZLim')~=axlimits(k,6:7)), set(axlimits(k,1),'ZLim',axlimits(k,6:7)); end;
    end;



function out = arrow_WarpToFill(notstretched,manualcamera,curax)
% check if we are in "WarpToFill" mode.
    out = strcmp(get(curax,'WarpToFill'),'on');
    % 'WarpToFill' is undocumented, so may need to replace this by
    % out = ~( any(notstretched) & any(manualcamera) );



function out = arrow_warnlimits(axlimits,narrows)
% create a warning message if we've changed the axis limits
    msg = '';
    switch (size(axlimits,1))
        case 1, msg='';
        case 2, msg='on two axes ';
        otherwise, msg='on several axes ';
    end;
    msg = [upper(mfilename) ' changed the axis limits ' msg ...
           'when adding the arrow'];
    if (narrows>1), msg=[msg 's']; end;
    out = [msg '.' sprintf('\n') '         Call ' upper(mfilename) ...
           ' FIXLIMITS to reset them now.'];



function arrow_copyprops(fm,to)
% copy line properties to patches
    props  = {'EraseMode','LineStyle','LineWidth','Marker','MarkerSize',...
              'MarkerEdgeColor','MarkerFaceColor','ButtonDownFcn',      ...
              'Clipping','DeleteFcn','BusyAction','HandleVisibility',   ...
              'Selected','SelectionHighlight','Visible'};
    lineprops  = {'Color',    props{:}};
    patchprops = {'EdgeColor',props{:}};
    patch2props = {'FaceColor',patchprops{:}};
    fmpatch = strcmp(get(fm,'Type'),'patch');
    topatch = strcmp(get(to,'Type'),'patch');
    set(to( fmpatch& topatch),patch2props,get(fm( fmpatch& topatch),patch2props)); %p->p
    set(to(~fmpatch&~topatch),lineprops,  get(fm(~fmpatch&~topatch),lineprops  )); %l->l
    set(to( fmpatch&~topatch),lineprops,  get(fm( fmpatch&~topatch),patchprops )); %p->l
    set(to(~fmpatch& topatch),patchprops, get(fm(~fmpatch& topatch),lineprops)  ,'FaceColor','none'); %l->p



function arrow_props
% display further help info about ARROW properties
    c = sprintf('\n');
    disp([c ...
    'ARROW Properties:  Default values are given in [square brackets], and other' c ...
    '                   acceptable equivalent property names are in (parenthesis).' c c ...
    '  Start           The starting points. For N arrows,            B' c ...
    '                  this should be a Nx2 or Nx3 matrix.          /|\           ^' c ...
    '  Stop            The end points. For N arrows, this          /|||\          |' c ...
    '                  should be a Nx2 or Nx3 matrix.             //|||\\        L|' c ...
    '  Length          Length of the arrowhead (in pixels on     ///|||\\\       e|' c ...
    '                  screen, points on a page). [16] (Len)    ////|||\\\\      n|' c ...
    '  BaseAngle       Angle (degrees) of the base angle       /////|D|\\\\\     g|' c ...
    '                  ADE.  For a simple stick arrow, use    ////  |||  \\\\    t|' c ...
    '                  BaseAngle=TipAngle. [90] (Base)       ///    |||    \\\   h|' c ...
    '  TipAngle        Angle (degrees) of tip angle ABC.    //<----->||      \\   |' c ...
    '                  [16] (Tip)                          /   base |||        \  V' c ...
    '  Width           Width of the base in pixels.  Not  E   angle ||<-------->C' c ...
    '                  the ''LineWidth'' prop. [0] (Wid)            |||tipangle' c ...
    '  Page            If provided, non-empty, and not NaN,         |||' c ...
    '                  this causes ARROW to use hardcopy            |||' c ...
    '                  rather than onscreen proportions.             A' c ...
    '                  This is important if screen aspect        -->   <-- width' c ...
    '                  ratio and hardcopy aspect ratio are    ----CrossDir---->' c ...
    '                  vastly different. []' c...
    '  CrossDir        A vector giving the direction towards which the fletches' c ...
    '                  on the arrow should go.  [computed such that it is perpen-' c ...
    '                  dicular to both the arrow direction and the view direction' c ...
    '                  (i.e., as if it was pasted on a normal 2-D graph)]  (Note' c ...
    '                  that CrossDir is a vector.  Also note that if an axis is' c ...
    '                  plotted on a log scale, then the corresponding component' c ...
    '                  of CrossDir must also be set appropriately, i.e., to 1 for' c ...
    '                  no change in that direction, >1 for a positive change, >0' c ...
    '                  and <1 for negative change.)' c ...
    '  NormalDir       A vector normal to the fletch direction (CrossDir is then' c ...
    '                  computed by the vector cross product [Line]x[NormalDir]). []' c ...
    '                  (Note that NormalDir is a vector.  Unlike CrossDir,' c ...
    '                  NormalDir is used as is regardless of log-scaled axes.)' c ...
    '  Ends            Set which end has an arrowhead.  Valid values are ''none'',' c ...
    '                  ''stop'', ''start'', and ''both''. [''stop''] (End)' c...
    '  ObjectHandles   Vector of handles to previously-created arrows to be' c ...
    '                  updated or line objects to be converted to arrows.' c ...
    '                  [] (Object,Handle)' c ]);



function out = arrow_demo
    % demo
    % create the data
    [x,y,z] = peaks;
    [ddd,out.iii]=max(z(:));
    out.axlim = [min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z(:)) max(z(:))];
    
    % modify it by inserting some NaN's
    [m,n] = size(z);
    m = floor(m/2);
    n = floor(n/2);
    z(1:m,1:n) = nan(m,n);
    
    % graph it
    clf('reset');
    out.hs=surf(x,y,z);
    out.x=x; out.y=y; out.z=z;
    xlabel('x'); ylabel('y');
            
function h = arrow_demo3(in)
    % set the view
    axlim = in.axlim;
    axis(axlim);
    zlabel('z');
    %set(in.hs,'FaceColor','interp');
    view(viewmtx(-37.5,30,20));
    title(['Demo of the capabilities of the ARROW function in 3-D']);
    
    % Normal blue arrow
    h1 = feval(mfilename,[axlim(1) axlim(4) 4],[-.8 1.2 4], ...
               'EdgeColor','b','FaceColor','b');
    
    % Normal white arrow, clipped by the surface
    h2 = feval(mfilename,axlim([1 4 6]),[0 2 4]);
    t=text(-2.4,2.7,7.7,'arrow clipped by surf');
    
    % Baseangle<90
    h3 = feval(mfilename,[3 .125 3.5],[1.375 0.125 3.5],30,50);
    t2=text(3.1,.125,3.5,'local maximum');
    
    % Baseangle<90, fill and edge colors different
    h4 = feval(mfilename,axlim(1:2:5)*.5,[0 0 0],36,60,25, ...
               'EdgeColor','b','FaceColor','c');
    t3=text(axlim(1)*.5,axlim(3)*.5,axlim(5)*.5-.75,'origin');
    set(t3,'HorizontalAlignment','center');
    
    % Baseangle>90, black fill
    h5 = feval(mfilename,[-2.9 2.9 3],[-1.3 .4 3.2],30,120,[],6, ...
               'EdgeColor','r','FaceColor','k','LineWidth',2);
    
    % Baseangle>90, no fill
    h6 = feval(mfilename,[-2.9 2.9 1.3],[-1.3 .4 1.5],30,120,[],6, ...
               'EdgeColor','r','FaceColor','none','LineWidth',2);
    
    % Stick arrow
    h7 = feval(mfilename,[-1.6 -1.65 -6.5],[0 -1.65 -6.5],[],16,16);
    t4=text(-1.5,-1.65,-7.25,'global mininum');
    set(t4,'HorizontalAlignment','center');
    
    % Normal, black fill
    h8 = feval(mfilename,[-1.4 0 -7.2],[-1.4 0 -3],'FaceColor','k');
    t5=text(-1.5,0,-7.75,'local minimum');
    set(t5,'HorizontalAlignment','center');
    
    % Gray fill, crossdir specified, 'LineStyle' --
    h9 = feval(mfilename,[-3 2.2 -6],[-3 2.2 -.05],36,[],27,6,[],[0 -1 0], ...
               'EdgeColor','k','FaceColor',.75*[1 1 1],'LineStyle','--');
    
    % a series of normal arrows, linearly spaced, crossdir specified
    h10y=(0:4)'/3;
    h10 = feval(mfilename,[-3*ones(size(h10y)) h10y -6.5*ones(size(h10y))], ...
                [-3*ones(size(h10y)) h10y -.05*ones(size(h10y))], ...
                12,[],[],[],[],[0 -1 0]);
    
    % a series of normal arrows, linearly spaced
    h11x=(1:.33:2.8)';
    h11 = feval(mfilename,[h11x -3*ones(size(h11x)) 6.5*ones(size(h11x))], ...
                [h11x -3*ones(size(h11x)) -.05*ones(size(h11x))]);
    
    % series of magenta arrows, radially oriented, crossdir specified
    h12x=2; h12y=-3; h12z=axlim(5)/2; h12xr=1; h12zr=h12z; ir=.15;or=.81;
    h12t=(0:11)'/6*pi;
    h12 = feval(mfilename,                                           ...
                [h12x+h12xr*cos(h12t)*ir h12y*ones(size(h12t))       ...
                 h12z+h12zr*sin(h12t)*ir],[h12x+h12xr*cos(h12t)*or   ...
                 h12y*ones(size(h12t)) h12z+h12zr*sin(h12t)*or],     ...
                10,[],[],[],[],                                      ...
                [-h12xr*sin(h12t) zeros(size(h12t)) h12zr*cos(h12t)],...
                'FaceColor','none','EdgeColor','m');
    
    % series of normal arrows, tangentially oriented, crossdir specified
    or13=.91; h13t=(0:.5:12)'/6*pi;
    locs = [h12x+h12xr*cos(h13t)*or13 h12y*ones(size(h13t)) h12z+h12zr*sin(h13t)*or13];
    h13 = feval(mfilename,locs(1:end-1,:),locs(2:end,:),6);
    
    % arrow with no line ==> oriented downwards
    h14 = feval(mfilename,[3 3 .100001],[3 3 .1],30);
    t6=text(3,3,3.6,'no line'); set(t6,'HorizontalAlignment','center');
    
    % arrow with arrowheads at both ends
    h15 = feval(mfilename,[-.5 -3 -3],[1 -3 -3],'Ends','both','FaceColor','g', ...
                'Length',20,'Width',3,'CrossDir',[0 0 1],'TipAngle',25);
    
    h=[h1;h2;h3;h4;h5;h6;h7;h8;h9;h10;h11;h12;h13;h14;h15];

function h = arrow_demo2(in)
    axlim = in.axlim;
    dolog = 1;
    if (dolog), set(in.hs,'YData',10.^get(in.hs,'YData')); end;
    shading('interp');
    view(2);
    title(['Demo of the capabilities of the ARROW function in 2-D']);
    hold on; [C,H]=contour(in.x,in.y,in.z,20,'-'); hold off;
    for k=H',
        set(k,'ZData',(axlim(6)+1)*ones(size(get(k,'XData'))),'Color','k');
        if (dolog), set(k,'YData',10.^get(k,'YData')); end;
    end;
    if (dolog), axis([axlim(1:2) 10.^axlim(3:4)]); set(gca,'YScale','log');
    else,       axis(axlim(1:4)); end;
    
    % Normal blue arrow
    start = [axlim(1) axlim(4) axlim(6)+2];
    stop  = [in.x(in.iii) in.y(in.iii) axlim(6)+2];
    if (dolog), start(:,2)=10.^start(:,2); stop(:,2)=10.^stop(:,2); end;
    h1 = feval(mfilename,start,stop,'EdgeColor','b','FaceColor','b');
    
    % three arrows with varying fill, width, and baseangle
    start = [-3   -3   10; -3   -1.5 10; -1.5 -3   10];
    stop  = [-.03 -.03 10; -.03 -1.5 10; -1.5 -.03 10];
    if (dolog), start(:,2)=10.^start(:,2); stop(:,2)=10.^stop(:,2); end;
    h2 = feval(mfilename,start,stop,24,[90;60;120],[],[0;0;4],'Ends',str2mat('both','stop','stop'));
    set(h2(2),'EdgeColor',[0 .35 0],'FaceColor',[0 .85 .85]);
    set(h2(3),'EdgeColor','r','FaceColor',[1 .5 1]);
    h=[h1;h2];

function out = trueornan(x)
if isempty(x),
    out=x;
else,
    out = isnan(x);
    out(~out) = x(~out);
end;

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            