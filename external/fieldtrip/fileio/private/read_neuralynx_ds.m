function [dat] = read_neuralynx_ds(dirname, hdr, begsample, endsample, chanindx)

% READ_NEURALYNX_DS reads multiple single-channel Neuralynx files that are
% all contained in a single directory. Each file is treated as a single
% channel of a combined multi-channel dataset.
%
% Use as
%   [hdr] = read_neuralynx_ds(dirname)
%   [dat] = read_neuralynx_ds(dirname, hdr, begsample, endsample, chanindx)
%
% A Neuralynx dataset consists of a directory containing separate files,
% one for each channel. All Neuralynx datafiles starts with a 16k header
% (in ascii format), followed by an arbitrary number of data records. The
% format of the data records depend on the type of data contained in the
% channel (e.g. continuous or spike data).
%
% To read the timestamps of spike waveforms (nse) or clustered spikes (nts),
% the header should contain the fields
%   hdr.FirstTimeStamp
%   hdr.TimeStampPerSample
% These can only be obtained from the corresponding simultaneous LFP
% and/or MUA recordings.
%
% See also READ_NEURALYNX_NCS, READ_NEURALYNX_NSE, READ_NEURALYNX_NTS

% Copyright (C) 2006-2007, Robert Oostenveld
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

needhdr = (nargin==1);
needdat = (nargin>=2);

if needhdr
  % get the list of filenames
  ls = dir(dirname);
  ls = ls(~cell2mat({ls.isdir}));
  fname = {};
  for i=1:length(ls)
    fname{i} = fullfile(dirname, ls(i).name);
  end

  ftype = zeros(length(fname), 1);
  for i=1:length(fname)
    if     ft_filetype(fname{i}, 'neuralynx_ncs')
      ftype(i) = 1;
    elseif ft_filetype(fname{i}, 'neuralynx_nse')
      ftype(i) = 2;
    elseif ft_filetype(fname{i}, 'neuralynx_nts')
      ftype(i) = 3;
    end
  end

  % only remember the filenames that are relevant
  fname = fname(ftype>0);
  ftype = ftype(ftype>0);
  ftype_ncs = find(ftype==1);
  ftype_nse = find(ftype==2);
  ftype_nts = find(ftype==3);

  if length(fname)==0
    ft_error('the dataset directory contains no supported files');
  end

  for i=1:length(fname)
    % this will only work if all files within a dataset return a similar header structure
    switch ftype(i)
      case 1
        orig(i) = read_neuralynx_ncs(fname{i}, 0, 0);
      case 2
        orig(i) = read_neuralynx_nse(fname{i}, 0, 0);
      case 3
        orig(i) = read_neuralynx_nts(fname{i}, 0, 0);
      otherwise
        ft_error('unsupported');
    end
  end

  % combine the information from the different files in a single header
  for i=1:length(orig)
    if isfield(orig(i).hdr, 'NLX_Base_Class_Name')
      label{i}           = orig(i).hdr.NLX_Base_Class_Name;
    else
      label{i}           = orig(i).hdr.AcqEntName;
    end
    if isfield(orig(i).hdr, 'SubSamplingInterleave')
      SamplingFrequency(i) = orig(i).hdr.SamplingFrequency / orig(i).hdr.SubSamplingInterleave;
    else
      SamplingFrequency(i) = orig(i).hdr.SamplingFrequency;
    end
    ADBitVolts(i)        = orig(i).hdr.ADBitVolts;
    FirstTimeStamp(i)    = orig(i).hdr.FirstTimeStamp;
    LastTimeStamp(i)     = orig(i).hdr.LastTimeStamp;
    NRecords(i)          = orig(i).NRecords;
    % Note that the last timestamp corresponds with the first sample of the last
    % record and not with the last sample in the file.
  end

  for i=1:length(orig)
    % timestamps are measured in units of approximately 1us
    % in case of 32556 Hz sampling rate, there are approximately 30.7 timestamps per sample
    % in case of 1000 Hz sampling rate,  there are approximately 1000 timestamps per sample
    % note that the last timestamp in the original header corresponds with the
    % first sample of the last record, and not with the last sample
    switch ftype(i)
      case 1
        % ensure that the last timestamp refers to the last sample
        recordsize = 512; % each record contains 512 samples
      case 2
        recordsize = 32; % each record contains 32 samples
      case 3
        ft_error('this has not been implemented yet');
      otherwise
        ft_error('unsupported');
    end
    % ensure that the last timestamp refers to the last sample
    TimeStampPerSample(i) = double(LastTimeStamp(i)-FirstTimeStamp(i))/((NRecords(i)-1)*recordsize);  % this should be in double precision, since it can be fractional
    LastTimeStamp(i)      = LastTimeStamp(i) + uint64((recordsize-1)*TimeStampPerSample(i));          % this should be in uint64 precision
  end % for length(orig)

  if any(SamplingFrequency~=SamplingFrequency(1))
    ft_warning('not all channels have the same sampling rate');
  end

  if ~isempty(ftype_ncs)
    if any(FirstTimeStamp(ftype_ncs)~=FirstTimeStamp(ftype_ncs(1)))
      % there seems to be a matlab bug (observed in Matlab75 on windows) which causes this uint64 comparison to fail if there are exactly 8 files
      % therefore check once more after converting them to double
      if any(double(FirstTimeStamp(ftype_ncs))~=double(FirstTimeStamp(ftype_ncs(1))))
        ft_error('not all continuous channels start at the same time');
      end
    end
    if any(LastTimeStamp(ftype_ncs)~=LastTimeStamp(ftype_ncs(1)))
      % there seems to be a matlab bug (observed in Matlab75 on windows) which causes this uint64 comparison to fail if there are exactly 8 files
      % therefore check once more after converting them to double
      if any(double(LastTimeStamp(ftype_ncs))~=double(LastTimeStamp(ftype_ncs(1))))
        ft_warning('not all continuous channels end at the same time');
      end
    end
    if any(NRecords(ftype_ncs)~=NRecords(ftype_ncs(1)))
      ft_warning('not all continuous channels have the same number of records');
    end
  end % if ftype_ncs

  % construct the header that applies to all channels combined
  hdr.nChans         = length(label);
  hdr.label          = label;
  hdr.filename       = fname;
  hdr.nTrials        = 1;                           % it is continuous
  hdr.Fs             = SamplingFrequency(1);
  hdr.nSamplesPre    = 0;                           % it is continuous

  if ~isempty(ftype_ncs)
    % these elements are only relevant for continuously sampled channels
    hdr.nSamples           = NRecords(1) * 512;
    hdr.FirstTimeStamp     = FirstTimeStamp(1);
    hdr.LastTimeStamp      = LastTimeStamp(1);
    hdr.TimeStampPerSample = TimeStampPerSample(1);
  end

  % remember the original header details
  hdr.orig = orig(:);

  % return the header
  dat = hdr;

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the data of the selected channels (i.e. files)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if nargin<5
    % select all channels
    chanindx = 1:length(hdr.label);
  end
  nchan   = length(chanindx);
  nsample = endsample-begsample+1;
  dat     = zeros(nchan, nsample);

  for i=1:nchan
    thischan = chanindx(i);
    thisfile = hdr.filename{thischan};
    switch ft_filetype(thisfile)
    case 'neuralynx_ncs'
      % determine the records that contain the sample numbers of the requested segment
      begrecord  = ceil(begsample/512);
      endrecord  = ceil(endsample/512);
      ncs = read_neuralynx_ncs(thisfile, begrecord, endrecord);
      % copy the selected samples into the output
      begsel = begsample - (begrecord-1)*512;
      endsel = endsample - (begrecord-1)*512;
      dat(i,:) = ncs.dat(begsel:endsel);

    case 'neuralynx_nse'
      % read all spike waveforms and timestamps
      nse = read_neuralynx_nse(thisfile);
      % convert the timestamps into samples
      fprintf('%d timstamps\n', length(nse.TimeStamp));
      sample = double(nse.TimeStamp-hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
      sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
      dat(i,sample) = dat(i,sample) + 1;

    case 'neuralynx_nts'
      % read all timestamps
      nts = read_neuralynx_nts(thisfile);
      % convert the timestamps into samples
      sample = double(nse.TimeStamp-hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
      sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
      dat(i,sample) = dat(i,sample) + 1;

    end % switch ft_filetype
  end % for nchan
end % reading data
                                                                                                                                                                                                                                                                                                                                                       for only one channel
  fid = fopen(filename, 'rb', 'ieee-le');
  datoffset = (begsample-1)*blocksize*4;
  fseek(fid, hdroffset + datoffset, 'bof');  % skip the header and samples that do not have to be read
  fseek(fid, (chanindx-1)*4, 'cof');         % skip to the channel of interest
  % Note that the last block with 274*4 bytes can sometimes be incomplete, which is relevant when endsample=inf
  dat = fread(fid, [1 (endsample-begsample+1)], 'int32=>int32', (nboards*32+18-1)*4);
  if size(dat,2)<(endsample-begsample+1)
    ft_error('could not read all samples');
  end
  fclose(fid);
end

                                                                                                                                                                                                                                                                                                                                                                                                                                             ock is the pointer
% to the next block. Total number of entries is in obj.Qi.nrEntries

MainIndex = struct();
curIdx = 0;
nextIndexPointer = indexIdx;
curIdx2 = 1;
while curIdx < nrEntries
    
    fseek(h, nextIndexPointer, 'bof');
    nrIdx = fread(h,1, 'uint64');
    MainIndex(curIdx + nrIdx).sectionIdx = 0;   % Preallocate next set of indices
    var = fread(h,3*nrIdx, 'uint64');
    for i = 1: nrIdx
        MainIndex(curIdx + i).sectionIdx = var(3*(i-1)+1);
        MainIndex(curIdx + i).offset = var(3*(i-1)+2);
        MainIndex(curIdx + i).blockL = mod(var(3*(i-1)+3),2^32);
        MainIndex(curIdx + i).sectionL = round(var(3*(i-1)+3)/2^32);
    end
    nextIndexPointer = fread(h,1, 'uint64');
    curIdx = curIdx + i;
    curIdx2=curIdx2+1;
end
end

function infoGuids = read_nervus_header_infoGuids(h, StaticPackets, MainIndex)


infoIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'InfoGuids'),1)).index;
indexInstance = MainIndex(find([MainIndex.sectionIdx]==infoIdx,1));
nrInfoGuids = indexInstance.sectionL/16;
infoGuids = struct();
fseek(h, indexInstance.offset,'bof');
for i = 1:nrInfoGuids
    guidmixed = fread(h,16, 'uint8')';
    guidnonmixed = [guidmixed(04), guidmixed(03), guidmixed(02), guidmixed(01), ...
        guidmixed(06), guidmixed(05), guidmixed(08), guidmixed(07), ...
        guidmixed(09), guidmixed(10), guidmixed(11), guidmixed(12), ...
        guidmixed(13), guidmixed(15), guidmixed(15), guidmixed(16)];
    infoGuids(i).guid = num2str(guidnonmixed,'%02X');
end
end


function dynamicPackets = read_nervus_header_dynamicpackets(h, StaticPackets, MainIndex)
dynamicPackets = struct();
indexIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'InfoChangeStream'),1)).index;
offset = MainIndex(indexIdx).offset;
nrDynamicPackets = MainIndex(indexIdx).sectionL / 48;
fseek(h, offset, 'bof');

%Read first only the dynamic packets structure without actual data
for i = 1: nrDynamicPackets
    dynamicPackets(i).offset = offset+i*48;
    guidmixed = fread(h,16, 'uint8')';
    guidnonmixed = [guidmixed(04), guidmixed(03), guidmixed(02), guidmixed(01), ...
        guidmixed(06), guidmixed(05), guidmixed(08), guidmixed(07), ...
        guidmixed(09), guidmixed(10), guidmixed(11), guidmixed(12), ...
        guidmixed(13), guidmixed(14), guidmixed(15), guidmixed(16)];
    dynamicPackets(i).guid = num2str(guidnonmixed, '%02X');
    dynamicPackets(i).guidAsStr = sprintf('{%02X%02X%02X%02X-%02X%02X-%02X%02X-%02X%02X-%02X%02X%02X%02X%02X%02X}', guidnonmixed);
    dynamicPackets(i).date = datenum(1899,12,31) + fread(h,1,'double');
    dynamicPackets(i).datefrac = fread(h,1,'double');
    dynamicPackets(i).internalOffsetStart = fread(h,1, 'uint64')';
    dynamicPackets(i).packetSize = fread(h,1, 'uint64')';
    dynamicPackets(i).data = zeros(0, 1,'uint8');
    
    switch dynamicPackets(i).guid
        case 'BF7C95EF6C3B4E709E11779BFFF58EA7'
            dynamicPackets(i).IDStr = 'CHANNELGUID';
        case '8A19AA48BEA040D5B89F667FC578D635'
            dynamicPackets(i).IDStr = 'DERIVATIONGUID';
        case 'F824D60C995E4D949578893C755ECB99'
            dynamicPackets(i).IDStr = 'FILTERGUID';
        case '0295036135BB4A229F0BC78AAA5DB094'
            dynamicPackets(i).IDStr = 'DISPLAYGUID';
        case '782B34E88E514BB997013227BB882A23'
            dynamicPackets(i).IDStr = 'ACCINFOGUID';
        case 'A271CCCB515D4590B6A1DC170C8D6EE2'
            dynamicPackets(i).IDStr = 'TSGUID';
        case 'D01B34A09DBD11D393D300500400C148'
            dynamicPackets(i).IDStr = 'AUDIOINFOGUID';
        otherwise
            dynamicPackets(i).IDStr = 'UNKNOWN';
    end
end

%Then read the actual data from the pointers above
for i = 1: nrDynamicPackets
    %Look up the GUID of this dynamic packet in the static packets
    % to find the section index
    
    infoIdx = StaticPackets(find(strcmp({StaticPackets.tag},dynamicPackets(i).guidAsStr),1)).index;
    
    %Matching index segments
    indexInstances = MainIndex([MainIndex.sectionIdx] == infoIdx);
    
    %Then, treat all these sections as one contiguous memory block
    % and grab this packet across these instances
    
    internalOffset = 0;
    remainingDataToRead = dynamicPackets(i).packetSize;
    %disp(['Target packet ' dynamicPackets(i).IDStr ' : ' num2str(dynamicPackets(i).internalOffsetStart) ' to ' num2str(dynamicPackets(i).internalOffsetStart+dynamicPackets(i).packetSize) ' target read length ' num2str(remainingDataToRead)]);
    currentTargetStart = dynamicPackets(i).internalOffsetStart;
    for j = 1: size(indexInstances,2)
        currentInstance = indexInstances(j);
        
        %hitInThisSegment = '';
        if (internalOffset <= currentTargetStart) && (internalOffset+currentInstance.sectionL) >= currentTargetStart
            
            startAt = currentTargetStart;
            stopAt =  min(startAt+remainingDataToRead, internalOffset+currentInstance.sectionL);
            readLength = stopAt-startAt;
            
            filePosStart = currentInstance.offset+startAt-internalOffset;
            fseek(h,filePosStart, 'bof');
            dataPart = fread(h,readLength,'uint8=>uint8');
            dynamicPackets(i).data = cat(1, dynamicPackets(i).data, dataPart);
            
            %hitInThisSegment = ['HIT at  ' num2str(startAt) ' to ' num2str(stopAt)];
            %if (readLength < remainingDataToRead)
            %    hitInThisSegment = [hitInThisSegment ' (partial ' num2str(readLength) ' )'];
            %else
            %    hitInThisSegment = [hitInThisSegment ' (finished - this segment contributed ' num2str(readLength) ' )'];
            %end
            %hitInThisSegment = [hitInThisSegment ' abs file pos ' num2str(filePosStart) ' - ' num2str(filePosStart+readLength)];
            
            remainingDataToRead = remainingDataToRead-readLength;
            currentTargetStart = currentTargetStart + readLength;
            
        end
        %disp(['    Index ' num2str(j) ' Offset: ' num2str(internalOffset) ' to ' num2str(internalOffset+currentInstance.sectionL) ' ' num2str(hitInThisSegment)]);
        
        internalOffset = internalOffset + currentInstance.sectionL;
    end
end
end

function PatientInfo = read_nervus_header_patient(h, StaticPackets, Index)
%% Get PatientGUID
PatientInfo = struct();

infoProps = { 'patientID', 'firstName','middleName','lastName',...
    'altID','mothersMaidenName','DOB','DOD','street','sexID','phone',...
    'notes','dominance','siteID','suffix','prefix','degree','apartment',...
    'city','state','country','language','height','weight','race','religion',...
    'maritalStatus'};

infoIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'PATIENTINFOGUID'),1)).index;
indexInstance = Index(find([Index.sectionIdx]==infoIdx,1));
fseek(h, indexInstance.offset,'bof');
guid = fread(h, 16, 'uint8');
lSection = fread(h, 1, 'uint64');
% reserved = fread(h, 3, 'uint16');
nrValues = fread(h,1,'uint64');
nrBstr = fread(h,1,'uint64');

for i = 1:nrValues
    id = fread(h,1,'uint64');
    switch id
        case {7,8}
            unix_time = (fread(h,1, 'double')*(3600*24)) - 2209161600; % 2208988800; %8
            obj.segments(i).dateStr = datestr(unix_time/86400 + datenum(1970,1,1));
            value = datevec( obj.segments(i).dateStr );
            value = value([3 2 1]);
        case {23,24}
            value = fread(h,1,'double');
        otherwise
            value = 0;
    end
    PatientInfo.(infoProps{id}) = value;
end

strSetup = fread(h,nrBstr*2,'uint64');

for i=1:2:(nrBstr*2)
    id  = strSetup(i);
    value = deblank(cast(fread(h, strSetup(i+1) + 1, 'uint16'),'char')');
    info.(infoProps{id}) = value;
end

end

function sigInfo = read_nervus_header_SignalInfo(h, StaticPackets, Index, ITEMNAMESIZE, LABELSIZE, UNITSIZE)
infoIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'InfoGuids'),1)).index;
indexInstance = Index(find([Index.sectionIdx]==infoIdx,1));
fseek(h, indexInstance.offset,'bof');

sigInfo = struct();
SIG_struct = struct();
sensorIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'SIGNALINFOGUID'),1)).index;
indexInstance = Index(find([Index.sectionIdx]==sensorIdx,1));
fseek(h, indexInstance.offset,'bof');
SIG_struct.guid = fread(h, 16, 'uint8');
SIG_struct.name = fread(h, ITEMNAMESIZE, '*char');
unkown = fread(h, 152, '*char');         %#ok<NASGU>
fseek(h, 512, 'cof');
nrIdx = fread(h,1, 'uint16');  %783
misc1 = fread(h,3, 'uint16'); %#ok<NASGU>

for i = 1: nrIdx
    sigInfo(i).sensorName = deblank(cast(fread(h, LABELSIZE, 'uint16'),'char')');
    sigInfo(i).transducer = deblank(cast(fread(h, UNITSIZE, 'uint16'),'char')');
    sigInfo(i).guid = fread(h, 16, '*uint8');
    sigInfo(i).bBiPolar = logical(fread(h, 1 ,'uint32'));
    sigInfo(i).bAC = logical(fread(h, 1 ,'uint32'));
    sigInfo(i).bHighFilter = logical(fread(h, 1 ,'uint32'));
    sigInfo(i).color =  fread(h, 1 ,'uint32');
    reserved = fread(h, 256, '*char'); %#ok<NASGU>
end
end


function channelInfo = read_nervus_header_ChannelInfo(h, StaticPackets, Index, ITEMNAMESIZE, LABELSIZE)
%% Get CHANNELINFO (CHANNELGUID)
CH_struct = struct();
sensorIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'CHANNELGUID'),1)).index;
indexInstance = Index(find([Index.sectionIdx]==sensorIdx,1));
fseek(h, indexInstance.offset,'bof');
CH_struct.guid = fread(h, 16, 'uint8');
CH_struct.name = fread(h, ITEMNAMESIZE, '*char');
fseek(h, 152, 'cof');
CH_struct.reserved = fread(h, 16, 'uint8');
CH_struct.deviceID = fread(h, 16, 'uint8');
fseek(h, 488, 'cof');

nrIdx = fread(h,2, 'uint32');  %783
channelInfo = struct();
for i = 1: nrIdx(2)
    channelInfo(i).sensor = deblank(cast(fread(h, LABELSIZE, 'uint16'),'char')');
    channelInfo(i).samplingRate = fread(h,1,'double');
    channelInfo(i).bOn = logical(fread(h, 1 ,'uint32'));
    channelInfo(i).lInputID = fread(h, 1 ,'uint32');
    channelInfo(i).lInputSettingID = fread(h,1,'uint32');
    channelInfo(i).reserved = fread(h,4,'char');
    fseek(h, 128, 'cof');
end

curIdx = 0;
for i = 1: length(channelInfo)
    if channelInfo(i).bOn
        channelInfo(i).indexID = curIdx;
        curIdx = curIdx+1;
    else
        channelInfo(i).indexID = -1;
    end
end
end

function [TSInfo] = read_nervus_header_TSInfo(DynamicPackets, TSLABELSIZE, LABELSIZE)
tsPackets = DynamicPackets(strcmp({DynamicPackets.IDStr},'TSGUID'));

if isempty(tsPackets)
    ft_error(['No TSINFO found']);
end    

tsPacket = tsPackets(1);
TSInfo = read_nervus_header_one_TSInfo(tsPacket, TSLABELSIZE, LABELSIZE);

if length(tsPackets) > 1
    allEqual = 1;
    for i = 2: size(tsPackets,2)
        nextTsPacket = tsPackets(i);
        nextTSInfo = read_nervus_header_one_TSInfo(nextTsPacket, TSLABELSIZE, LABELSIZE);       
        areEqual = compareTsInfoPackets(TSInfo, nextTSInfo);
        if (areEqual == 0)
            allEqual = 0;
            break;
        end
    end    
    if (allEqual == 0)            
        ft_error('Multiple TSInfo packets found and they are not the same.');
    end
end
end

function [TSInfo] = read_nervus_header_one_TSInfo(tsPacket, TSLABELSIZE, LABELSIZE)    
    TSInfo = struct();
    elems = typecast(tsPacket.data(753:756),'uint32');
    %alloc = typecast(tsPacket.data(757:760),'uint32');
    
    offset = 761;
    for i = 1:elems
        internalOffset = 0;
        TSInfo(i).label = deblank(char(typecast(tsPacket.data(offset:(offset+TSLABELSIZE-1))','uint16')));
        internalOffset = internalOffset + TSLABELSIZE*2;
        TSInfo(i).activeSensor = deblank(char(typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+LABELSIZE))','uint16')));
        internalOffset = internalOffset + TSLABELSIZE;
        TSInfo(i).refSensor = deblank(char(typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','uint16')));
        internalOffset = internalOffset + 8;
        internalOffset = internalOffset + 56;
        TSInfo(i).lowcut = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        internalOffset = internalOffset + 8;
        TSInfo(i).hiCut = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        internalOffset = internalOffset + 8;
        TSInfo(i).samplingRate = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        internalOffset = internalOffset + 8;
        TSInfo(i).resolution = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        internalOffset = internalOffset + 8;
        TSInfo(i).specialMark = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+2))','uint16');
        internalOffset = internalOffset + 2;
        TSInfo(i).notch = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+2))','uint16');
        internalOffset = internalOffset + 2;
        TSInfo(i).eeg_offset = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        offset = offset + 552;
        %disp([num2str(i) ' : ' TSInfo(i).label ' : ' TSInfo(i).activeSensor ' : ' TSInfo(i).refSensor ' : ' num2str(TSInfo(i).samplingRate)]);
    end

end

function areEqual = compareTsInfoPackets(TSInfo1, TSInfo2)    
    areEqual = 1;
    if (size(TSInfo1,2) ~= size(TSInfo2,2))
        areEqual = 0;
    else
        for i = 1:size(TSInfo1,2)
            if (~strcmp(TSInfo1(i).label,TSInfo2(i).label))
                areEqual = 0;
                break;
            end
            if (~strcmp(TSInfo1(i).activeSensor,TSInfo2(i).activeSensor))
                areEqual = 0;
                break;
            end
            if (~strcmp(TSInfo1(i).refSensor,TSInfo2(i).refSensor))
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).lowcut ~= TSInfo2(i).lowcut)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).hiCut ~= TSInfo2(i).hiCut)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).samplingRate ~= TSInfo2(i).samplingRate)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).resolution ~= TSInfo2(i).resolution)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).specialMark ~= TSInfo2(i).specialMark)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).notch ~= TSInfo2(i).notch)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).eeg_offset ~= TSInfo2(i).eeg_offset)
                areEqual = 0;
                break;
            end
        end
    end
        
end

function [segments] = read_nervus_header_Segments(h, StaticPackets, Index, TSInfo)
%% Get Segment Start Times
segmentIdx = StaticPackets(find(strcmp({StaticPackets.IDStr}, 'SegmentStream'),1)).index;
indexIdx = find([Index.sectionIdx] == segmentIdx, 1);
segmentInstance = Index(indexIdx);

nrSegments = segmentInstance.sectionL/152;
fseek(h, segmentInstance.offset,'bof');
segments = struct();
for i = 1: nrSegments
    dateOLE = fread(h,1, 'double');
    segments(i).dateOLE = dateOLE;
    unix_time = (dateOLE*(3600*24)) - 2209161600; % 2208988800; %8
    segments(i).dateStr = datestr(unix_time/86400 + datenum(1970,1,1));
    datev = datevec( segments(i).dateStr );
    segments(i).startDate = datev(1:3);
    segments(i).startTime = datev(4:6);
    fseek(h, 8 , 'cof'); %16
    segments(i).duration = fread(h,1, 'double'); %24
    fseek(h, 128 , 'cof'); %152       
end

% Get nrValues per segment and channel
for iSeg = 1:length(segments)
    % Add Channel Names to segments
    segments(iSeg).chName = {TSInfo.label};
    segments(iSeg).refName = {TSInfo.refSensor};
    segments(iSeg).samplingRate = [TSInfo.samplingRate];
    segments(iSeg).scale = [TSInfo.resolution];
    segments(iSeg).sampleCount = max(segments(iSeg).samplingRate*segments(iSeg).duration);
end
end

function [eventMarkers] = read_nervus_header_events(h, StaticPackets, Index)
%% Get events  - Andrei Barborica, Dec 2015
% Find sequence of events, that are stored in the section tagged 'Events'
eventsSection = strcmp({StaticPackets.tag}, 'Events');
idxSection = find(eventsSection);
indexIdx = find([Index.sectionIdx] == StaticPackets(idxSection).index);
offset = Index(indexIdx).offset;

ePktLen = 272;    % Event packet length, see EVENTPACKET definition
eMrkLen = 240;    % Event marker length, see EVENTMARKER definition
evtPktGUID = hex2dec({'80', 'F6', '99', 'B7', 'A4', '72', 'D3', '11', '93', 'D3', '00', '50', '04', '00', 'C1', '48'}); % GUID for event packet header
HCEVENT_ANNOTATION        =  '{A5A95612-A7F8-11CF-831A-0800091B5BDA}';
HCEVENT_SEIZURE           =  '{A5A95646-A7F8-11CF-831A-0800091B5BDA}';
HCEVENT_FORMATCHANGE      =  '{08784382-C765-11D3-90CE-00104B6F4F70}';
HCEVENT_PHOTIC            =  '{6FF394DA-D1B8-46DA-B78F-866C67CF02AF}';
HCEVENT_POSTHYPERVENT     =  '{481DFC97-013C-4BC5-A203-871B0375A519}';
HCEVENT_REVIEWPROGRESS    =  '{725798BF-CD1C-4909-B793-6C7864C27AB7}';
HCEVENT_EXAMSTART         =  '{96315D79-5C24-4A65-B334-E31A95088D55}';
HCEVENT_HYPERVENTILATION  =  '{A5A95608-A7F8-11CF-831A-0800091B5BDA}';                            
HCEVENT_IMPEDANCE         =  '{A5A95617-A7F8-11CF-831A-0800091B5BDA}';
DAYSECS = 86400.0;  % From nrvdate.h

fseek(h,offset,'bof');
pktGUID = fread(h,16,'uint8');
pktLen  = fread(h,1,'uint64');
eventMarkers = struct();
i = 0;    % Event counter
while (pktGUID == evtPktGUID)
    i = i + 1;
    % Please refer to EVENTMARKER structure in the Nervus file documentation
    fseek(h,8,'cof'); % Skip eventID, not used
    evtDate           = fread(h,1,'double');
    evtDateFraction   = fread(h,1,'double');
    eventMarkers(i).dateOLE = evtDate;
    eventMarkers(i).dateFraction = evtDateFraction;
    evtPOSIXTime = evtDate*DAYSECS + evtDateFraction - 2209161600; % 2208988800; %8
    eventMarkers(i).dateStr = datestr(evtPOSIXTime/DAYSECS + datenum(1970,1,1),'dd-mmmm-yyyy HH:MM:SS.FFF'); % Save fractions of seconds, as well
    eventMarkers(i).duration  = fread(h,1,'double');
    fseek(h,48,'cof');
    evtUser                       = fread(h,12,'uint16');
    eventMarkers(i).user      = deblank(char(evtUser).');
    evtTextLen                    = fread(h,1,'uint64');
    evtGUID                       = fread(h,16,'uint8');
    eventMarkers(i).GUID      = sprintf('{%.2X%.2X%.2X%.2X-%.2X%.2X-%.2X%.2X-%.2X%.2X-%.2X%.2X%.2X%.2X%.2X%.2X}',evtGUID([4 3 2 1 6 5 8 7 9:16]));
    fseek(h,16,'cof');    % Skip Reserved4 array
    evtLabel                      = fread(h,32,'uint16'); % LABELSIZE = 32;
    evtLabel                      = deblank(char(evtLabel).');    % Not used
    eventMarkers(i).label = evtLabel;
    
    % Only a subset of all event types are dealt with
    switch eventMarkers(i).GUID
        case HCEVENT_SEIZURE
            eventMarkers(i).IDStr = 'Seizure';
            %disp(' Seizure event');
        case HCEVENT_ANNOTATION
            eventMarkers(i).IDStr = 'Annotation';
            fseek(h,32,'cof');    % Skip Reserved5 array
            evtAnnotation = fread(h,evtTextLen,'uint16');
            eventMarkers(i).annotation = deblank(char(evtAnnotation).');
            %disp(sprintf(' Annotation:%s',evtAnnotation));
        case HCEVENT_FORMATCHANGE
            eventMarkers(i).IDStr = 'Format change';
        case HCEVENT_PHOTIC
            eventMarkers(i).IDStr = 'Photic';
        case HCEVENT_POSTHYPERVENT
            eventMarkers(i).IDStr = 'Posthyperventilation';
        case HCEVENT_REVIEWPROGRESS 
            eventMarkers(i).IDStr = 'Review progress';
        case HCEVENT_EXAMSTART
            eventMarkers(i).IDStr = 'Exam start';
        case HCEVENT_HYPERVENTILATION
            eventMarkers(i).IDStr = 'Hyperventilation';
        case HCEVENT_IMPEDANCE
            eventMarkers(i).IDStr = 'Impedance';
        otherwise
            eventMarkers(i).IDStr = 'UNKNOWN';
    end
    
    % Next packet
    offset = offset + pktLen;
    fseek(h,offset,'bof');
    pktGUID = fread(h,16,'uint8');
    pktLen  = fread(h,1,'uint64');
end
end

function [montage] = read_nervus_header_montage(h, StaticPackets, Index)
%% Get montage  - Andrei Barborica, Dec 2015
% Derivation (montage)
mtgIdx  = StaticPackets(find(strcmp({StaticPackets.IDStr},'DERIVATIONGUID'),1)).index;
indexIdx      = find([Index.sectionIdx]==mtgIdx,1);
fseek(h,Index(indexIdx(1)).offset + 40,'bof');    % Beginning of current montage name
mtgName       = deblank(char(fread(h,32,'uint16')).');
fseek(h,640,'cof');                             % Number of traces in the montage
numDerivations = fread(h,1,'uint32');
numDerivations2 = fread(h,1,'uint32');

montage = struct();
for i = 1:numDerivations
    montage(i).derivationName = deblank(char(fread(h,64,'uint16')).');
    montage(i).signalName1    = deblank(char(fread(h,32,'uint16')).');
    montage(i).signalName2    = deblank(char(fread(h,32,'uint16')).');
    fseek(h,264,'cof');         % Skip additional info
end

% Display properties
dispIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'DISPLAYGUID'),1)).index;
indexIdx  = find([Index.sectionIdx]==dispIdx,1);
fseek(h,Index(indexIdx(1)).offset + 40,'bof');    % Beginning of current montage name
displayName = deblank(char(fread(h,32,'uint16')).');
fseek(h,640,'cof'); % Number of traces in the montage
numTraces = fread(h,1,'uint32');
numTraces2 = fread(h,1,'uint32');

if (numTraces == numDerivations)
    for i = 1:numTraces
        fseek(h,32,'cof');
        montage(i).color = fread(h,1,'uint32'); % Use typecast(uint32(montage(i).color),'uint8') to convert to RGB array
        fseek(h,136-4,'cof');
    end
else
    disp('Could not match montage derivations with display color table');
end
end


function [montages] = read_nervus_header_dynamic_montages(DynamicPackets)
%% Get montages from the dynamic packets  - Jan Brogger - September 2016

montagePackets = DynamicPackets(strcmp({DynamicPackets.IDStr},'DERIVATIONGUID'));
montages = struct();
for numMontage = 1:size(montagePackets,2)
    montagePacket = montagePackets(numMontage);
    
    guidmixed = montagePacket.data(1:16);
    guidnonmixed = [guidmixed(04), guidmixed(03), guidmixed(02), guidmixed(01), ...
        guidmixed(06), guidmixed(05), guidmixed(08), guidmixed(07), ...
        guidmixed(09), guidmixed(10), guidmixed(11), guidmixed(12), ...
        guidmixed(13), guidmixed(15), guidmixed(15), guidmixed(16)];
    montages(numMontage).guid1 = num2str(guidnonmixed,'%02X');
        
    montages(numMontage).packetSize = typecast(montagePacket.data(17:24),'uint64');
    
    guidmixed2 = montagePacket.data(25:40);
    guidnonmixed2 = [guidmixed2(04), guidmixed2(03), guidmixed2(02), guidmixed2(01), ...
        guidmixed2(06), guidmixed2(05), guidmixed2(08), guidmixed2(07), ...
        guidmixed2(09), guidmixed2(10), guidmixed2(11), guidmixed2(12), ...
        guidmixed2(13), guidmixed2(15), guidmixed2(15), guidmixed2(16)];
    montages(numMontage).guid2 = num2str(guidnonmixed2,'%02X');
    
    montages(numMontage).itemName = deblank(cast(montagePacket.data(41:104),'char')');    
    montages(numMontage).elements = typecast(montagePacket.data(745:748),'uint32');
    montages(numMontage).alloc = typecast(montagePacket.data(749:752),'uint32');
    
    offset = 753;
    for i=1:montages(numMontage).elements        
        montages(numMontage).channel(i).name = deblank(cast(montagePacket.data(offset:(offset+127)),'char')');
        offset = offset + 128;
        montages(numMontage).channel(i).active = deblank(cast(montagePacket.data(offset:(offset+63)),'char')');
        offset = offset + 64;
        montages(numMontage).channel(i).reference = deblank(cast(montagePacket.data(offset:(offset+63)),'char')');
        offset = offset + 64;
        montages(numMontage).channel(i).isDerived = montagePacket.data(offset);
        offset = offset + 1;
        montages(numMontage).channel(i).isSpecial = montagePacket.data(offset);
        offset = offset + 1;
        offset = offset + 256;        
        offset = offset + 6;
    end        
end

end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     002222375895663792223 -0.001883949592848320779
'MRT23-1706'    -0.08340999993678881175 -0.06747303878644381192 -0.01504468196754400397 0.002958614061142303495 -0.006189050101080968008    -0.005640230797651546793    0.0002591916738018565738    0.002917317441589548782
'MRT24-1706'    0.03171168279119977035  0.04365824586675463165  0.04667283211934530884  0.000216492168100109177 0.00274607520382899508  0.003346553492217337576 9.575290406123761601e-05    -0.002637571559936767104
'MRT25-1706'    -0.09910903306266298507 -0.01001841365376395995 -0.07960848962092809122 0.002469112155534166145 -0.003822984385791183912    -0.008347447632039494553    0.001898502888129895627 0.003940294515481054978
'MRT26-1706'    0.0007734276721666665444    -0.06700129664722917489 -0.04038920592138294657 0.0003854113951465349048    -0.001563770739720348911    8.564504304294472126e-05    -0.00232512996713209099 0.00269018265925205155
'MRT31-1706'    0.04050681318299218653  0.03125316201216151585  0.001506916185139704854 -0.003533473004486763474    0.001414057153773838746 0.004750884268308786136 0.0005100418475276040104    -0.002746408068342049997
'MRT32-1706'    0.08217026326770715539  0.02448282558704611978  -0.07416495503046548199 -0.007868065332407993967    0.003291684021555728572 0.002770759820770513863 0.002439647833194805767 9.25362291549575519e-05
'MRT33-1706'    0.007427308952660549059 -0.01154170313289417348 0.01369371728950028247  -0.0001564389549000813249   -0.0009190814209578675002   0.0001141384526413533919    -0.0007738174915289131603   -7.83291616818462713e-05
'MRT34-1706'    0.03811942565064375954  0.09709976587065491382  -0.1165401771927899177  -0.007682893268070568882    0.005641741760780463397 -0.001578302192996052476    0.0101381861130977538   -0.002297401560141792889
'MRT35-1706'    -0.03387641660861703125 -0.07724995311540595877 -0.01147930860206294    0.0003949886581202345583    -0.003508316223875716615    -0.003308568367685030005    -0.002868273518264801477    0.006596507278721489183
'MRT41-1706'    -0.1758567692921810532  -0.03515977917188119889 -0.01519512535015960994 0.01134698594789500248  -0.00276902816365899888 -0.01341278181279694642 6.535853401789164057e-05    0.002465644192952843282
'MRT42-1706'    -0.02152130408445456292 -0.03205011379447719744 0.0009744124587214717447    0.0006857491666887961584    -0.002992507465501631752    -0.002249384464828016027    0.0002026534437485450564    0.001707255060233871247
'MRT43-1706'    0.09889795875194561103  0.003888088698199017312 0.01979261571006893025  -0.005854555830175177736    0.002550630691211766554 0.007814062754916960632 0.001534605528027466453 -0.001740738529028019663
'MRT44-1706'    -0.05691008159632153507 -0.1645107481226831536  -0.04779496534737437408 -0.004393966613594979531    -0.008692954536789214773    -0.00741345907509721215 -0.001635450184394400455    0.01400487427236647882
'MZC01-1706'    -0.01494286518344861876 0.09399604291418038604  0.001514042424394920847 -0.0002662129026606222925   0.00272306583656161446  -0.001462974992030099041    0.000286755318048692229 -0.002427846703108969481
'MZC02-1706'    0.01048523110613878463  0.0009020160280337854216    -0.05016703820507201206 -0.0002998510518858664388   -0.000515313659549525337    -0.003264881453735892902    0.00102562807721430743  0.0008618217559544430603
'MZF01-1706'    0.01968281285050175675  -0.06023107763614371502 0.06121508202898146705  0.0001622292073713806819    -0.003467093553142917405    0.00172970768133278211  0.0002289916838627152932    0.0003681591649551428712
'MZF02-1706'    0.06678984472783246196  -0.182440558735761621   -0.09178606635743639941 -0.007584419234038608168    -0.01038819388683161551 -0.002998489417295726635    0.007930818537335527704 0.004545979240750585605
'MZF03-1706'    0.05325156675761270192  0.1174738411772354302   0.02751730704591508897  -0.001211262822488221903    0.003903677346005184935 0.002496012419996640041 -0.003155274508161381476    -0.001828277363740007846
'MZO01-1706'    -0.02439720069461531055 0.06592929524778615158  -0.07484527118930360545 -0.001991400967305742824    -0.002043057608432549601    0.0009466362350490286834    0.002905927484150935911 0.0003895711563017285189
'MZO02-1706'    -0.004193409439129229879    -0.03588347742750992719 -0.01056765234081206463 0.000828269113155421957 6.357854193919900885e-05    0.0005194924160432156995    -0.0005895915780497043351   0.002171980319106391489
'MZP01-1706'    0.03381324548225054377  -0.0969433821121011563  0.08366330975129604441  0.002633418695614468154 0.0001724388743565483367    0.001263530525531823259 -0.004448313057072220827    0.0004818693201072772776
'MZP02-1706'    0.002633998031293571709 0.1230252286459674327   -0.1066414172307003361  -0.005316138877684267912    0.0008464339200343246094    0.0006126163070335105586    0.007819364143145711085 -0.001852859759265136543
};
coef(end).label = temp(:,1);
coef(end).value = cell2mat(temp(:,2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coef(end+1).name = 'G2BR';
coef(end).reference = {'BP1-1706', 'BP2-1706', 'BP3-1706', 'P11-1706', 'P12-1706', 'P13-1706', 'P22-1706', 'P23-1706'};
temp = {
'MLC11-1706'    -0.1382228139973263881  -0.0662533451338533641  -0.08514893299283234074 0.2886713463945499991   -0.03145826574391359653 -0.1032094036312776741  -0.2677078788695917422  -0.1726664652951542811
'MLC12-1706'    -0.1725242892252873694  -0.02316799131720484115 -0.09864268844080101295 0.2708002617960945102   0.004966079429890845123 0.006366494614789523429 -0.2078166773781876531  -0.2548113407111339268
'MLC13-1706'    -0.06433513939400752057 -0.02469987208653662281 -0.1547419125782673044  0.1418973878845695979   0.1154807762403902011   0.1188750122484098504   -0.02036827737249559053 -0.3102419760200211818
'MLC14-1706'    -0.0101725826309732939  0.07652184257691851954  0.06285597934187520108  0.02169242670745954657  -0.2582685301006519274  -0.1754755617201763096  -0.0770185861151197193  0.2127978003944720475
'MLC15-1706'    -0.0234806731320204784  0.02967070404938945694  0.04419021236823218041  -0.1472778785094536191  0.3208338160851761933   0.1096039676317745737   0.1155024694090843046   -0.0999885415781458925
'MLC21-1706'    0.08743625635265907248  0.04198538522939099238  0.002595763142494787734 -0.2948187247460258842  0.008232970318121625442 0.02056178412134233163  0.2672628912666229484   0.1725954757795352679
'MLC22-1706'    -0.1860618882723089618  0.01115782510601683909  -0.07663727931105152047 0.2479503631283714593   0.04773637890704240938  0.08225074763523326782  -0.1824246449791908642  -0.2695854761088272711
'MLC23-1706'    0.09932052461823512313  0.02459404557479431819  0.1243893787964647085   -0.1214521353695743461  -0.1467888516667933962  -0.1954673507175130331  0.0827014306951008582   0.2565056823171180667
'MLC24-1706'    0.09542051075649532477  -0.09846890209426244212 0.1029865957655926867   0.04568292306962887933  -0.2578532484170544814  -0.212186495553021659   -0.026742951757970361   0.20175432492741463
'MLC31-1706'    0.05360363065430327062  -0.0290658497723569817  0.08444762605308975267  -0.2746546031050620673  -0.02850152883999820461 -0.08464483138426905084 0.2478561306001653697   0.1935575302880424842
'MLC32-1706'    0.01417751084994100727  -0.009609644382386633293    0.1893266810670302669   -0.1537872939301765385  -0.1114998313699658156  -0.2384590696179459679  0.1708144723346288196   0.2018385584725512771
'MLC33-1706'    -0.07741406054111449331 0.140602962564341849    -0.1613843832943822887  0.1657525588942669881   0.01936887116091059666  0.2843437265599630948   -0.2427637099503201223  -0.03446355440781331131
'MLC41-1706'    -0.05205615869574087978 0.04122180525868625323  -0.2040984775909389648  0.306050945666863472    -0.004845563660082926466    -0.02656975327913183027 -0.2965742741062059595  -0.1030378866416584993
'MLC42-1706'    -0.02389402919451669161 0.07844658973459740781  -0.1838052010794586355  0.2877669212562368739   0.01362667991036912302  0.1231018548508429622   -0.2956834717016642577  -0.06702306755529127691
'MLC43-1706'    -0.005209111846768275933    0.03384509470499676154  -0.198775722190652343   0.233533809231346895    -0.04116800311516879973 0.204014512418244004    -0.257621450547025066   0.1104885236548849847
'MLF11-1706'    0.06358350010679529596  -0.0784445226147908875  0.07465856520074432012  -0.06599042869640468767 -0.2771210257315528946  -0.006568449735539011285    0.2547869704665173263   -0.02066346292026897263
'MLF12-1706'    0.05193664472756427675  0.007433386036594913024 -0.004354256429245254717    -0.0104594934129525613  -0.1185032398077054266  -0.004922107449456929358    0.3094693534941974078   -0.03052891107477028951
'MLF21-1706'    -0.094404553516449774   0.05713796765517231008  0.08266203354542170367  0.09905781386387088605  0.3179322781523764574   0.06462239355808441221  -0.2089744688446816823  0.09630895339785545961
'MLF22-1706'    -0.005401466460099793067    0.07799861220626723435  0.03423511561012176907  0.02010940732427979327  0.1913544074894210911   0.03221201176237215025  -0.2826473303848106022  0.09805244526818387985
'MLF23-1706'    0.01742441361442586059  0.008819523871757835309 0.09118678122227265315  -0.00035833527450891895 0.05733828285980938771  0.001937081080728298614 0.3197619531956262806   -0.02787631604312304501
'MLF31-1706'    0.04654583340223795995  -0.047033385721240481   -0.03064091442360981216 0.01509156098077189812  -0.2121415682595462204  -0.1112084492731538476  0.1523790428804583652   -0.240462429525251653
'MLF32-1706'    0.04746952247722049933  0.02258704023958208601  0.02194903989996015853  -0.0414873715701736534  0.0450056070330318303   0.02014421752427554557  -0.239552750474094106   0.2240642801760762137
'MLF33-1706'    0.05869476665463646547  -0.08361178423739487608 0.1122550523316090815   0.01652918328603633635  -0.2009617190249260721  -0.03349511899062078413 -0.2660305009062145509  0.1264355272430096322
'MLF34-1706'    -0.0218949579017706196  0.04946250378941413617  -0.1116108859354522365  -0.0491185333376887856  0.2506667734195127117   0.01359099838777923404  0.2722985984247927593   -0.04207272814955062556
'MLF41-1706'    -0.03942545908135971361 -0.02407661327905865484 -0.1064493759388815408  0.1087447642448455681   -0.1595908287940409898  -0.2032976679415572085  -0.06620409544320542128 -0.2639041423904764372
'MLF42-1706'    0.05363422328221816027  0.07967140221233734443  0.08951534003717770416  -0.128937214577220377   0.08810343302557513412  0.07615585441940531108  -0.04098182085352847226 0.3213350654240125781
'MLF43-1706'    -0.06803485675123928378 0.02096578177021040518  -0.07395732983967347185 0.05434491795285490562  0.1411017643030943458   0.07179673034290874611  0.1643525586201984023   -0.2682898473340409873
'MLF44-1706'    -0.06260930679826015588 0.04295959449824775928  -0.07156891208540423766 0.006612428169781000946 0.2081223075433809855   0.08387857375424689899  0.1932466340951884409   -0.2172583991596181707
'MLF45-1706'    -0.05273449228516488563 -0.008524618559531807471    0.02247780394187684841  -0.1073240959812722661  0.3276726083020715108   0.07029252766680937903  0.1898373308456635378   -0.09345150444010934254
'MLF51-1706'    0.1624317912595521218   0.08610496267367034884  0.2761220436858867755   -0.2291957786065656477  0.06302140921111969163  0.1150472021390444161   0.1585149441735389519   0.2768588909419009703
'MLF52-1706'    0.07813537475270850019  0.01258956549997311894  0.0704459322986105696   -0.1925680107901989457  -0.04418983798888989623 -0.05323237435557155728 0.0706631643303960999   0.315963044331685039
'MLO11-1706'    0.001623072866654050317 -0.0505657566297947822  0.08349523534917566159  0.1933396783320261991   0.3210745787150182107   -0.10751155424449374    -0.1046291748235521901  -0.07963442963719634393
'MLO12-1706'    -0.0932543347501357256  0.0816296319231611528   -0.03166616724170764202 -0.270890689883850011   -0.204853914567168699   0.1259003486000364325   0.01927732627214839151  0.04117169739736232376
'MLO21-1706'    -0.01017560646354364814 0.01985941627204387133  0.05514742758507345383  -0.2571367677595031709  -0.2793730885099177663  0.02122191927113283028  0.06579356396974349319  0.01006261887051268511
'MLO22-1706'    0.04228339699326733037  0.01212366091993068162  0.08879376494519121421  0.3260521580611809855   0.07925669214609561142  -0.01819740163416505921 -0.009591517815856168139    -0.003271903174406161872
'MLO31-1706'    0.01883014372343630102  -0.1412157130247556902  -0.03076257068120249322 -0.2120135195069595135  -0.3234457946035887432  0.01496756686803598843  0.114501508697362378    0.0212206696608342113
'MLO32-1706'    -0.06075735201504169763 -0.05558799106456796535 0.1118396305158921855   0.3025245377843828187   0.1916345478112770895   -0.0186420866884403659  -0.03524606438421783866 -0.002627656002018843341
'MLO33-1706'    0.04693372383741226117  0.08044034771046811028  0.0653618770327214077   0.3384757636250301793   -0.0258698388823904446  -0.009578034291444910767    -0.005871694612373942572    -0.005290584850912462972
'MLO41-1706'    0.03923384480471365371  -0.01525029131105908736 0.02714999417571031803  0.2546115164329926239   0.2690482753680357164   -0.009767390421729766486    -0.06730807681393313757 -0.004156501435087131199
'MLO42-1706'    0.001428077299708022108 0.003844081474509456094 0.01959334188366524537  0.3214212593155210129   0.09069845192517340948  -0.009534190431742349883    -0.01124799407868778967 -0.003469786724388063733
'MLO43-1706'    -0.07099233623833486884 -0.02880924754568231147 0.08856351332414597044  -0.3060528484596556353  0.128839616138026275    0.003165219732926745702 0.00994938826448356177  0.0006090051591933717291
'MLP11-1706'    -0.03249607645747909501 -0.03205139619401013718 -0.02829538554438305181 0.07148907765356921074  -0.1146178887559892168  0.2960008321297108447   -0.1589520302232046711  0.140935032443428615
'MLP12-1706'    0.07426184295142532199  0.1492185896179328841   -0.05111344519783470042 -0.06136746941840764275 0.001724787166957722582 0.3282088611462577465   -0.130134171293704165   -0.004936585729251252545
'MLP13-1706'    -0.1550668393682321788  0.1438313884956307909   -0.1629892686117093148  -0.0394500456953563311  0.1701208339783343848   0.2852099729622109003   -0.08313554870547654185 -0.14560583845332567
'MLP21-1706'    0.02654843949261885203  0.04836841501601622101  0.1462326737279008859   0.03725154422432774837  0.2241878770131820464   -0.2473136379978512212  0.02254832068650799645  -0.1907638261900793664
'MLP22-1706'    -0.0686176390587671936  -0.02320309969255752308 0.001400987193962166446 -0.09754388151307355481 -0.1665353760105017233  0.2813326097644027191   -0.06635495897286450284 0.1135272370178197487
'MLP31-1706'    0.00410461651847717629  -0.05408423003028114961 0.1377390235016080489   0.1692504652513706909   0.2497298125801123425   -0.2136456344231150439  -0.02055520016833128211 -0.1160272495532611442
'MLP32-1706'    0.06336749823488235789  -0.05042684966902965388 -0.01654673930951551519 0.236575211873756186    0.09534624686072780408  -0.2177310220670032048  0.03421893273639377286  -0.03397234014568384364
'MLP33-1706'    -0.03323736611209159725 0.03002182472506493074  0.02689757791333211023  -0.2207079389025662641  0.128711273619539962    0.2266871851504149993   -0.03736291801907636045 -0.05088352815167572485
'MLP34-1706'    0.01804484165856728647  0.01031036791405627942  0.07549176520119624256  -0.1699046426367207818  0.2899420145878607857   0.1446179218179321768   0.06650532942325219909  -0.0992577900696323262
'MLT11-1706'    -0.04704491941066894084 -0.03797941156058894097 -0.07737316306722230586 -0.0371319234728271913  0.2231138592719367209   0.003212535244369423872 0.277973060208098588    -0.02185444193729014356
'MLT12-1706'    -0.08447872404705902838 -0.03003933274360280972 -0.08427636400014384965 -0.08897256645086894233 0.309660363262441618    0.007748120488902505119 0.2347308524974138499   -0.02053834692366506356
'MLT13-1706'    -0.04155250655432526918 -0.01767708282357837238 -0.005058913529922428806    -0.1744857916782741325  0.3460584289539624936   0.01625148968725365073  0.1543297160145147839   -0.01757610927126036088
'MLT14-1706'    0.02473601543887470749  0.09954298179551401837  0.12087844240992629 0.1964107424914252886   -0.3303216413858801603  -0.01522956155732175729 -0.126241152957685121   0.009533979930416870002
'MLT15-1706'    0.09975884082223460125  -0.05496818271985834392 0.04251169561728222746  0.2673987405566053965   -0.2452043039081176767  -0.03954824050132876162 -0.04696591748658365123 0.02039530731989156603
'MLT16-1706'    -0.01809030105728349699 0.03788238383809146292  0.0911475289757280871   -0.3155033569353283829  0.04735155161616331349  0.08297294674179991891  -0.005384807996269442947    -0.006461459956578846642
'MLT21-1706'    -0.1288538542114405494  0.002481952925272324847 0.03485703117358385922  -0.009551033352974021978    0.1422676410576970696   -0.004939041608370217835    0.3031750340360210294   -0.02265133277841845094
'MLT22-1706'    -0.0237058575720274084  0.04512844725750929376  -0.006987253626539095472    0.05103612862324776595  -0.2530604354950670154  -0.009800139935999267055    -0.2678564659799040348  0.01978457309890234092
'MLT23-1706'    0.01363484305779106079  -0.03864465828962739685 0.02261442090285450263  -0.1509838934441197489  0.3446303859724349494   0.0201338703042538876   0.1770152141325894435   -0.01924419428211009142
'MLT24-1706'    0.01393201296458278921  0.0305976414413392378   -0.01880820659807747042 0.1839983223451140404   -0.3361884590426196251  -0.01389905054059030476 -0.1375696812691714344  0.01182448961353790157
'MLT25-1706'    0.06463614949607772442  -0.009091585275117273698    0.01958275607012524913  0.2400990650349234534   -0.2953646988618118252  -0.01301339264680832254 -0.08030090445001912547 0.008743103831794342071
'MLT26-1706'    0.08681852359305899935  -0.04323256250900678332 -0.07327380800820051943 0.3092413470854838975   -0.1484321325726151719  -0.01074461506969169935 -0.0151504182588766325  0.003048914048373080506
'MLT31-1706'    0.003035953703011757404 0.0646580228454244621   0.06415362135829164036  0.04362663534979851349  -0.2322643812263140428  -0.004476755280500705793    -0.2705967667099912766  0.01602903942775381341
'MLT32-1706'    -0.03739804724905297639 -0.07804379670670040514 -0.07195296403668691165 -0.09604611440909098519 0.3052418652000546406   0.006742533954966011857 0.2313814044002251435   -0.01405236708920687763
'MLT33-1706'    0.01828050731997953038  0.07289809325464820244  0.00657024311860331979  0.1769381056866486668   -0.3549014141525700383  -0.009858498387001077279    -0.1586977712752264957  0.00536515252067098266
'MLT34-1706'    -0.0006301953427265624685   0.03576772048592252035  -0.01198610988375524086 0.2152738178167998007   -0.3307915708951239542  -0.01134930826208157238 -0.1127700732372116377  0.005556709078458120163
'MLT35-1706'    0.05334568942542123465  0.01680234893714786792  0.1533550635366545123   0.2992700536244244369   -0.2357057131142628326  -0.0111909195602401014  -0.04976395803997243211 0.004429677610622913264
'MLT41-1706'    0.08651884553088011465  0.01301954871435762398  0.0461874236875364913   0.05370947177516879889  -0.2573876221729367031  -0.0004022934679234792255   -0.2578139585112749765  0.0181297648629859115
'MLT42-1706'    0.03185784436153658861  -0.04740810871268212601 -0.0697960927546261295  -0.1473082410533657349  0.3355281876026416366   0.01029525218719577404  0.1834483212803902363   -0.009071629251501650559
'MLT43-1706'    0.05372747600040422844  0.02026929145742300403  0.03885607146340949031  0.1801580011202781884   -0.3276590704776112251  -0.007393763903438013881    -0.132827908291525465   0.008161154239529398377
'MLT44-1706'    -0.004786981027231705654    0.0351094683864824203   -0.05814207825039392069 -0.2452282195270626131  0.281805606880582149    0.009707600877972163977 0.07841898469102287494  -0.008444893768065512424
'MRC11-1706'    -0.0007879057259526575652   0.0562766694771630327   -0.1606509791041782764  0.2625214587867037985   -0.02943882930227824379 -0.1765816276048170208  -0.2781150695457602628  -0.1126389196907389145
'MRC12-1706'    0.005041312833811165128 -0.06072117216660291583 -0.05206699265535509602 0.2086495590651812615   -0.003767579200974270279    -0.2391770468357752033  -0.2618243726129056981  -0.004602653674354887856
'MRC13-1706'    0.0210067559732514085   0.01087679324347014254  -0.05049581465957574122 0.01484813121440720085  0.121673549590438107    -0.3081245755028024225  -0.1336690921949419963  0.1230624412539903928
'MRC14-1706'    -0.1626377966401037489  0.0734525194434082207   0.1031786964687588237   0.08959637570570466725  -0.2602790116005987975  0.218979276625638819    -0.02874054868554443143 -0.1788435128919033945
'MRC15-1706'    -0.1843432827081813996  -0.06793325587954066047 0.0206622197920528633   0.1270255143763397709   -0.3312232160428612571  0.1019088645291856693   -0.1542841425969307123  -0.1151075077535118479
'MRC21-1706'    -0.1624387619387309689  -0.05665387335815671999 -0.1797641205371675044  0.275760674553589713    -0.006878584310087037731    -0.180999150896493266   -0.2951779406261597072  -0.01604416014587845343
'MRC22-1706'    0.008115207652171683994 -0.02539517701411761474 -0.1537991189617910759  0.1807428031374978261   0.04658037739398068755  -0.2681088102671762563  -0.240935352738834152   0.08616231667606881472
'MRC23-1706'    0.1373696362710403462   -0.006754250254273726044    -0.02322382248563889875 0.08395796000441761719  0.1489590287865106  -0.2543618894601860636  -0.1241346125955813534  0.1985899454334250913
'MRC24-1706'    -0.1322683152887810243  -0.1253515894290049559  0.1446953244022337604   0.02769056553012722716  -0.2584609775592148773  0.1966702959277972507   -0.04996274691785355243 -0.2170447783805828645
'MRC31-1706'    -0.03351178862986126106 -0.02502474764367358315 -0.1354577267289230336  0.2487483637503605793   0.02800565614940635217  -0.1888168006782719177  -0.2715657568930669452  0.08482959003495658956
'MRC32-1706'    0.09422505621642351947  0.0974268035931084575   0.01434750330881125191  0.1759475490781015061   0.1147891059968733918   -0.1882355818634162681  -0.1525823291225343881  0.2385474019058428652
'MRC33-1706'    -0.02743513551602546133 -0.001346386555345108205    0.1040961934368644787   -0.2389427230023080662  -0.0210474394704532955  0.03072307783726403035  0.160168753914748635    -0.2825274975756738027
'MRC41-1706'    0.004069615769342165013 0.01524889356876139163  -0.2066967427690987735  0.2944317960398715739   -0.005128244934604534074    -0.1038586236046362271  -0.3021520579721332367  -0.02904603565774256657
'MRC42-1706'    -0.09261069535108770734 -0.05391463190176774833 0.07393599397715618082  -0.2967904418944345424  -0.01839631920055650668 0.0676362606430673946   0.2893330539464480933   -0.1241483698974569644
'MRC43-1706'    0.02121592901906437303  0.002736008405995041248 0.05867749995983034123  -0.2599204191479295711  0.04181595136570286086  -0.1132943588742273827  0.2334756974634588367   -0.2026394493973507904
'MRF11-1706'    -0.08237912727941071322 0.1930170658866355027   -0.02999511644876571778 0.2560554149778885846   0.2816132118458873901   0.02238669585085677899  -0.06731832491388607675 0.003857697094186068008
'MRF12-1706'    -0.1452633744918867542  0.0413703610092215554   -0.03966250163167530512 0.3163134685732207907   0.1070690457748723207   0.02720564582295028136  -0.007841875063684807315    0.00389592610458493584
'MRF21-1706'    0.09178898921567094082  -0.07326030686495163535 0.08279190922802177888  -0.2031894200380436055  -0.3079692909353750618  -0.0900654060700137038  0.08916715451558471228  -0.06269989707734738293
'MRF22-1706'    -0.1201093945905412491  0.1666597350767906605   -0.1266918845073792821  0.2800702632171455653   0.1963824608020995344   0.09204699582438576055  -0.01677482812683377122 0.02336377569770818619
'MRF23-1706'    0.1020462997170876462   -0.04476148723999442275 0.1125374056494373165   -0.3225737173313791839  0.0573288654152874777   -0.02539360876236581438 0.002355324897689480695 0.004255860165688138755
'MRF31-1706'    0.04460755158466561676  0.03889906756587595271  0.04177857900706357513  -0.1552087321145515231  -0.2052526905362134779  -0.2388507361215514435  -0.02003904756496036055 -0.110713764729872155
'MRF32-1706'    -0.1130533484546892981  0.03040892813732698485  -0.007850927498421034212    0.2429085798779209715   0.04094478195904468348  0.2200857368082461407   0.03784565422790188427  0.01452989250527078439
'MRF33-1706'    0.05205094682228531 0.004683024081423316914 -0.04904272170391472258 -0.2709531288968970997  0.2013610811022400471   -0.1260709754696639762  0.01923096933089796581  0.04400575957673252131
'MRF34-1706'    -0.1976976839368163597  0.008096467905396718956 -0.07168205328873483717 0.2845413282209070527   -0.2524060671096146069  0.03109573688763932892  -0.0492872162319630755  -0.01768523987644998047
'MRF41-1706'    -0.01709039766296914162 -0.1242241645463833533  0.05760983901736607682  0.06728767597205795314  -0.1646797689443417756  -0.260967487800225495   -0.1096413467969140909  -0.1927757931200057839
'MRF42-1706'    0.01785010709892954212  -0.02780013680753913585 -0.1056340188579661943  -0.04193121114594819399 -0.07854366550232408373 -0.3264883544679060834  -0.1283969803927337838  -0.06592678160634457551
'MRF43-1706'    -0.1460013779933039824  -0.008416346706667056365    -0.03687283923864802943 0.163355218596175078    -0.1412543970141160177  0.2692111348867542997   0.05720260391292623925  -0.07662781704671199623
'MRF44-1706'    0.05257004161810163723  0.09384741603594291826  -0.011268845423295204   -0.190943968764359373   0.2101659028426117737   -0.2186085251700451426  -0.005747623259752606892    0.08681586790408382659
'MRF45-1706'    0.07199443390726908976  0.007034553280036182424 -0.06714184252894038474 -0.1840608704768749992  0.3227559652076311125   -0.09959805645941838259 0.107333823150499913    0.07921069084164619745
'MRF51-1706'    -0.00913693632605947631 0.03027867128393653406  0.1652189968253011398   -0.168070378189555808   0.06081735921405818923  0.2681085377762136845   0.2261921790313926039   0.1069252315962646721
'MRF52-1706'    0.07649395040977598204  -0.05309537920131012639 -0.147506915141427275   0.06143852183256792721  0.04835541207699078914  -0.3203646898990325509  -0.1854305377398010046  0.06147484919150039695
'MRO11-1706'    0.06910029908177262037  0.03545983650013832811  -0.02158888727785019529 -0.1020278751344255996  -0.3175560286413296462  0.08403824520574382229  0.1853070260950459291   0.1048268365378309114
'MRO12-1706'    -0.08793773769676554997 -0.0171440733558369876  -0.1298907460213891141  -0.02076220236820512655 -0.2090988154049477077  0.038203607832962512    0.2714144291478555182   0.1368635910393601052
'MRO21-1706'    0.0438694744555706076   -0.1092273371107565227  0.07793492706636623302  0.06955358762672385342  0.2774552828328116227   -0.005840745098776823144    -0.2636923358538589834  -0.01961517797204404695
'MRO22-1706'    0.03997476187009500748  0.06355656022146444206  -0.07888175138367259365 -0.01117412174041480799 -0.08469174072528504182 0.0003682080590293678422    0.3274628421764786967   0.02126712273219566809
'MRO31-1706'    -0.03620589722182501746 -0.02987930659140425327 -0.06213334914691368005 0.118407085133361456    0.3259932373850928933   -0.01070798936638819171 -0.1976188461788445094  -0.009302204064925326948
'MRO32-1706'    -0.09094560252909791864 0.03158194181827086999  0.0230778135580117158   -0.02822735310259225777 -0.2018575532181628718  -0.003045371570343902157    0.2842089016990146377   0.01019822220958127579
'MRO33-1706'    -0.06683738042582462735 -0.08101546494756467487 -0.09737747036912622334 -0.004012321966511324216    -0.01612191771294267661 -0.004108697475914506285    -0.3127554768216350145  -0.007036500232096088228
'MRO41-1706'    -0.1225363074102404359  -0.006166923769823740714    -0.03533575305960776763 -0.07596378520017350866 -0.2911935506207178892  -0.00211488118868098391 0.2517530950482770313   0.01466944197950798695
'MRO42-1706'    -0.1056737531664330432  -0.1792631923173793884  -0.1301742665729932191  0.001986541658350563823 0.09561490736448920169  -0.0066325803334242317  -0.3050403144744178174  0.0006606828106451717586
'MRO43-1706'    0.01417578775829909064  0.09833592717382419468  0.04164056376765437606  -0.006287842766885024318    0.1151264732907840288   -0.002040487800555591986    0.320971392498041852    0.007277238339984764409
'MRP11-1706'    0.07937136737516901908  -0.08194471897082386547 0.00911682031052066591  0.1623852560418392044   -0.111313111170271109   0.1443868143757177369   -0.07907093939012151129 0.2945295122235159768
'MRP12-1706'    -0.06666497082207931135 -0.1396905727416886678  0.1547207712218562903   -0.1268877644416627481  -0.002864878921118730291    -0.002126475673897641316    -0.07116602560667395494 -0.3271894550938625468
'MRP13-1706'    -0.07270631909909819335 -0.03602690840902320751 -0.0164057944139923105  -0.08076570241445495124 -0.173085186979941491   0.1360476641012798305   -0.0500671400106301942  -0.290216626621441931
'MRP21-1706'    -0.0313210310086282015  -0.1285000062876426319  0.04184803149379558856  -0.01984220168867284473 0.2310851181898425089   -0.1972613619388189354  -0.03632270041376398084 -0.2445844123905571987
'MRP22-1706'    -0.002997104464619428435    0.07090456386015700496  -0.05100851650664366133 0.06356271010231347163  -0.1725637712268349444  0.1194692718967358541   0.09891304721593767446  0.2830716257311182149
'MRP31-1706'    0.04405419856043098215  0.0641848636923310345   -0.01402512920915510719 -0.01486844389148727144 -0.2464546592395191749  0.1185010264911498246   0.1632689096493861458   0.2096961493572045143
'MRP32-1706'    -0.06877166500509926395 0.006233115323576100819 0.09947401374409642338  -0.03169758388264812266 0.0937861738658565397   -0.03664328933815391709 -0.2348250697569170542  -0.2271607089027971038
'MRP33-1706'    -0.08348282514167326696 -0.1002621734012028726  -0.03903049761368913156 -0.03219885497932736124 -0.1309134925510340997  0.04797249600258152169  -0.2240468287700749517  -0.2235606873069985734
'MRP34-1706'    0.1102335224006878445   0.0674893320915589956   -0.0362946234611954896  -0.06044850207409657034 0.2918716097737355519   -0.1104362190697510127  0.1708385718994873403   0.166809603947957602
'MRT11-1706'    -0.09724202022874177398 -0.0266421507343533176  -0.07346252924992914546 0.2814838508913077297   -0.2369801585956532453  0.02309369494145910653  -0.04459200365448148795 -0.01055945544159708413
'MRT12-1