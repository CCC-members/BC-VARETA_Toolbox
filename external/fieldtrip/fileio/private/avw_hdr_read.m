function [ avw, machine ] = avw_hdr_read(fileprefix, machine, verbose)

% avw_hdr_read - read Analyze format data header (*.hdr)
%
% [ avw, machine ] = avw_hdr_read(fileprefix, [machine], [verbose])
%
% fileprefix - string filename (without .hdr); the file name
%              can be given as a full path or relative to the
%              current directory.
%
% machine - a string, see machineformat in fread for details.
%           The default here is 'ieee-le' but the routine
%           will automatically switch between little and big
%           endian to read any such Analyze header.  It
%           reports the appropriate machine format and can
%           return the machine value.
%
% avw.hdr - a struct, all fields returned from the header.
%           For details, find a good description on the web
%           or see the Analyze File Format pdf in the
%           mri_toolbox doc folder or read this .m file.
%
% verbose - the default is to output processing information to the command
%           window.  If verbose = 0, this will not happen.
%
% This function is called by avw_img_read
%
% See also avw_hdr_write, avw_hdr_make, avw_view_hdr, avw_view
%

% $Revision$ $Date: 2009/01/14 09:24:45 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%                    The Analyze format and c code below is copyright
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('verbose','var'), verbose = 1; end

if verbose,
    version = '[$Revision$]';
    fprintf('\nAVW_HDR_READ [v%s]\n',version(12:16));  tic;
end

if ~exist('fileprefix','var'),
    ft_error('...no input fileprefix - see help avw_hdr_read\n\n');
end
if ~exist('machine','var'), machine = 'ieee-le'; end


if contains(fileprefix, '.hdr')
    % fprintf('...removing .hdr extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.hdr','');
end
if contains(fileprefix, '.img')
    % fprintf('...removing .img extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.img','');
end
file = sprintf('%s.hdr',fileprefix);

if exist(file),
    if verbose,
        fprintf('...reading %s Analyze format',machine);
    end
    fid = fopen_or_error(file,'r',machine);
    avw.hdr = read_header(fid,verbose);
    avw.fileprefix = fileprefix;
    fclose(fid);

    if ~isequal(avw.hdr.hk.sizeof_hdr,348),
        if verbose, fprintf('...failed.\n'); end
        % first try reading the opposite endian to 'machine'
        switch machine,
        case 'ieee-le', machine = 'ieee-be';
        case 'ieee-be', machine = 'ieee-le';
        end
        if verbose, fprintf('...reading %s Analyze format',machine); end
        fid = fopen_or_error(file,'r',machine);
        avw.hdr = read_header(fid,verbose);
        avw.fileprefix = fileprefix;
        fclose(fid);
    end
    if ~isequal(avw.hdr.hk.sizeof_hdr,348)
        % Now throw an error
        if verbose, fprintf('...failed.\n'); end
        ft_error('...size of header not equal to 348 bytes!\n\n');
    end
else
    ft_error('...cannot find file %s.hdr\n\n',file);
end

if verbose
    t=toc; fprintf('...done (%5.2f sec).\n',t);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dsr ] = read_header(fid,verbose)

    % Original header structures - ANALYZE 7.5
    %struct dsr
    %       {
    %       struct header_key hk;            /*   0 +  40       */
    %       struct image_dimension dime;     /*  40 + 108       */
    %       struct data_history hist;        /* 148 + 200       */
    %       };                               /* total= 348 bytes*/
    dsr.hk   = header_key(fid);
    dsr.dime = image_dimension(fid,verbose);
    dsr.hist = data_history(fid);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hk] = header_key(fid)

    % The required elements in the header_key substructure are:
    %
    % int sizeof_header   Must indicate the byte size of the header file.
    % int extents         Should be 16384, the image file is created as
    %                     contiguous with a minimum extent size.
    % char regular        Must be 'r' to indicate that all images and
    %                     volumes are the same size.

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

    hk.sizeof_hdr    = fread(fid, 1,'*int32');  % should be 348!
    hk.data_type     = fread(fid,10,'*char')';
    hk.db_name       = fread(fid,18,'*char')';
    hk.extents       = fread(fid, 1,'*int32');
    hk.session_error = fread(fid, 1,'*int16');
    hk.regular       = fread(fid, 1,'*char')'; % might be uint8
    hk.hkey_un0      = fread(fid, 1,'*uint8')';

    % check if this value was a char zero
    if hk.hkey_un0 == 48,
        hk.hkey_un0 = 0;
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dime ] = image_dimension(fid,verbose)

    %struct image_dimension
    %       {                                /* off + size      */
    %       short int dim[8];                /* 0 + 16          */
    %           /*
    %           dim[0]      Number of dimensions in database; usually 4.
    %           dim[1]      Image X dimension;  number of *pixels* in an image row.
    %           dim[2]      Image Y dimension;  number of *pixel rows* in slice.
    %           dim[3]      Volume Z dimension; number of *slices* in a volume.
    %           dim[4]      Time points; number of volumes in database
    %           */
    %       char vox_units[4];               /* 16 + 4          */
    %       char cal_units[8];               /* 20 + 8          */
    %       short int unused1;               /* 28 + 2          */
    %       short int datatype;              /* 30 + 2          */
    %       short int bitpix;                /* 32 + 2          */
    %       short int dim_un0;               /* 34 + 2          */
    %       float pixdim[8];                 /* 36 + 32         */
    %           /*
    %               pixdim[] specifies the voxel dimensions:
    %               pixdim[1] - voxel width, mm
    %               pixdim[2] - voxel height, mm
    %               pixdim[3] - slice thickness, mm
    %               pixdim[4] - volume timing, in msec
    %                   ..etc
    %           */
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

    dime.dim        = fread(fid,8,'*int16')';
    dime.vox_units  = fread(fid,4,'*char')';
    dime.cal_units  = fread(fid,8,'*char')';
    dime.unused1    = fread(fid,1,'*int16');
    dime.datatype   = fread(fid,1,'*int16');
    dime.bitpix     = fread(fid,1,'*int16');
    dime.dim_un0    = fread(fid,1,'*int16');
    dime.pixdim     = fread(fid,8,'*float')';
    dime.vox_offset = fread(fid,1,'*float');
    dime.roi_scale  = fread(fid,1,'*float');
    dime.funused1   = fread(fid,1,'*float');
    dime.funused2   = fread(fid,1,'*float');
    dime.cal_max    = fread(fid,1,'*float');
    dime.cal_min    = fread(fid,1,'*float');
    dime.compressed = fread(fid,1,'*int32');
    dime.verified   = fread(fid,1,'*int32');
    dime.glmax      = fread(fid,1,'*int32');
    dime.glmin      = fread(fid,1,'*int32');

    if dime.dim(1) < 4, % Number of dimensions in database; usually 4.
        if verbose,
            fprintf('...ensuring 4 dimensions in avw.hdr.dime.dim\n');
        end
        dime.dim(1) = int16(4);
    end
    if dime.dim(5) < 1, % Time points; number of volumes in database
        if verbose,
            fprintf('...ensuring at least 1 volume in avw.hdr.dime.dim(5)\n');
        end
        dime.dim(5) = int16(1);
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ hist ] = data_history(fid)

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

    hist.descrip     = fread(fid,80,'*char')';
    hist.aux_file    = fread(fid,24,'*char')';
    hist.orient      = fread(fid, 1,'*uint8');  % see note below on char
    hist.originator  = fread(fid,10,'*char')';
    hist.generated   = fread(fid,10,'*char')';
    hist.scannum     = fread(fid,10,'*char')';
    hist.patient_id  = fread(fid,10,'*char')';
    hist.exp_date    = fread(fid,10,'*char')';
    hist.exp_time    = fread(fid,10,'*char')';
    hist.hist_un0    = fread(fid, 3,'*char')';
    hist.views       = fread(fid, 1,'*int32');
    hist.vols_added  = fread(fid, 1,'*int32');
    hist.start_field = fread(fid, 1,'*int32');
    hist.field_skip  = fread(fid, 1,'*int32');
    hist.omax        = fread(fid, 1,'*int32');
    hist.omin        = fread(fid, 1,'*int32');
    hist.smax        = fread(fid, 1,'*int32');
    hist.smin        = fread(fid, 1,'*int32');

    % check if hist.orient was saved as ascii char value
    switch hist.orient,
        case 48, hist.orient = uint8(0);
        case 49, hist.orient = uint8(1);
        case 50, hist.orient = uint8(2);
        case 51, hist.orient = uint8(3);
        case 52, hist.orient = uint8(4);
        case 53, hist.orient = uint8(5);
    end

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


% Comments
% The header format is flexible and can be extended for new
% user-defined data types. The essential structures of the header
% are the header_key and the image_dimension.
%

% The required elements in the header_key substructure are:
%
% int sizeof_header   Must indicate the byte size of the header file.
% int extents         Should be 16384, the image file is created as
%                     contiguous with a minimum extent size.
% char regular        Must be 'r' to indicate that all images and
%                     volumes are the same size.
%

% The image_dimension substructure describes the organization and
% size of the images. These elements enable the database to reference
% images by volume and slice number. Explanation of each element follows:
%
% short int dim[ ];      /* Array of the image dimensions */
%
% dim[0]      Number of dimensions in database; usually 4.
% dim[1]      Image X dimension; number of pixels in an image row.
% dim[2]      Image Y dimension; number of pixel rows in slice.
% dim[3]      Volume Z dimension; number of slices in a volume.
% dim[4]      Time points; number of volumes in database.
% dim[5]      Undocumented.
% dim[6]      Undocumented.
% dim[7]      Undocumented.
%
% char vox_units[4]     Specifies the spatial units of measure for a voxel.
% char cal_units[8]      Specifies the name of the calibration unit.
% short int unused1      /* Unused */
% short int datatype      /* Datatype for this image set */
% /*Acceptable values for datatype are*/
% #define DT_NONE             0
% #define DT_UNKNOWN          0    /*Unknown data type*/
% #define DT_BINARY           1    /*Binary             ( 1 bit per voxel)*/
% #define DT_UNSIGNED_CHAR    2    /*Unsigned character ( 8 bits per voxel)*/
% #define DT_SIGNED_SHORT     4    /*Signed short       (16 bits per voxel)*/
% #define DT_SIGNED_INT       8    /*Signed integer     (32 bits per voxel)*/
% #define DT_FLOAT           16    /*Floating point     (32 bits per voxel)*/
% #define DT_COMPLEX         32    /*Complex (64 bits per voxel; 2 floating point numbers)/*
% #define DT_DOUBLE          64    /*Double precision   (64 bits per voxel)*/
% #define DT_RGB            128    /*A Red-Green-Blue datatype*/
% #define DT_ALL            255    /*Undocumented*/
%
% short int bitpix;    /* Number of bits per pixel; 1, 8, 16, 32, or 64. */
% short int dim_un0;   /* Unused */
%
% float pixdim[];     Parallel array to dim[], giving real world measurements in mm and ms.
%       pixdim[0];    Pixel dimensions?
%       pixdim[1];    Voxel width in mm.
%       pixdim[2];    Voxel height in mm.
%       pixdim[3];    Slice thickness in mm.
%       pixdim[4];    timeslice in ms (ie, TR in fMRI).
%       pixdim[5];    Undocumented.
%       pixdim[6];    Undocumented.
%       pixdim[7];    Undocumented.
%
% float vox_offset;   Byte offset in the .img file at which voxels start. This value can be
%                     negative to specify that the absolute value is applied for every image
%                     in the file.
%
% float roi_scale; Specifies the Region Of Interest scale?
% float funused1; Undocumented.
% float funused2; Undocumented.
%
% float cal_max; Specifies the upper bound of the range of calibration values.
% float cal_min; Specifies the lower bound of the range of calibration values.
%
% int compressed; Undocumented.
% int verified;   Undocumented.
%
% int glmax;    The maximum pixel value for the entire database.
% int glmin;    The minimum pixel value for the entire database.
%
%
                                                                                                                                                                                                                                                                                                                                                                                                                                                               'crc'
      };

    dirname = filename;
    clear filename
    [path, file] = fileparts(dirname);
    for i=1:hdr.nChans
      downscale(i) = 0;
      if ~isempty(strmatch(hdr.label{i}, statuschannel))
        format{i} = 'uint32';
      else
        format{i} = 'int32';
      end
      filename{i} = fullfile(dirname, [file '.' hdr.label{i} '.bin']);
    end

    if ~isfolder(dirname)
      mkdir(dirname);
    end

    % open and write to the output files, one for each selected channel
    fid = zeros(hdr.nChans,1);
    for j=1:hdr.nChans

      if append==false
        fid(j) = fopen_or_error(filename{j}, 'wb', 'ieee-le'); % open the file
        magic = format{j};                               % this used to be the channel name
        magic((end+1):8) = 0;                            % pad with zeros
        magic(8) = downscale(j);                         % number of bits to shift
        fwrite(fid(j), magic(1:8));                      % write the 8-byte file header
      else
        fid(j) = fopen_or_error(filename{j}, 'ab', 'ieee-le');    % open the file for appending
      end % if append

      % convert the data into the correct class
      buf = dat(j,:);
      if ~strcmp(class(buf), format{j})
        switch format{j}
          case 'int16'
            buf = int16(buf);
          case 'int32'
            buf = int32(buf);
          case 'single'
            buf = single(buf);
          case 'double'
            buf = double(buf);
          case 'uint32'
            buf = uint32(buf);
          otherwise
            ft_error('unsupported format conversion');
        end
      end

      % apply the scaling, this corresponds to bit shifting
      buf = buf ./ (2^downscale(j));

      % write the segment of data to the output file
      fwrite(fid(j), buf, format{j}, 'ieee-le');

      fclose(fid(j));
    end % for each channel

  case {'flac' 'm4a' 'mp4' 'oga' 'ogg' 'wav'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This writes data Y to a Windows WAVE file specified by the file name
    % WAVEFILE, with a sample rate of FS Hz and with NBITS number of bits.
    % NBITS must be 8, 16, 24, or 32.  For NBITS < 32, amplitude values
    % outside the range [-1,+1] are clipped
    %
    % Supported extensions for AUDIOWRITE are:
    %   .flac
    % 	.m4a
    % 	.mp4
    % 	.oga
    % 	.ogg
    % 	.wav
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not supported for this data format');
    end

    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end

    if nchans~=1
      ft_error('this format only supports single channel continuous data');
    end

    [p, f, x] = fileparts(filename);
    if isempty(x)
      % append the format as extension
      filename = [filename '.' dataformat];
    end

    audiowrite(filename, dat, hdr.Fs, 'BitsPerSample', nbits);

  case 'plexon_nex'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single or mulitple channel Plexon NEX file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end

    [path, file] = fileparts(filename);
    filename = fullfile(path, [file, '.nex']);
    if nchans~=1
      ft_error('only supported for single-channel data');
    end
    % construct a NEX structure with  the required parts of the header
    nex.hdr.VarHeader.Type       = 5; % continuous
    nex.hdr.VarHeader.Name       = hdr.label{1};
    nex.hdr.VarHeader.WFrequency = hdr.Fs;
    if isfield(hdr, 'FirstTimeStamp')
      nex.hdr.FileHeader.Frequency = hdr.Fs * hdr.TimeStampPerSample;
      nex.var.ts = hdr.FirstTimeStamp;
    else
      ft_warning('no timestamp information available');
      nex.hdr.FileHeader.Frequency  = nan;
      nex.var.ts = nan;
    end
    nex.var.indx = 0;
    nex.var.dat  = dat;

    write_plexon_nex(filename, nex);

    if 0
      % the following code snippet can be used for testing
      [nex2.var, nex2.hdr] = read_plexon_nex(filename, 'channel', 1);
    end

  case 'neuralynx_ncs'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single channel Neuralynx NCS file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end

    if nchans>1
      ft_error('only supported for single-channel data');
    end

    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.ncs']);

    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      % WARNING the AD channel index assumes that the data was read from a DMA or SDMA file
      % the first 17 channels contain status info, this number is zero-offset
      ADCHANNEL  = chanindx - 17 - 1;
      LABEL      = hdr.label{chanindx};
    elseif hdr.nChans==1
      ADCHANNEL  = -1;            % unknown
      LABEL      = hdr.label{1};  % single channel
    else
      ft_error('cannot determine channel label');
    end

    FSAMPLE    = hdr.Fs;
    RECORDNSMP = 512;
    RECORDSIZE = 1044;

    % cut the downsampled LFP data into record-size pieces
    nrecords = ceil(nsamples/RECORDNSMP);
    fprintf('construct ncs with %d records\n', nrecords);

    % construct a ncs structure with all header details and the data in it
    ncs                = [];
    ncs.NumValidSamp   = ones(1,nrecords) * RECORDNSMP;   % except for the last block
    ncs.ChanNumber     = ones(1,nrecords) * ADCHANNEL;
    ncs.SampFreq       = ones(1,nrecords) * FSAMPLE;
    ncs.TimeStamp      = zeros(1,nrecords,'uint64');

    if rem(nsamples, RECORDNSMP)>0
      % the data length is not an integer number of records, pad the last record with zeros
      dat = cat(2, dat, zeros(nchans, nrecords*RECORDNSMP-nsamples));
      ncs.NumValidSamp(end) = rem(nsamples, RECORDNSMP);
    end

    ncs.dat = reshape(dat, RECORDNSMP, nrecords);

    for i=1:nrecords
      % timestamps should be 64 bit unsigned integers
      ncs.TimeStamp(i) = uint64(hdr.FirstTimeStamp) + uint64((i-1)*RECORDNSMP*hdr.TimeStampPerSample);
    end

    % add the elements that will go into the ascii header
    ncs.hdr.CheetahRev            = '4.23.0';
    ncs.hdr.NLX_Base_Class_Type   = 'CscAcqEnt';
    ncs.hdr.NLX_Base_Class_Name   = LABEL;
    ncs.hdr.RecordSize            = RECORDSIZE;
    ncs.hdr.ADChannel             = ADCHANNEL;
    ncs.hdr.SamplingFrequency     = FSAMPLE;

    % write it to a file
    fprintf('writing to %s\n', filename);
    write_neuralynx_ncs(filename, ncs);

    if 0
      % the following code snippet can be used for testing
      ncs2 = read_neuralynx_ncs(filename, 1, inf);
    end

  case 'gdf'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiple channel GDF file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    if ~isempty(chanindx)
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    write_gdf(filename, hdr, dat);

  case 'edf'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiple channel European Data Format file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    if ~isempty(chanindx)
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    write_edf(filename, hdr, dat);

  case 'anywave_ades'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see http://meg.univ-amu.fr/wiki/AnyWave:ADES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    if ~isempty(chanindx)
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label     = hdr.label(chanindx);
      hdr.chantype  = hdr.chantype(chanindx);
      hdr.chanunit  = hdr.chanunit(chanindx);
      hdr.nChans    = length(chanindx);
    end

    dattype = unique(hdr.chantype);
    datunit = cell(size(dattype));
    for i=1:numel(dattype)
      unit = hdr.chanunit(strcmp(hdr.chantype, dattype{i}));
      if ~all(strcmp(unit, unit{1}))
        ft_error('channels of the same type with different units are not supported');
      end
      datunit{i} = unit{1};
    end

    % only change these after checking channel types and units
    chantype = adestype(hdr.chantype);
    dattype  = adestype(dattype);

    % ensure that all channels have the right scaling
    for i=1:size(dat,1)
      switch chantype{i}
        case 'MEG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'T');
        case 'Reference'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'T');
        case 'GRAD'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'T/m');
        case 'EEG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'V');
        case 'SEEG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'V');
        case 'EMG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'V');
        case 'ECG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'V');
        otherwise
          % FIXME I am not sure what scaling to apply
      end
    end

    [p, f, x] = fileparts(filename);
    filename = fullfile(p, f); % without extension
    mat2ades(dat, filename, hdr.Fs, hdr.label, chantype, dattype, datunit);

  otherwise
    ft_error('unsupported data format');
end % switch dataformat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function type = adestype(type)
for i=1:numel(type)
  switch lower(type{i})
    case 'meggrad'
      type{i} = 'MEG'; % this is for CTF and BTi/4D
    case 'megmag'
      type{i} = 'MEG'; % this is for Neuromag and BTi/4D
    case 'megplanar'
      type{i} = 'GRAD'; % this is for Neuromag
    case {'refmag' 'refgrad'}
      type{i} = 'Reference';
    case 'eeg'
      type{i} = 'EEG';
    case {'seeg' 'ecog' 'ieeg'}
      type{i} = 'SEEG'; % all intracranial channels
    case 'ecg'
      type{i} = 'ECG';
    case 'emg'
      type{i} = 'EMG';
    case 'trigger'
      type{i} = 'Trigger';
    case 'source'
      type{i} = 'Source'; % virtual channel
    otherwise
      type{i} = 'Other';
  end
end
                                                                                                                                                                                                                                                                                                                                                                                                                                               \n', BrainStructurelabel{i}, surffile);
      
      % also add metadata to gifti, which avoids wb_view to ask for it
      % interactively upon opening the file
      metadata.name = 'AnatomicalStructurePrimary';
      metadata.value = uppercase2lowercase(BrainStructurelabel{i});
      
      ft_write_headshape(surffile, mesh, 'format', 'gifti', 'metadata', metadata);
    end
    
  else
    mesh.pnt = source.pos;
    mesh.tri = source.tri;
    mesh.unit = source.unit;
    
    [p, f, x] = fileparts(filename);
    filetok = tokenize(f, '.');
    surffile = fullfile(p, [filetok{1} '.surf.gii']);
    ft_write_headshape(surffile, mesh, 'format', 'gifti');
  end
  
end % if writesurface and isfield tri


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to print lists of numbers with appropriate whitespace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = printwithspace(x)
x = x(:)'; % convert to vector
if all(round(x)==x)
  % print as integer value
  s = sprintf('%d ', x);
else
  % print as floating point value
  s = sprintf('%f ', x);
end
s = s(1:end-1);

function s = printwithcomma(x)
x = x(:)'; % convert to vector
if all(round(x)==x)
  % print as integer value
  s = sprintf('%d,', x);
else
  % print as floating point value
  s = sprintf('%f,', x);
end
s = s(1:end-1);

function s = stringpad(s, n)
while length(s)<n
  s = [' ' s];
end

function s = uppercase2lowercase(s)
sel = [0 strfind(s,'_') numel(s)+1];
sout = '';
for m = 1:numel(sel)-1
  sout = [sout, s(sel(m)+1) lower(s((sel(m)+2):(sel(m+1)-1)))];
end
s = sout;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION from roboos/matlab/triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pntR, triR] = remove_vertices(pnt, tri, removepnt)

npnt = size(pnt,1);
ntri = size(tri,1);

if all(removepnt==0 | removepnt==1)
  removepnt = find(removepnt);
end

% remove the vertices and determine the new numbering (indices) in numb
keeppnt = setdiff(1:npnt, removepnt);
numb    = zeros(1,npnt);
numb(keeppnt) = 1:length(keeppnt);

% look for triangles referring to removed vertices
removetri = false(ntri,1);
removetri(ismember(tri(:,1), removepnt)) = true;
removetri(ismember(tri(:,2), removepnt)) = true;
removetri(ismember(tri(:,3), removepnt)) = true;

% remove the vertices and triangles
pntR = pnt(keeppnt, :);
triR = tri(~removetri,:);

% renumber the vertex indices for the triangles
triR = numb(triR);
                                                                                                                                                                                                                                                                                                                                                                                                                                             a precomputed triangulation of some sort
      shape = tmp.bnd;
    elseif isfield(tmp, 'mesh')
      % the variable in the file is most likely a precomputed triangulation of some sort
      shape = tmp.mesh;
    elseif isfield(tmp, 'elec')
      % represent the electrodes as headshape
      tmp.elec        = ft_datatype_sens(tmp.elec);
      shape.fid.pos   = tmp.elec.chanpos;
      shape.fid.label = tmp.elec.label;
    elseif isfield(tmp, 'Vertices')
      % this applies to BrainStorm cortical meshes
      shape.pos = tmp.Vertices;
      % copy some optional fields over with a new name
      shape = copyfields(tmp, shape, {'Faces', 'Curvature', 'SulciMap'});
      shape = renamefields(shape, {'Faces', 'Curvature', 'SulciMap'}, {'tri', 'curv', 'sulc'});
    elseif numel(fieldnames(tmp))==1
      fn = fieldnames(tmp);
      shape = tmp.(fn{1});
      % check that it has vertices and triangles
      assert(isfield(shape, 'pos') && isfield(shape, 'tri'), 'no headshape found in MATLAB file')
    else
      ft_error('no headshape found in MATLAB file');
    end
    
  case {'freesurfer_triangle_binary', 'freesurfer_quadrangle'}
    % the freesurfer toolbox is required for this
    ft_hastoolbox('freesurfer', 1);
    
    [pos, tri] = read_surf(filename);
    
    if min(tri(:)) == 0
      % start counting from 1
      tri = tri + 1;
    end
    shape.pos = pos;
    shape.tri = tri;
    
    % for the left and right
    [path,name,ext] = fileparts(filename);
    
    if strcmp(ext, '.inflated') % does the shift only for inflated surface
      if strcmp(name, 'lh')
        % assume freesurfer inflated mesh in mm, mni space
        % move the mesh a bit to the left, to avoid overlap with the right
        % hemisphere
        shape.pos(:,1) = shape.pos(:,1) - max(shape.pos(:,1)) - 10;
        
      elseif strcmp(name, 'rh')
        % id.
        % move the mesh a bit to the right, to avoid overlap with the left
        % hemisphere
        shape.pos(:,1) = shape.pos(:,1) - min(shape.pos(:,1)) + 10;
      end
    end
    
    if exist(fullfile(path, [name,'.sulc']), 'file'), shape.sulc = read_curv(fullfile(path, [name,'.sulc'])); end
    if exist(fullfile(path, [name,'.curv']), 'file'), shape.curv = read_curv(fullfile(path, [name,'.curv'])); end
    if exist(fullfile(path, [name,'.area']), 'file'), shape.area = read_curv(fullfile(path, [name,'.area'])); end
    if exist(fullfile(path, [name,'.thickness']), 'file'), shape.thickness = read_curv(fullfile(path, [name,'.thickness'])); end
    
  case 'stl'
    [pos, tri, nrm] = read_stl(filename);
    shape.pos = pos;
    shape.tri = tri;
    
  case 'obj'
    ft_hastoolbox('wavefront', 1);
    % Only tested for structure.io .obj thus far
    [vertex, faces, texture, ~] = read_obj_new(filename);
    
    shape.pos   = vertex;
    shape.pos   = shape.pos - repmat(sum(shape.pos)/length(shape.pos),...
        [length(shape.pos),1]); %centering vertices
    shape.tri   = faces(1:end-1,:,:); % remove the last row which is zeros
    
    if hasimage      
      % Refines the mesh and textures to increase resolution of the colormapping
      [shape.pos, shape.tri, texture] = refine(shape.pos, shape.tri,...
          'banks', texture);
      
      picture = imread(image);
      color   = (zeros(length(shape.pos),3));
      for i=1:length(shape.pos)
        color(i,1:3) = picture(floor((1-texture(i,2))*length(picture)),...
            1+floor(texture(i,1)*length(picture)),1:3);
      end
      
      % If color is specified as 0-255 rather than 0-1 correct by dividing
      % by 255
      if range(color(:)) > 1
          color = color./255;
      end
      
      shape.color = color;

    elseif size(vertex,2)==6
      % the vertices also contain RGB colors
      
      color = vertex(:,4:6);
      % If color is specified as 0-255 rather than 0-1 correct by dividing
      % by 255
      if range(color(:)) > 1
          color = color./255;
      end
      
      shape.color = color;
    end
    
  case 'vtk'
    [pos, tri] = read_vtk(filename);
    shape.pos = pos;
    shape.tri = tri;
  
  case 'vtk_xml'
    data = read_vtk_xml(filename);
    shape.orig = data;
    shape.pos  = data.Points;
    if isfield(data, 'Lines')
      shape.line = data.Lines;
    end
  
  case 'mrtrix_tck'
    ft_hastoolbox('mrtrix', 1);
    shape = read_tck(filename);
  
  case 'trackvis_trk'
    shape = read_trk(filename);
  
  case 'off'
    [pos, plc] = read_off(filename);
    shape.pos  = pos;
    shape.tri  = plc;
    
  case 'mne_tri'
    % FIXME this should be implemented, consistent with ft_write_headshape
    keyboard
    
  case 'mne_pos'
    % FIXME this should be implemented, consistent with ft_write_headshape
    keyboard
    
  case 'netmeg'
    hdr = ft_read_header(filename);
    if isfield(hdr.orig, 'headshapedata')
      shape.pos = hdr.orig.Var.headshapedata;
    else
      ft_error('the NetMEG file "%s" does not contain headshape data', filename);
    end
    
  case 'vista'
    ft_hastoolbox('simbio', 1);
    [nodes,elements,labels] = read_vista_mesh(filename);
    shape.pos     = nodes;
    if size(elements,2)==8
      shape.hex     = elements;
    elseif size(elements,2)==4
      shape.tet = elements;
    else
      ft_error('unknown elements format')
    end
    % representation of data is compatible with ft_datatype_parcellation
    shape.tissue = zeros(size(labels));
    numlabels = size(unique(labels),1);
    shape.tissuelabel = {};
    for i = 1:numlabels
      ulabel = unique(labels);
      shape.tissue(labels == ulabel(i)) = i;
      shape.tissuelabel{i} = num2str(ulabel(i));
    end
    
  case 'tet'
    % the toolbox from Gabriel Peyre has a function for this
    ft_hastoolbox('toolbox_graph', 1);
    [vertex, face] = read_tet(filename);
    %     'vertex' is a '3 x nb.vert' array specifying the position of the vertices.
    %     'face' is a '4 x nb.face' array specifying the connectivity of the tet mesh.
    shape.pos = vertex';
    shape.tet = face';
    
  case 'tetgen_ele'
    % reads in the tetgen format and rearranges according to FT conventions
    % tetgen files also return a 'faces' field, which is not used here
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, f); % without the extension
    IMPORT = importdata([filename '.ele'],' ',1);
    shape.tet = IMPORT.data(:,2:5);
    if size(IMPORT.data,2)==6
      labels = IMPORT.data(:,6);
      % representation of tissue type is compatible with ft_datatype_parcellation
      numlabels    = size(unique(labels),1);
      ulabel       = unique(labels);
      shape.tissue = zeros(size(labels));
      shape.tissuelabel = {};
      for i = 1:numlabels
        shape.tissue(labels == ulabel(i)) = i;
        shape.tissuelabel{i} = num2str(ulabel(i));
      end
    end
    IMPORT = importdata([filename '.node'],' ',1);
    shape.pos = IMPORT.data(:,2:4);
    
  case 'brainsuite_dfs'
    % this requires the readdfs function from the BrainSuite MATLAB utilities
    ft_hastoolbox('brainsuite', 1);
    
    dfs = readdfs(filename);
    % these are expressed in MRI dimensions
    shape.pos  = dfs.vertices;
    shape.tri  = dfs.faces;
    shape.unit = 'unkown';
    
    % the filename is something like 2467264.right.mid.cortex.svreg.dfs
    % whereas the corresponding MRI is 2467264.nii and might be gzipped
    [p, f, x] = fileparts(filename);
    while ~isempty(x)
      [junk, f, x] = fileparts(f);
    end
    
    if exist(fullfile(p, [f '.nii']), 'file')
      fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
      mri = ft_read_mri(fullfile(p, [f '.nii']));
      transform = eye(4);
      transform(1:3,4) = mri.transform(1:3,4); % only use the translation
      shape.pos  = ft_warp_apply(transform, shape.pos);
      shape.unit = mri.unit;
    elseif exist(fullfile(p, [f '.nii.gz']), 'file')
      fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
      mri = ft_read_mri(fullfile(p, [f '.nii.gz']));
      transform = eye(4);
      transform(1:3,4) = mri.transform(1:3,4); % only use the translation
      shape.pos  = ft_warp_apply(transform, shape.pos);
      shape.unit = mri.unit;
    else
      ft_warning('could not find accompanying MRI file, returning vertices in voxel coordinates');
    end
    
  case 'brainvisa_mesh'
    % this requires the loadmesh function from the BrainVISA MATLAB utilities
    ft_hastoolbox('brainvisa', 1);
    [shape.pos, shape.tri, shape.nrm] = loadmesh(filename);
    shape.tri = shape.tri + 1; % they should be 1-offset, not 0-offset
    shape.unit = 'unkown';
    
    if exist([filename '.minf'], 'file')
      minffid = fopen_or_error([filename '.minf']);
      hdr=fgetl(minffid);
      tfm_idx = strfind(hdr,'''transformations'':') + 21;
      transform = sscanf(hdr(tfm_idx:end),'%f,',[4 4])';
      fclose(minffid);
      if ~isempty(transform)
        shape.pos = ft_warp_apply(transform, shape.pos);
        shape = rmfield(shape, 'unit'); % it will be determined later on, based on the size
      end
    end
    
    if isempty(transform)
      % the transformation was not present in the minf file, try to get it from the MRI
      
      % the filename is something like subject01_Rwhite_inflated_4d.mesh
      % and it is accompanied by subject01.nii
      [p, f, x] = fileparts(filename);
      f = tokenize(f, '_');
      f = f{1};
      
      if exist(fullfile(p, [f '.nii']), 'file')
        fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
        mri = ft_read_mri(fullfile(p, [f '.nii']));
        shape.pos  = ft_warp_apply(mri.transform, shape.pos);
        shape.unit = mri.unit;
        transform = true; % used for feedback
      elseif exist(fullfile(p, [f '.nii.gz']), 'file')
        fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii.gz']));
        mri = ft_read_mri(fullfile(p, [f '.nii.gz']));
        shape.pos  = ft_warp_apply(mri.transform, shape.pos);
        shape.unit = mri.unit;
        transform = true; % used for feedback
      end
    end
    
    if isempty(transform)
      ft_warning('cound not determine the coordinate transformation, returning vertices in voxel coordinates');
    end
    
  case 'brainvoyager_srf'
    [pos, tri, srf] = read_bv_srf(filename);
    shape.pos = pos;
    shape.tri = tri;
    
    % FIXME add details from srf if possible
    % FIXME do transform
    % FIXME remove vertices that are not in a triangle
    % FIXME add unit
    
  case 'besa_sfp'
    [lab, pos] = read_besa_sfp(filename, 0);
    shape.pos = pos;
    
    % assume that all non-'headshape' points are fiducial markers
    hs = strmatch('headshape', lab);
    lab(hs) = [];
    pos(hs, :) = [];
    shape.fid.label = lab;
    shape.fid.pos = pos;
    
  case 'asa_elc'
    elec = ft_read_sens(filename);
    
    shape.fid.pos   = elec.chanpos;
    shape.fid.label = elec.label;
    
    npos = read_asa(filename, 'NumberHeadShapePoints=', '%d');
    if ~isempty(npos) && npos>0
      origunit = read_asa(filename, 'UnitHeadShapePoints', '%s', 1);
      pos = read_asa(filename, 'HeadShapePoints', '%f', npos, ':');
      pos = ft_scalingfactor(origunit, 'mm')*pos;
      
      shape.pos = pos;
    end
    
  case 'neuromag_mesh'
    fid = fopen_or_error(filename, 'rt');
    npos = fscanf(fid, '%d', 1);
    pos = fscanf(fid, '%f', [6 npos])';
    ntri = fscanf(fid, '%d', 1);
    tri = fscanf(fid, '%d', [3 ntri])';
    fclose(fid);
    
    shape.pos = pos(:,1:3); % vertex positions
    shape.nrm = pos(:,4:6); % vertex normals
    shape.tri = tri;
    
  otherwise
    % try reading it from an electrode of volume conduction model file
    success = false;
    
    if ~success
      % try reading it as electrode positions
      % and treat those as fiducials
      try
        elec = ft_read_sens(filename);
        if ~ft_senstype(elec, 'eeg')
          ft_error('headshape information can not be read from MEG gradiometer file');
        else
          shape.fid.pos   = elec.chanpos;
          shape.fid.label = elec.label;
          success = 1;
        end
      catch
        success = false;
      end % try
    end
    
    if ~success
      % try reading it as volume conductor
      % and treat the skin surface as headshape
      try
        headmodel = ft_read_headmodel(filename);
        if ~ft_headmodeltype(headmodel, 'bem')
          ft_error('skin surface can only be extracted from boundary element model');
        else
          if ~isfield(headmodel, 'skin')
            headmodel.skin = find_outermost_boundary(headmodel.bnd);
          end
          shape.pos = headmodel.bnd(headmodel.skin).pos;
          shape.tri = headmodel.bnd(headmodel.skin).tri; % also return the triangulation
          success = 1;
        end
      catch
        success = false;
      end % try
    end
    
    if ~success
      ft_error('unknown fileformat "%s" for head shape information', fileformat);
    end
end % switch fileformat

if isfield(shape, 'label')
  % ensure that it is a column
  shape.label = shape.label(:);
end

if isfield(shape, 'fid') && isfield(shape.fid, 'label')
  % ensure that it is a column
  shape.fid.label = shape.fid.label(:);
end

% this will add the units to the head shape and optionally convert
if ~isempty(unit)
  shape = ft_convert_units(shape, unit);
else
  try
    % ft_determine_units will fail for triangle-only gifties.
    shape = ft_determine_units(shape);
  end
end

% ensure that vertex positions are given in pos, not in pnt
shape = fixpos(shape);

% ensure that the numerical arrays are represented in double precision and not as integers
shape = ft_struct2double(shape);
end
                                                                                                                                                                     )
            for iSens = length(hdr.label)+1 : orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
              hdr.label{iSens} = ['s2_unknown', num2str(iSens)];
            end
          else
            ft_warning('found more lables in xml.pnsSet than channels in signal 2, thus can not use info in pnsSet, and labeling with s2_eN instead')
            for iSens = orig.signal(1).blockhdr(1).nsignals+1 : orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
              hdr.label{iSens} = ['s2_E', num2str(iSens)];
            end
          end
        else % signal2 is not PIBbox
          ft_warning('creating channel labels for signal 2 on the fly')
          for iSens = 1:orig.signal(2).blockhdr(1).nsignals
            hdr.label{end+1} = ['s2_E', num2str(iSens)];
          end
        end
      elseif length(orig.signal) > 2
        % loop over signals and label channels accordingly
        ft_warning('creating channel labels for signal 2 to signal N on the fly')
        for iSig = 2:length(orig.signal)
          for iSens = 1:orig.signal(iSig).blockhdr(1).nsignals
            if iSig == 1 && iSens == 1
              hdr.label{1} = ['s',num2str(iSig),'_E', num2str(iSens)];
            else
              hdr.label{end+1} = ['s',num2str(iSig),'_E', num2str(iSens)];
            end
          end
        end
      end
    else % no xml.sensorLayout present
      ft_warning('no sensorLayout found in xml files, creating channel labels on the fly')
      for iSig = 1:length(orig.signal)
        for iSens = 1:orig.signal(iSig).blockhdr(1).nsignals
          if iSig == 1 && iSens == 1
            hdr.label{1} = ['s',num2str(iSig),'_E', num2str(iSens)];
          else
            hdr.label{end+1} = ['s',num2str(iSig),'_E', num2str(iSens)];
          end
        end
      end
    end
    
    % check if multiple epochs are present
    if isfield(orig.xml,'epochs')
      % add info to header about which sample correspond to which epochs, becasue this is quite hard for user to get...
      epochdef = zeros(length(orig.xml.epochs),3);
      for iEpoch = 1:length(orig.xml.epochs)
        if iEpoch == 1
          epochdef(iEpoch,1) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./(1000000./hdr.Fs))+1;
          epochdef(iEpoch,2) = round(str2double(orig.xml.epochs(iEpoch).epoch.endTime  )./(1000000./hdr.Fs));
          epochdef(iEpoch,3) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./(1000000./hdr.Fs)); % offset corresponds to timing
        else
          NbSampEpoch = round(str2double(orig.xml.epochs(iEpoch).epoch.endTime)./(1000000./hdr.Fs) - str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./(1000000./hdr.Fs));
          epochdef(iEpoch,1) = epochdef(iEpoch-1,2) + 1;
          epochdef(iEpoch,2) = epochdef(iEpoch-1,2) + NbSampEpoch;
          epochdef(iEpoch,3) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./(1000000./hdr.Fs)); % offset corresponds to timing
        end
      end
      
      if epochdef(end,2) ~= hdr.nSamples
        % check for NS 4.5.4 picosecond timing
        if (epochdef(end,2)/1000) == hdr.nSamples
          for iEpoch=1:size(epochdef,1)
            epochdef(iEpoch,1) = ((epochdef(iEpoch,1)-1)/1000)+1;
            epochdef(iEpoch,2) = epochdef(iEpoch,2)/1000;
            epochdef(iEpoch,3) = epochdef(iEpoch,3)/1000;
          end
          ft_warning('mff apparently generated by NetStation 4.5.4.  Adjusting time scale to microseconds from nanoseconds.');
        else
          ft_error('number of samples in all epochs do not add up to total number of samples')
        end
      end
      
      epochLengths = epochdef(:,2)-epochdef(:,1)+1;
      if ~any(diff(epochLengths))
        hdr.nSamples = epochLengths(1);
        hdr.nTrials  = length(epochLengths);
        
      else
        ft_warning('the data contains multiple epochs with variable length, possibly causing discontinuities in the data')
        % sanity check
        if epochdef(end,2) ~= hdr.nSamples
          % check for NS 4.5.4 picosecond timing
          if (epochdef(end,2)/1000) == hdr.nSamples
            for iEpoch=1:size(epochdef,1)
              epochdef(iEpoch,1)=((epochdef(iEpoch,1)-1)/1000)+1;
              epochdef(iEpoch,2)=epochdef(iEpoch,2)/1000;
              epochdef(iEpoch,3)=epochdef(iEpoch,3)/1000;
            end
            disp('mff apparently generated by NetStation 4.5.4.  Adjusting time scale to microseconds from nanoseconds.');
          else
            ft_error('number of samples in all epochs do not add up to total number of samples')
          end
        end
      end
      orig.epochdef = epochdef;
    end
    hdr.orig = orig;
    
  case 'egi_mff_v2'
    % ensure that the EGI_MFF_V2 toolbox is on the path
    ft_hastoolbox('egi_mff_v2', 1);
    
    %%%%%%%%%%%%%%%%%%%%%%
    %workaround for MATLAB bug resulting in global variables being cleared
    globalTemp=cell(0);
    globalList=whos('global');
    varList=whos;
    for i=1:length(globalList)
      eval(['global ' globalList(i).name ';']);
      eval(['globalTemp{end+1}=' globalList(i).name ';']);
    end
    %%%%%%%%%%%%%%%%%%%%%%
    
    % ensure that the JVM is running and the jar file is on the path
    mff_setup;
    
    %%%%%%%%%%%%%%%%%%%%%%
    %workaround for MATLAB bug resulting in global variables being cleared
    varNames={varList.name};
    for i=1:length(globalList)
      eval(['global ' globalList(i).name ';']);
      eval([globalList(i).name '=globalTemp{i};']);
      if ~any(strcmp(globalList(i).name,varNames)) %was global variable originally out of scope?
        eval(['clear ' globalList(i).name ';']); %clears link to global variable without affecting it
      end
    end
    clear globalTemp globalList varNames varList;
    %%%%%%%%%%%%%%%%%%%%%%
    
    if isunix && filename(1)~=filesep
      % add the full path to the dataset directory
      filename = fullfile(pwd, filename);
    elseif ispc && ~any(strcmp(filename(2),{':','\'}))
      % add the full path, including drive letter or slashes as needed.
      filename = fullfile(pwd, filename);
    end
    
    hdr = read_mff_header(filename);
    
  case {'egi_mff_v3' 'egi_mff'} % this is the default
    ft_hastoolbox('mffmatlabio', 1);
    hdr = mff_fileio_read_header(filename);
    
  case 'fcdc_buffer'
    % read from a networked buffer for realtime analysis
    [host, port] = filetype_check_uri(filename);
    
    if retry
      orig = [];
      while isempty(orig)
        try
          % try reading the header, catch the error and retry
          orig = buffer('get_hdr', [], host, port);
        catch
          ft_warning('could not read header from %s, retrying in 1 second', filename);
          pause(1);
        end
      end % while
    else
      % try reading the header only once, give error if it fails
      orig = buffer('get_hdr', [], host, port);
    end % if retry
    
    % construct the standard header elements
    hdr.Fs          = orig.fsample;
    hdr.nChans      = orig.nchans;
    hdr.nSamples    = orig.nsamples;
    hdr.nSamplesPre = 0;  % since continuous
    hdr.nTrials     = 1;  % since continuous
    hdr.orig        = []; % this will contain the chunks (if present)
    
    % add the contents of attached NEUROMAG_HEADER chunk after decoding to MATLAB structure
    if isfield(orig, 'neuromag_header')
      if isempty(cachechunk)
        % this only needs to be decoded once
        cachechunk = decode_fif(orig);
      end
      
      % convert to FieldTrip format header
      hdr.label       = cachechunk.ch_names(:);
      hdr.nChans      = cachechunk.nchan;
      hdr.Fs          = cachechunk.sfreq;
      
      % add a gradiometer structure for forward and inverse modelling
      try
        [grad, elec] = mne2grad(cachechunk, true, coilaccuracy); % the coordsys is 'dewar'
        if ~isempty(grad)
          hdr.grad = grad;
        end
        if ~isempty(elec)
          hdr.elec = elec;
        end
      catch
        disp(lasterr);
      end
      
      % store the original details
      hdr.orig = cachechunk;
    end
    
    % add the contents of attached CTF_RES4 chunk after decoding to MATLAB structure
    if isfield(orig, 'ctf_res4')
      if isempty(cachechunk)
        % this only needs to be decoded once
        cachechunk = decode_res4(orig.ctf_res4);
      end
      % copy the gradiometer details
      hdr.grad = cachechunk.grad;
      hdr.orig = cachechunk.orig;
      if isfield(orig, 'channel_names')
        % get the same selection of channels from the two chunks
        [selbuf, selres4] = match_str(orig.channel_names, cachechunk.label);
        if length(selres4)<length(orig.channel_names)
          ft_error('the res4 chunk did not contain all channels')
        end
        % copy some of the channel details
        hdr.label     = cachechunk.label(selres4);
        hdr.chantype  = cachechunk.chantype(selres4);
        hdr.chanunit  = cachechunk.chanunit(selres4);
        % add the channel names chunk as well
        hdr.orig.channel_names = orig.channel_names;
      end
      % add the raw chunk as well
      hdr.orig.ctf_res4 = orig.ctf_res4;
    end
    
    % add the contents of attached NIFTI_1 chunk after decoding to MATLAB structure
    if isfield(orig, 'nifti_1')
      hdr.nifti_1 = decode_nifti1(orig.nifti_1);
      % add the raw chunk as well
      hdr.orig.nifti_1 = orig.nifti_1;
    end
    
    % add the contents of attached SiemensAP chunk after decoding to MATLAB structure
    if isfield(orig, 'siemensap') && exist('sap2matlab')==3 % only run this if MEX file is present
      hdr.siemensap = sap2matlab(orig.siemensap);
      % add the raw chunk as well
      hdr.orig.siemensap = orig.siemensap;
    end
    
    if ~isfield(hdr, 'label')
      % prevent overwriting the labels that we might have gotten from a RES4 chunk
      if isfield(orig, 'channel_names')
        hdr.label = orig.channel_names;
      else
        hdr.label = cell(hdr.nChans,1);
        if hdr.nChans < 2000 % don't do this for fMRI etc.
          ft_warning('creating fake channel names');        % give this warning only once
          for i=1:hdr.nChans
            hdr.label{i} = sprintf('%d', i);
          end
        else
          ft_warning('skipping fake channel names');        % give this warning only once
          checkUniqueLabels = false;
        end
      end
    end
    
    if ~isfield(hdr, 'chantype')
      % prevent overwriting the chantypes that we might have gotten from a RES4 chunk
      hdr.chantype = cell(hdr.nChans,1);
      if hdr.nChans < 2000 % don't do this for fMRI etc.
        hdr.chantype = repmat({'unknown'}, 1, hdr.nChans);
      end
    end
    
    if ~isfield(hdr, 'chanunit')
      % prevent overwriting the chanunits that we might have gotten from a RES4 chunk
      hdr.chanunit = cell(hdr.nChans,1);
      if hdr.nChans < 2000 % don't do this for fMRI etc.
        hdr.chanunit = repmat({'unknown'}, 1, hdr.nChans);
      end
    end
    
    hdr.orig.bufsize = orig.bufsize;
    
    
  case 'fcdc_buffer_offline'
    [hdr, nameFlag] = read_buffer_offline_header(headerfile);
    switch nameFlag
      case 0
        % no labels generated (fMRI etc)
        checkUniqueLabels = false; % no need to check these
      case 1
        % has generated fake channels
        % give this warning only once
        ft_warning('creating fake channel names');
        checkUniqueLabels = false; % no need to check these
      case 2
        % got labels from chunk, check those
        checkUniqueLabels = true;
    end
    
  case 'fcdc_matbin'
    % this is multiplexed data in a *.bin file, accompanied by a MATLAB file containing the header
    load(headerfile, 'hdr');
    
  case 'fcdc_mysql'
    % check that the required low-level toolbox is available
    ft_hastoolbox('mysql', 1);
    % read from a MySQL server listening somewhere else on the network
    db_open(filename);
    if db_blob
      hdr = db_select_blob('fieldtrip.header', 'msg', 1);
    else
      hdr = db_select('fieldtrip.header', {'nChans', 'nSamples', 'nSamplesPre', 'Fs', 'label'}, 1);
      hdr.label = mxDeserialize(hdr.label);
    end
    
  case 'gtec_hdf5'
    % check that the required low-level toolbox is available
    ft_hastoolbox('gtec', 1);
    % there is only a precompiled *.p reader that reads the whole file at once
    orig = ghdf5read(filename);
    for i=1:numel(orig.RawData.AcquisitionTaskDescription.ChannelProperties.ChannelProperties)
      lab = orig.RawData.AcquisitionTaskDescription.ChannelProperties.ChannelProperties(i).ChannelName;
      typ = orig.RawData.AcquisitionTaskDescription.ChannelProperties.ChannelProperties(1).ChannelType;
      if isnumeric(lab)
        hdr.label{i} = num2str(lab);
      else
        hdr.label{i} = lab;
      end
      if ischar(typ)
        hdr.chantype{i} = lower(typ);
      else
        hdr.chantype{i} = 'unknown';
      end
    end
    hdr.Fs          = orig.RawData.AcquisitionTaskDescription.SamplingFrequency;
    hdr.nChans      = size(orig.RawData.Samples, 1);
    hdr.nSamples    = size(orig.RawData.Samples, 2);
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1; % assume continuous data, not epoched
    assert(orig.RawData.AcquisitionTaskDescription.NumberOfAcquiredChannels==hdr.nChans, 'inconsistent number of channels');
    % remember the complete data upon request
    if cache
      hdr.orig = orig;
    end
    
  case 'gtec_mat'
    % this is a simple MATLAB format, it contains a log and a names variable
    tmp = load(headerfile);
    log   = tmp.log;
    names = tmp.names;
    
    hdr.label = cellstr(names);
    hdr.nChans = size(log,1);
    hdr.nSamples = size(log,2);
    hdr.nSamplesPre = 0;
    hdr.nTrials = 1; % assume continuous data, not epoched
    
    % compute the sampling frequency from the time channel
    sel = strcmp(hdr.label, 'Time');
    time = log(sel,:);
    
    hdr.Fs = 1./(time(2)-time(1));
    
    % also remember the complete data upon request
    if cache
      hdr.orig.log = log;
      hdr.orig.names = names;
    end
    
  case 'gdf'
    % this requires the biosig toolbox
    ft_hastoolbox('BIOSIG', 1);
    % In the case that the gdf files are written by one of the FieldTrip
    % realtime applications, such as biosig2ft, the gdf recording can be
    % split over multiple 1GB files. The sequence of files is then
    %   filename.gdf   <- this is the one that should be specified as the filename/dataset
    %   filename_1.gdf
    %   filename_2.gdf
    %   ...
    
    [p, f, x] = fileparts(filename);
    if exist(sprintf('%s_%d%s', fullfile(p, f), 1, x), 'file')
      % there are multiple files, count the number of additional files (excluding the first one)
      count = 0;
      while exist(sprintf('%s_%d%s', fullfile(p, f), count+1, x), 'file')
        count = count+1;
      end
      hdr = read_biosig_header(filename);
      for i=1:count
        hdr(i+1) = read_biosig_header(sprintf('%s_%d%s', fullfile(p, f), i, x));
        % do some sanity checks
        if hdr(i+1).nChans~=hdr(1).nChans
          ft_error('multiple GDF files detected that should be appended, but the channel count is inconsistent');
        elseif hdr(i+1).Fs~=hdr(1).Fs
          ft_error('multiple GDF files detected that should be appended, but the sampling frequency is inconsistent');
        elseif ~isequal(hdr(i+1).label, hdr(1).label)
          ft_error('multiple GDF files detected that should be appended, but the channel names are inconsistent');
        end
      end % for count
      % combine all headers into one
      combinedhdr             = [];
      combinedhdr.Fs          = hdr(1).Fs;
      combinedhdr.nChans      = hdr(1).nChans;
      combinedhdr.nSamples    = sum([hdr.nSamples].*[hdr.nTrials]);
      combinedhdr.nSamplesPre = 0;
      combinedhdr.nTrials     = 1;
      combinedhdr.label       = hdr(1).label;
      combinedhdr.orig        = hdr; % include all individual file details
      hdr = combinedhdr;
      
    else
      % there is only a single file
      hdr = read_biosig_header(filename);
      % the GDF format is always continuous
      hdr.nSamples = hdr.nSamples * hdr.nTrials;
      hdr.nTrials = 1;
      hdr.nSamplesPre = 0;
    end % if single or multiple gdf files
    
  case {'homer_nirs'}
    % Homer files are MATLAB files in disguise
    orig = load(filename, '-mat');
    
    hdr.label       = {};
    hdr.nChans      = size(orig.d,2);
    hdr.nSamples    = size(orig.d,1);
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1; % assume continuous data, not epoched
    hdr.Fs          = 1/median(diff(orig.t));
    
    % number of wavelengths times sources times detectors
    assert(numel(orig.SD.Lambda)*orig.SD.nSrcs*orig.SD.nDets >= hdr.nChans);
    
    for i=1:hdr.nChans
      hdr.label{i} = num2str(i);
    end
    
    hdr.chantype = repmat({'nirs'}, hdr.nChans, 1);
    hdr.chanunit = repmat({'unknown'}, hdr.nChans, 1);
    
    % convert the measurement configuration details to an optode structure
    try
    end
    hdr.opto = homer2opto(orig.SD);
    
    % keep the header details
    hdr.orig.SD = orig.SD;
    
  case {'itab_raw' 'itab_mhd'}
    % read the full header information frtom the binary header structure
    header_info = read_itab_mhd(headerfile);
    
    % these are the channels that are visible to FieldTrip
    chansel = 1:header_info.nchan;
    
    % convert the header information into a FieldTrip compatible format
    hdr.nChans      = length(chansel);
    hdr.label       = {header_info.ch(chansel).label};
    hdr.label       = hdr.label(:);  % should be column vector
    hdr.Fs          = header_info.smpfq;
    % it will always be continuous data
    hdr.nSamples    = header_info.ntpdata;
    hdr.nSamplesPre = 0; % it is a single continuous trial
    hdr.nTrials     = 1; % it is a single continuous trial
    % keep the original details AND the list of channels as used by FieldTrip
    hdr.orig         = header_info;
    hdr.orig.chansel = chansel;
    % add the gradiometer definition
    hdr.grad         = itab2grad(header_info);
    
  case 'jaga16'
    % this is hard-coded for the Jinga-Hi JAGA16 system with 16 channels
    packetsize = (4*2 + 6*2 + 16*43*2); % in bytes
    % read the first packet
    fid  = fopen_or_error(filename, 'r');
    buf  = fread(fid, packetsize/2, 'uint16');
    fclose(fid);
    
    if buf(1)==0
      % it does not have timestamps, i.e. it is the raw UDP stream
      packetsize = packetsize - 8; % in bytes
      packet     = jaga16_packet(buf(1:(packetsize/2)), false);
    else
      % each packet starts with a timestamp
      packet = jaga16_packet(buf, true);
    end
    
    % determine the number of packets from the file size
    info     = dir(filename);
    npackets = floor((info.bytes)/packetsize/2);
    
    hdr             = [];
    hdr.Fs          = packet.fsample;
    hdr.nChans      = packet.nchan;
    hdr.nSamples    = 43;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = npackets;
    hdr.label       = cell(hdr.nChans,1);
    hdr.chantype    = cell(hdr.nChans,1);
    hdr.chanunit    = cell(hdr.nChans,1);
    for i=1:hdr.nChans
      hdr.label{i} = sprintf('%d', i);
      hdr.chantype{i} = 'eeg';
      hdr.chanunit{i} = 'uV';
    end
    
    % store some low-level details
    hdr.orig.offset     = 0;
    hdr.orig.packetsize = packetsize;
    hdr.orig.packet     = packet;
    hdr.orig.info      