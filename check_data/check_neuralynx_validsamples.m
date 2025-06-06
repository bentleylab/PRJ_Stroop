function check_neuralynx_validsamples(directory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK_NEURALYNX_VALIDSAMPLES searches and inspects neuralynx (.ncs) files
% for invalid samples in a directory and all of its subdirectories. It then
% prints out in the command window how many 512-sample 'reports' containing 
% invalid samples were detected in each neuralynx recording. 
%
% Use as:
%   check_neuralynx_validsamples(directory)
%
% Arjen Stolk, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

% Colin paths:
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

% Search NLX raw data dirs
if ~exist(directory,'dir')
    switch directory
        case 'IR67'
            directory = '/home/knight/ecog/DATA_FOLDER/Irvine/IR67/Datafiles/IR67_SU/';
        case 'IR69'
            directory = '/home/knight/ecog/DATA_FOLDER/Irvine/IR69/Datafiles/IR69_Task_Data/';
        case 'IR75'
            directory = '/home/knight/ecog/DATA_FOLDER/Irvine/IR75/Data/IR75/';
        case 'IR82'
            directory = '/home/knight/ecog/DATA_FOLDER/Irvine/IR82/Datafiles/SU/';
        case 'IR84'
            directory = '/home/knight/ecog/DATA_FOLDER/Irvine/IR84/Datafiles/SU/';
        case 'IR85'
            directory = '/home/knight/ecog/DATA_FOLDER/Irvine/IR85/Datafiles/IR85_SUA_NLX/';
        case 'IR86'
            directory = '/home/knight/ecog/DATA_FOLDER/Irvine/IR86/Datafiles/';
        case 'IR87'
            directory = '/home/knight/ecog/DATA_FOLDER/Irvine/IR87/Datafiles/';
        otherwise
            error(['Unknown directory is not valid dir or SBJ: ' directory]);
    end
    fprintf('============== %s =========\n',directory);
end

ft_warning off % requires ft toolbox

% total_len = zeros([1 0]);
% invalid_n = zeros([1 0]);
% srate     = zeros([1 0]);
% dset_ix   = 0;
list = dir(directory);
for l = 1:numel(list) % list loop
  
  if (strcmp(list(l).name, '.') || strcmp(list(l).name, '..') || strcmp(list(l).name, '.DS_Store')) % ignore '.' and '..' cases
    continue; % skip this 'file'
  end
  
  full_directory = fullfile(directory, list(l).name);
  if isequal(list(l).isdir, 1) % recurse down
    check_neuralynx_validsamples(full_directory);
  elseif isequal(list(l).isdir, 0) && contains(list(l).name, '.ncs')
%     fprintf('======================================================\n');
%     fprintf('\tStarting to check %s\n',directory);
%     fprintf('======================================================\n');
    try
      ncs = read_neuralynx_ncs(full_directory);
      idx = find(ncs.NumValidSamp(1:end-1)<512); % find 'records' with fewer than 512 valid samples (excluding the final record)
%       dset_ix = dset_ix+1;
      total_len = ncs.NRecords*512;
      invalid_n = zeros(size(idx));
      srate = ncs.hdr.SamplingFrequency;
      if ~isempty(idx)
          for i = 1:numel(idx)
              invalid_n(i) = 512-ncs.NumValidSamp(idx(i));
          end
%         fprintf(['>> ' full_directory ' contains ' num2str(numel(idx)) ' records with invalid samples <<\n']);
%         fprintf('%i / %i datasets are invalid / tested (%.4f%%)\n',...
%             sum(invalid_n>0),numel(total_len),sum(invalid_n>0)/numel(total_len));
%         fprintf('invalid lengths in samples: %.2f mean; %i min - %i max\n',...
%             nanmean(invalid_n),min(invalid_n),max(invalid_n));
        invalid_len = invalid_n./srate;
%         fprintf('invalid lengths in time (s): %.4f mean; %.4f min - %.4f max\n',...
%             nanmean(invalid_len),min(invalid_len),max(invalid_len));
        fprintf('%-10i%-10.4f%-10.4f\n',sum(invalid_n),sum(invalid_len),total_len/srate);
%       else
%         fprintf(['>> ' full_directory ' checked out okay <<\n']);
      end
    catch
      fprintf(['>> was not able to check ' full_directory ' <<\n']);
    end
    
    % catch other files/recordings from that session for same channel
    other_rec = find(~cellfun(@isempty,strfind({list.name},[list(l).name(1:end-4) '_0'])));
    if ~isempty(other_rec)
        for ix = 1:numel(other_rec)
            new_dir = fullfile(directory, list(other_rec(ix)).name);
            try
                ncs = read_neuralynx_ncs(new_dir);
                idx = find(ncs.NumValidSamp(1:end-1)<512); % find 'records' with fewer than 512 valid samples (excluding the final record)
                dset_ix = dset_ix+1;
                total_len(dset_ix) = ncs.NRecords*512;
                invalid_n(dset_ix) = NaN;
                srate(dset_ix) = ncs.hdr.SamplingFrequency;
                if ~isempty(idx)
                    for i = 1:numel(idx)
                        invalid_n(dset_ix) = invalid_n(dset_ix)+ncs.NumValidSamp(idx(i))-512;
                    end
                    fprintf(['>> ' new_dir ' contains ' num2str(numel(idx)) ' records with invalid samples <<\n']);
                    fprintf('%i / %i datasets are invalid / tested (%.4f%%)\n',...
                        sum(invalid_n>0),numel(total_len),sum(invalid_n>0)/numel(total_len));
                    fprintf('invalid lengths in samples: %.2f mean; %i min - %i max\n',...
                        nanmean(invalid_n),min(invalid_n),max(invalid_n));
                    invalid_len = invalid_n/srate;
                    fprintf('invalid lengths in time (s): %.4f mean; %.4f min - %.4f max\n',...
                        nanmean(invalid_len),min(invalid_len),max(invalid_len));
%                 else
%                     fprintf(['>> ' new_dir ' checked out okay <<\n']);
                end
            catch
                fprintf(['>> was not able to check ' new_dir ' <<\n']);
            end
        end
    end
    return; % use only one .ncs files from each directory
  end
  
end % end of list loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ncs] = read_neuralynx_ncs(filename, begrecord, endrecord)

% READ_NEURALYNX_NCS reads a single continuous channel file
%
% Use as
%   [ncs] = read_neuralynx_ncs(filename)
%   [ncs] = read_neuralynx_ncs(filename, begrecord, endrecord)

% Copyright (C) 2005-2007, Robert Oostenveld
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

persistent mexWarning
if isempty(mexWarning)
  mexWarning = false;
end

if nargin<2
  begrecord = 1;
end
if nargin<3
  endrecord = inf;
end

% Using mex files will decrease reading times 2-10-fold on big files
% I had 0.5 sec vs 3.1 sec on 20kHz 600sec file

% determine whether precompiled Nlx2MatCSC from Neualynx is available

isMexv6 = false;
isMexv3 = false;
if ispc
  % first look for Neuralynx version 6 libs
  isMexv6 = ft_hastoolbox('neuralynx_v6', 2); % let's leave warnings for debug
  if ~isMexv6
    % look for Ueli's libs as alternative
    isMexv3 = ft_hastoolbox('neuralynx_v3', 2); % let's leave warnings for debug
  end
elseif ismac || isunix
  % look for Ueli's libs only
  isMexv3 = ft_hastoolbox('neuralynx_v3', 2); % let's leave warnings for debug
end

% if ~isMexv6 && ~mexWarning
%   ft_warning('Reading Neuralynx CSC files is faster if you install the MATLAB importer mex files, see http://neuralynx.com/research_software/file_converters_and_utilities/');
%   mexWarning = true;
% end

if isMexv6 || isMexv3
  % Neuralynx mex files use C-style flags, so let's name them for convinience
  READ_ALL = ones(1,5);
  flags = num2cell(diag(READ_ALL), [1,5]);
  [READ_TST, READ_CHAN, READ_FREQ, READ_VAL, READ_SAMP] = flags{:};
  % the request vector will look like:
  %     >> TST_FLAG+FREQ_FLAG+SAMP_FLAG
  %     >> ans =
  %              1
  %              0
  %              1
  %              0
  %              1
  HEADER_NO  = 0;
  HEADER_YES = 1;
  EXTRACT_RECORD_RANGE = 2;
end

% the file starts with a 16*1024 bytes header in ascii, followed by a number of records
hdr = neuralynx_getheader(filename);
fid = fopen(filename, 'rb', 'ieee-le');

% determine the length of the file
fseek(fid, 0, 'eof');
headersize = 16384;
recordsize = 1044;
NRecords   = floor((ftell(fid) - headersize)/recordsize);

if NRecords>0
  
  % read out part of the dataset to detect whether there were jumps
  NRecords_to_read = min(NRecords, 100); % read out maximum 100 blocks of data
  
  if isMexv6
    [TimeStamp, ChanNumber, SampFreq] = Nlx2MatCSC(filename, READ_TST+READ_CHAN+READ_FREQ, HEADER_NO, EXTRACT_RECORD_RANGE, [1, NRecords_to_read]);
    TimeStamp = uint64(TimeStamp); % to match signature of ft_read_... output, as mex gives us doubles
  elseif isMexv3
    % note that the indexing in the mex file is 0-offset (C++ style) rather than 1-offset (MATLAB style)
    [TimeStamp, ChanNumber, SampFreq] = Nlx2MatCSC_v3(filename, READ_TST+READ_CHAN+READ_FREQ, HEADER_NO, EXTRACT_RECORD_RANGE, [0, NRecords_to_read-1]);
    TimeStamp = uint64(TimeStamp); % to match signature of ft_read_... output, as mex gives us doubles      
  else
    TimeStamp        = zeros(1, NRecords_to_read, 'uint64');
    ChanNumber       = zeros(1, NRecords_to_read);
    SampFreq         = zeros(1, NRecords_to_read);
    for k=1:NRecords_to_read
      
      % set to the correct position
      status = fseek(fid, headersize + (k-1)*recordsize, 'bof');
      if status~=0
        ft_error('cannot jump to the requested record');
      end
      
      % read a single continuous data record
      TimeStamp(k)    = fread(fid,   1, 'uint64=>uint64');
      ChanNumber(k)   = fread(fid,   1, 'int32');
      SampFreq(k)     = fread(fid,   1, 'int32');
    end
  end % if isMexv6
  
  % explicitly sort the timestamps to deal with negative timestamp jumps that can occur
  ts1 = TimeStamp(1);
  dts = double(TimeStamp - TimeStamp(1));
  dts = unique(dts);
  dts = sort(dts);
  TimeStamp = uint64(dts) + ts1;
  
  % for this block of data: automatically detect the gaps;
  % there's a gap if no round off error of the sampling frequency could
  % explain the jump (which is always > one block)
  Fs       = mode(double(SampFreq));
  if abs(Fs/hdr.SamplingFrequency-1)>0.01
    ft_warning('the sampling frequency as read out from the header equals %2.2f and differs from the mode sampling frequency as read out from the data %2.2f\n', ...
      hdr.SamplingFrequency, Fs);
    
    % check which one was correct
    d = double(TimeStamp(2:end)-TimeStamp(1:end-1));
    fsEst = 1e6./mode(d);
    indx = nearest([Fs hdr.SamplingFrequency], fsEst);
    if indx==1
      ft_warning('correcting the header frequency from %2.2f to %2.2f', hdr.SamplingFrequency, Fs);
      hdr.SamplingFrequency = Fs;
    end
  end
  
  % detect the number of timestamps per block while avoiding influencce of gaps
  d = double(TimeStamp(2:end)-TimeStamp(1:end-1));
  maxJump  = ceil(10^6./(Fs-1))*512;
  gapCorrectedTimeStampPerSample =  nanmean(d(d<maxJump))/512;
  
  % read the timestamp from the first and last record
  if (ispc), fclose(fid); end
  ts1 = neuralynx_timestamp(filename, 1);
  tsE = neuralynx_timestamp(filename, inf);
  if (ispc), fid = fopen(filename, 'rb', 'ieee-le'); end
  
  hdr.FirstTimeStamp  = ts1;
  hdr.LastTimeStamp   = tsE;
  
  % compare whether there's at least a block missing
  minJump = min(d);
  ts_range_predicted = (NRecords-1)*512*gapCorrectedTimeStampPerSample;
  ts_range_observed  = double(tsE-ts1);
  if abs(ts_range_predicted-ts_range_observed)>minJump
    ft_warning('discontinuous recording, predicted number of timestamps and observed number of timestamps differ by %2.2f \n Please consult the wiki on http://www.fieldtriptoolbox.org/getting_started/neuralynx?&#discontinuous_recordings',...
      abs(ts_range_predicted-ts_range_observed) );
  end
  
else
  hdr.FirstTimeStamp = nan;
  hdr.LastTimeStamp  = nan;
end

if begrecord==0 && endrecord==0
  % only read the header
elseif begrecord<1
  ft_error('cannot read before the first record');
elseif begrecord>NRecords
  ft_error('cannot read beyond the last record')
elseif endrecord>NRecords
  endrecord = NRecords;
end

if begrecord>=1 && endrecord>=begrecord
  % leave numrecord information here for proper synchronisation
  numrecord    = (endrecord-begrecord+1);
  if isMexv6
    % ft_warning('Reading with neuralynx_v6 mex files');
    [TimeStamp, ChanNumber, SampFreq, NumValidSamp, Samp] = Nlx2MatCSC(filename, READ_ALL, HEADER_NO, EXTRACT_RECORD_RANGE, [begrecord, endrecord]);
    TimeStamp = uint64(TimeStamp); % to match signature of ft_read_... output
  elseif isMexv3
    % ft_warning('Reading with neuralynx_v3 mex files');
    % note that the indexing in the mex file is 0-offset (C++ style) rather than 1-offset (MATLAB style)
    [TimeStamp, ChanNumber, SampFreq, NumValidSamp, Samp] = Nlx2MatCSC_v3(filename, READ_ALL, HEADER_NO, EXTRACT_RECORD_RANGE, [begrecord-1, endrecord-1]);
    TimeStamp = uint64(TimeStamp); % to match signature of ft_read_... output
  else
    % ft_warning('Reading with native MATLAB code');
    % rewind to the first record to be read
    status = fseek(fid, headersize + (begrecord-1)*recordsize, 'bof');
    if status~=0
      ft_error('cannot jump to the requested record');
    end
    
    TimeStamp    = zeros(1,numrecord,'uint64');
    ChanNumber   = zeros(1,numrecord);
    SampFreq     = zeros(1,numrecord);
    NumValidSamp = zeros(1,numrecord);
    Samp         = zeros(512,numrecord);  % this allows easy reshaping into a 1xNsamples vector
    
    for k=1:numrecord
      % read a single continuous data record
      TimeStamp(k)    = fread(fid,   1, 'uint64=>uint64');
      ChanNumber(k)   = fread(fid,   1, 'int32');
      SampFreq(k)     = fread(fid,   1, 'int32');
      NumValidSamp(k) = fread(fid,   1, 'int32');
      Samp(:,k)       = fread(fid, 512, 'int16');
      % mark the invalid samples
      Samp((NumValidSamp+1):end,k) = nan;
    end
    
  end % if isMexv6
  
  ts1 = TimeStamp(1);
  dts = double(TimeStamp-ts1); % no problem with doubles here as numbers are small
  
  [val,indx] = sort(dts);
  [A,I] = unique(val); % consider only the unique values
  indx = indx(I);
  
  % store the record data in the output structure
  ncs.TimeStamp    = uint64(dts(indx)) + ts1; % convert back to original class
  ncs.ChanNumber   = ChanNumber(indx);
  ncs.SampFreq     = SampFreq(indx);
  ncs.NumValidSamp = NumValidSamp(indx);
  % apply the scaling factor from ADBitVolts and convert to uV
  ncs.dat          = Samp(:,indx) * hdr.ADBitVolts * 1e6;
end
fclose(fid);

% store the header info in the output structure
ncs.NRecords = NRecords;
ncs.hdr      = hdr;


function [hdr] = neuralynx_getheader(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading the 16384 byte header from any Neuralynx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2007, Robert Oostenveld
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

fid     = fopen(filename, 'rb', 'ieee-le');
buf     = fread(fid, 16*1024, 'uint8=>char');
fclose(fid);

buf     = buf(:)';
nl      = find(buf==10);    % determine the new-lines
cr      = find(buf==13);    % determine the carriage-returns
begchar = [1 nl(1:(end-1))];
endchar = nl - 1;
num     = length(nl);

hdr        = [];
hdr.Header = buf; % remember the full header in its original format

for i=1:num
  line = fliplr(deblank(fliplr(deblank(char(buf(begchar(i):endchar(i)))))));
  if numel(line)==0
    % line is empty
    continue
  elseif line(1)=='#'
    % line contains a comment
    continue
  else
    % strip the '-' sign
    while line(1)=='-'
      line = line(2:end);
    end
    % replace tabs with spaces
    line(find(line==9)) = ' ';
    % cut into pieces
    item = strread(line, '%s');
    if length(item)==2
      key = item{1};
      val = item{2};
      if any(val(1)=='-01234567989')
        % try to convert to number
        val = str2num(val);
        if isempty(val)
          % revert to the original text
          val = item{2};
        end
      end
      % remove unuseable characters from the variable name (key)
      key = key(key~=':');
      % assign the value to the header structure
      indx =  strfind(key,char(181)); % avoid problems with this field in the new Neuralynx header      
      if ~isempty(indx)
        key(indx) = 'm';
      end
      hdr = setfield(hdr, key, val);
    end
  end
end


function [t] = neuralynx_timestamp(filename, num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading a single timestamp of a single channel Neuralynx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2007, Robert Oostenveld
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

headersize = 16384;
switch ft_filetype(filename)
  case 'neuralynx_ncs'
    recordsize = 1044;  % in bytes
  case 'neuralynx_nse'
    recordsize = 112;   % in bytes
  case 'neuralynx_nst'
    recordsize = 304;   % in bytes
  case 'neuralynx_nts'
    recordsize = 8;     % in bytes
  case 'neuralynx_ntt'
    recordsize = 304;   % in bytes
end

fid = fopen(filename, 'rb', 'ieee-le');

if (ispc)
  % this is to fix a bug in the windwos version which does not want to do uint64=>uint64
  % however this code will fail if the MSB is set (only likely in very long recordings)
  if ~isinf(num)
    % read the timestamp of the indicated record
    fseek(fid, headersize + (num-1)*recordsize, 'bof');
    t = fread(fid, 1, 'uint64=>integer*8');
  else
    % read the timestamp of the last record
    fseek(fid, -recordsize, 'eof');
    t = fread(fid, 1, 'uint64=>integer*8');
  end
  t = uint64(t);
else
  if ~isinf(num)
    % read the timestamp of the indicated record
    fseek(fid, headersize + (num-1)*recordsize, 'bof');
    t = fread(fid, 1, 'uint64=>uint64');
  else
    % read the timestamp of the last record
    fseek(fid, -recordsize, 'eof');
    t = fread(fid, 1, 'uint64=>uint64');
  end
end

fclose(fid);
