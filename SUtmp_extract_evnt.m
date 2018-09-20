% function create_nk_photodiode(subjectm)
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(genpath('/Users/colinhoy/Code/Apps/wave_clus/'));
addpath(genpath('/Users/colinhoy/Code/Apps/UR_NLX2MAT_releaseDec2015/'));
addpath(ft_dir);
ft_defaults

%% Read photodiode
photo_label = 'Photo1_0010';
inverted = 1;
nlx_dir = [root_dir 'PRJ_Stroop/data/IR69/00_raw/SU_R1_2018-02-07_13-17-17/'];
photo = ft_read_neuralynx_interp({[nlx_dir 'photo/' photo_label '.ncs']});
photo.label = {'Photo1'};
if inverted
    photo.trial{1} = photo.trial{1}*-1;
end

%% Read, downsample, and save mic
mic_label = 'Mic1_0010';
mic = ft_read_neuralynx_interp({[nlx_dir 'mic/' mic_label '.ncs']});
mic.label = {'Mic1'};
if proc_vars.resample_mic
    cfg_dsmp = [];
    cfg_dsmp.resamplefs = photo.fsample;
    mic_dsmp = ft_resampledata(cfg_dsmp,mic);
end

%% Align Photodiode and Mic
fname = {[nlx_dir 'photo/' photo_label '.ncs'], [nlx_dir 'mic/' mic_label '.ncs']};
% get the minimum and maximum timestamp across all channels
for i = 1:numel(fname)
  ts    = ft_read_data(fname{i}, 'timestamp', true);  
  mn(i) = ts(1);
  mx(i) = ts(end); 
  ts1   = ts(1);
  md(i) = mode(diff(double(ts-ts1)));

  % get the minimum and maximum across all channels
  if i>1 && mn(i)<min_all
    min_all = mn(i);
  else
    min_all = mn(i);
  end
  
  if i>1 && mx(i)>max_all
    max_all = mx(i);
  else
    max_all = mx(i); 
  end
end

% take the mode of the modes and construct the interpolation axis
% the interpolation axis should be casted in doubles
mode_dts  = mode(md);
rng       = double(max_all-min_all); % this is small num, can be double
offset    = double(mn-min_all); % store the offset per channel
offsetmx  = double(max_all-mx);
tsinterp  = [0:mode_dts:rng]; % the timestamp interpolation axis


%% Process Microphone data
mic_data = mic.trial{1};
%rescale to prevent clipping, add 0.05 fudge factor
mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);

%% Concatenate photo and mic
cfga = [];
evnt = ft_appenddata(cfga,photo,mic_dsmp);

%% Save data out
evnt_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_evnt',block_suffix,'.mat');
save(evnt_out_filename, '-v7.3', 'evnt');

% Save Microphone data separately as .wav for listening
mic_data_filename = strcat(SBJ_vars.dirs.import,SBJ,'_mic_recording',block_suffix,'.wav');
audiowrite(mic_data_filename,mic_data_rescale,evnt.fsample);

% nlx_hdr = ft_read_header(nlx_photo_dir);
% nlx_photo_ix = find(~cellfun(@isempty, regexp(nlx_hdr.label, photo_label))== 1);
% nlx_photo = ft_read_data(nlx_photo_dir, 'header', nlx_hdr, 'chanidx', 1);
% t_photo = linspace(double(nlx_hdr.FirstTimeStamp),double(nlx_hdr.LastTimeStamp),nlx_hdr.nSamples);
% t_photo_s = t_photo/1000000;
% t_photo_s = t_photo_s-t_photo_s(1);
%
% % View photodiode in chunks to find the task data
% chunk_sec = 20*60;  % 20 min chunks
% n_chunks = round(nlx_hdr.nSamples/nlx_hdr.Fs/chunk_sec);
% chunk_sz = 20*60*nlx_hdr.Fs;
% chunks = zeros([n_chunks 2]);
% chunks(1,:) = [1 chunk_sz];
% for chunk_ix = 2:n_chunks
%     chunks(chunk_ix,:) = chunks(chunk_ix-1,:)+chunk_sz;
%     plot(t_photo_s(chunks(chunk_ix,:)), nlx_photo(chunks(chunk_ix,:)));
%     pause;
% end

% %% Load spike sorting clusters and align
% load([nlx_dir 'micro/times_' micro_label '.mat']);
% hdr = ft_read_header([nlx_dir 'micro/' micro_label '.ncs']);
% % Convert cluster_class spike times to ms relative to start of file
% cluster_class(:,2) = cluster_class(:,2)-double(hdr.FirstTimeStamp)/1000;    % adjust hdr from 1000 kHz to 1 kHz

% --------------------------------------------------------
% Synchronizes Neuralynx with Nihon Kohden recording, and creates a
% photodiode channel for the latter based on the first
% Tested on IR69 and IR75
%
% Arjen Stolk, 2018
% astolk@berkeley.edu
% --------------------------------------------------------


% subject-specific settings
chan              = 'RAM2'; % channel used for sync'ing
photo             = 'Photo1'; % Neuralynx photodiode
dset_nl           = '/Users/arjsto/Projects/Ecog/data/IR69/datafiles/2018-02-10_09-24-14'; % Neuralynx dataset
tcgver            = 1; % 1: TCG, 2: TCG kids, 3: dual with IR76 (IR75 only)
d                 = 1; % dataset

% header information
eval(subjectm);

% load nihon kohden
dset_nk           = [subj(tcgver).datadir subj(tcgver).eegfile{d}];
hdr_nk            = ft_read_header(dset_nk);
chanindx_nk       = find(~cellfun(@isempty, regexp(hdr_nk.label, chan))== 1);
coi_nk            = ft_read_data(dset_nk, 'header', hdr_nk, 'chanindx', chanindx_nk);

% load neuralynx
hdr_nl            = ft_read_header(dset_nl);
chanindx_nl       = find(~cellfun(@isempty, regexp(hdr_nl.label, chan))== 1);
coi_nl            = ft_read_data(dset_nl, 'header', hdr_nl, 'chanindx', chanindx_nl);
coi_nl            = coi_nl.*-1; % there might be a flip of sign (should end up with a positive corr)
photoindx_nl      = find(~cellfun(@isempty, regexp(hdr_nl.label, photo))== 1);
if ~isempty(photoindx_nl) % if photodiode is stored in the same folder
    hdr_nl_photo      = hdr_nl;
    photo_nl          = ft_read_data(dset_nl, 'header', hdr_nl, 'chanindx', photoindx_nl);
else % if photodiode is stored separately because of a separate sample frequency (folder_photo)
    hdr_nl_photo      = ft_read_header([dset_nl '_photo']);
    photoindx_nl      = find(~cellfun(@isempty, regexp(hdr_nl_photo.label, photo))== 1);
  photo_nl          = ft_read_data([dset_nl '_photo'], 'header', hdr_nl_photo, 'chanindx', photoindx_nl);
  % additionally correct for timestamp offsets between the two (note,
  % neuralynx thinks at a resolution of 1000kHz, meaning there are 125
  % samples in-between two consecutive samples of a timeseries sampled with
  % 8kHz)
  t_chan = linspace(double(hdr_nl.FirstTimeStamp),double(hdr_nl.LastTimeStamp),hdr_nl.nSamples);
  t_photo = linspace(double(hdr_nl_photo.FirstTimeStamp),double(hdr_nl_photo.LastTimeStamp),hdr_nl_photo.nSamples);
  t_shared = intersect(t_chan,t_photo);
  i_chan = ismember(t_chan, t_shared);
  i_photo = ismember(t_photo, t_shared);
  coi_nl = coi_nl(find(i_chan==1,1,'first'):find(i_chan==1,1,'last')); % coi_nl(i_chan), but allowing anti-aliasing when downsampling still
  photo_nl = photo_nl(find(i_photo==1,1,'first'):find(i_photo==1,1,'last')); % photo_nl(i_photo)
end
photo_nl          = photo_nl.*-1; % there might be a flip of sign (should end up with a positive corr)

% downsample neuralynx to nihon kohden samplerate
if hdr_nl.Fs>hdr_nk.Fs
  fprintf('downsampling Neuralynx from %d to %d Hz\n', hdr_nl.Fs, hdr_nk.Fs)
  coi_nl_ds         = ft_preproc_resample(coi_nl, hdr_nl.Fs, hdr_nk.Fs, 'resample');
else
  coi_nl_ds         = coi_nl;
end
if hdr_nl_photo.Fs>hdr_nk.Fs
  fprintf('downsampling Neuralynx photodiode from %d to %d Hz\n', hdr_nl_photo.Fs, hdr_nk.Fs)
  photo_nl_ds       = ft_preproc_resample(photo_nl, hdr_nl_photo.Fs, hdr_nk.Fs, 'resample');
else
  photo_nl_ds       = photo_nl;
end

% remove extreme values
coi_nk((coi_nk>median(coi_nk)+4*std(coi_nk))|(coi_nk<median(coi_nk)-4*std(coi_nk))) = median(coi_nk);
coi_nl_ds((coi_nl_ds>median(coi_nl_ds)+4*std(coi_nl_ds))|(coi_nl_ds<median(coi_nl_ds)-4*std(coi_nl_ds))) = median(coi_nl_ds);

% synchronize nihon kohden and neuralynx timeseries
[c, lags]         = xcov(coi_nk', coi_nl_ds');
[~, idx]          = max(c - smooth(c, hdr_nk.Fs*10)); % find the sharp peak

% doublecheck
figure; subplot(2,1,1);
hold on; plot(lags,c);
hold on; plot(lags(idx),c(idx),'k*');
ylabel('correlation');
xlabel('lag');
subplot(2,1,2);
t = 1:numel(coi_nk); 
hold on; plot(t, zscore(coi_nk));
t2 = lags(idx):lags(idx)+numel(coi_nl_ds)-1;
hold on; plot(t2, zscore(coi_nl_ds)+10);
t3 = lags(idx):lags(idx)+numel(photo_nl_ds)-1; % ignore any offset between photo and chan
hold on; plot(t3, zscore(photo_nl_ds)+20);
legend('NK', 'NL', 'NL photo');
print([subj(1).datadir 'datafiles/sync_nk-nl_' subjectm(9:end) '_' num2str(tcgver) '_' num2str(d)], '-dpdf');

% create nihon kohden photodiode channel
trig_ts = ones(1,numel(coi_nk)).*median(photo_nl_ds);
trig_ts(t3(t3>0 & t3<numel(trig_ts))) = photo_nl_ds(t3>0 & t3<numel(trig_ts));
save([subj(tcgver).datadir subj(tcgver).eegfile{d}(1:end-4) '_photo_' num2str(tcgver) '_' num2str(d) '.mat'], 'trig_ts');
