function A04_spectral_decompose(subject_id, freq_band, freq_spacing_ratio, frac_bandwidth, phase_cutoff_freq)






decimation_factor = 5;  % saved data will be at s_rate divided by this number
data_file_str = [subject_id '.mat'];

%% File paths
helper_function_dir_name = '/_TOOLBOXES';
data_dir_name = ['../' subject_id];
input_dir_name = '03_fulldata';
output_dir_name = '04_spectra';
if(~exist(fullfile(pwd,helper_function_dir_name)))
  fprintf('\nERROR: Helper function directory does not exist.\n');
  return;
end
if(~exist(fullfile(pwd,data_dir_name,input_dir_name,data_file_str)))
  fprintf('\nERROR: Cannot find data at [%s].\n', fullfile(pwd,data_dir_name,input_dir_name,data_file_str));
  return;
end

% Import helper functions to Matlab path
addpath(genpath(fullfile(pwd,helper_function_dir_name)));

% Define path to input data and output data
data_path = fullfile(pwd,data_dir_name,input_dir_name);
%   create directory if needed
output_path = fullfile(pwd,data_dir_name,output_dir_name);
if(exist(output_path, 'dir') == 0)
  mkdir(output_path);
end

%% Determine frequency bands
freq_band_centers(1) = freq_band(1);
while freq_band_centers(end) < freq_band(end)
  freq_band_centers(end+1) = freq_band_centers(end) * (1+freq_spacing_ratio);
end
n_bands = numel(freq_band_centers);

%% Load input data
fprintf('Loading %s\n',fullfile(data_path,data_file_str));
load(fullfile(data_path,data_file_str));
n_channels = header_ecog.n_channels;
n_samples = header_ecog.n_samples;
s_rate = header_ecog.sample_rate;
data_evnt = []; clear data_evnt; clear header_evnt; % This will take up too much space so deleting

%% Change header information to fit decimated data points
decimated_time_indices = 1:decimation_factor:n_samples;
header_ecog.n_samples = numel(decimated_time_indices);
header_ecog.sample_rate_time_domain = s_rate;
header_ecog.sample_rate = s_rate/decimation_factor;

%% Calculate FFT of signal
fprintf('\tTaking FFT ');
fft_data = complex(single(zeros(size(data_ecog))));
for channel_n = 1:n_channels
    fft_data(channel_n,:) = fft(squeeze(data_ecog(channel_n,:)),n_samples);
    fprintf('.');
end
fprintf('\n');

%% Spectrally decompose the data by frequency bins and save
fprintf('\tDecomposing data into %d frequency bands from %.1fHz to %.1fHz\n', n_bands, freq_band_centers(1), freq_band_centers(end));
for band_n = 1:n_bands
    
    % Define frequency band range
    low_freq = freq_band_centers(band_n)*(1.0-frac_bandwidth);
    hig_freq = freq_band_centers(band_n)*(1.0+frac_bandwidth);
    
    % Preallocate matrices
    if(freq_band_centers(band_n) < phase_cutoff_freq)
      phs_data = zeros(n_channels,n_samples);
    end
    amp_data = zeros(n_channels,n_samples);
    
    fprintf('\t%d of %d [%3.1fHz] ', band_n, n_bands, freq_band_centers(band_n));
    
    % Loop through all channels and get amplitude and phase information
    current_freq = freq_band_centers(band_n);
    fprintf('\n');
    parfor channel_n = 1:n_channels
        % Bandpass channel data
        bandpass_data = X_bandpassfft(squeeze(fft_data(channel_n,:)), s_rate, ...
            low_freq, hig_freq, n_samples);
        % Generate amplitude envelope and phase information
        %   Medfilt amplitude before downsampling. It also reduces spikeyness
        amp_data(channel_n,:) = medfilt1(double(abs(hilbert(bandpass_data)).^2),ceil(decimation_factor+1));
        if(current_freq < phase_cutoff_freq)
          phs_data(channel_n,:) = angle(hilbert(bandpass_data))';
        end
        fprintf('\b.\n');
    end
    fprintf('\b');
    
    % Downsample and change data type to single to save storage space
    ecog_data_amp = single(amp_data(:,decimated_time_indices));
    if(freq_band_centers(band_n) < phase_cutoff_freq)
      ecog_data_phs = single(phs_data(:,decimated_time_indices));
    end
    
    % Save data for this frequency band
    header_ecog.frequency_band.band_n = band_n;
    header_ecog.frequency_band.frac_bandwidth = frac_bandwidth;
    header_ecog.frequency_band.center = freq_band_centers(band_n);
    header_ecog.frequency_band.lowedge = low_freq;
    header_ecog.frequency_band.highedge = hig_freq;
    header_ecog.frequency_band.freq_band_centers = freq_band_centers;
    if(freq_band_centers(band_n) < phase_cutoff_freq)
      save(fullfile(output_path,[subject_id '_' num2str(band_n) '.mat']), 'header_ecog', 'ecog_data_amp', 'ecog_data_phs', 'trial_info', 'subject_id', '-v7');
    else
      save(fullfile(output_path,[subject_id '_' num2str(band_n) '.mat']), 'header_ecog', 'ecog_data_amp', 'trial_info', 'subject_id', '-v7');
    end
    fprintf('\n');
    
end
fprintf('\n');










