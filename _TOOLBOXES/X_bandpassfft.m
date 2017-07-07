function smoothdata = X_bandpassfft(data, s_rate, lowcut, highcut, nfft)
%
% Inputs:
%   data      = [n_samples] FFTed signal of data to filter
%   s_rate    = data sampling rate (Hz)
%   lowcut    = low-edge frequency in pass band (Hz)
%   highcut   = high-edge frequency in pass band (Hz)
%   nfft      = Number of FFT points
%
%
% Outputs:
%    smoothdata = filtered signal
%
%% Example:
% s_rate = 500;
% n_samples = 1000;
% sig_freq = 20;
% lowcut = 18;
% highcut = 22;
% noise_level = 3;
% rand_phase = rand*pi;
% data = sin(linspace(0-rand_phase,sig_freq*2*pi*(n_samples/s_rate)-rand_phase,n_samples))*1.0 + rand(1,n_samples)*noise_level-(noise_level/2);
% X_bandpassfft(data, s_rate, lowcut, highcut);
%%


plot_it = 0;
plot_freqs = [0 80];

% Get size of data
n_samples = numel(data);

% Frequency vector for plotting
freq_vector = s_rate/2*linspace(0,1,1+nfft/2);
freq_vector = [freq_vector(1:end) freq_vector((end-1):-1:2)];

if( (highcut-lowcut) < (freq_vector(9)-freq_vector(7)) )
    fprintf('\n\nERROR IN X_bandpassfft().  Difference between lowcut and highcut is too small.\n');
    fprintf('\tlowcut: %3.3f\n', lowcut);
    fprintf('\thighcut: %3.3f\n', highcut);
    fprintf('\tdifference: %1.3f\n', highcut-lowcut);
    fprintf('\tNFFT (%4.0f) makes frequency resolution %1.3f,\n', nfft, (freq_vector(9)-freq_vector(8)));
    fprintf('\tso difference between highcut and lowcut should be at least %1.3f\n\n', (freq_vector(9)-freq_vector(7)));
end

% find closest freq in fft decomposition to low and high cutoffs
[~, idxl_1] = min(abs(freq_vector(1:floor(length(freq_vector)/2))-lowcut));
[~, idxh_1] = min(abs(freq_vector(1:floor(length(freq_vector)/2))-highcut));
[~, idxl_2] = min(abs(freq_vector(ceil(length(freq_vector)/2):end)-lowcut));  idxl_2 = idxl_2+ceil(length(freq_vector)/2)-1;
[~, idxh_2] = min(abs(freq_vector(ceil(length(freq_vector)/2):end)-highcut)); idxh_2 = idxh_2+ceil(length(freq_vector)/2)-1;

if(plot_it == 1)
    figure('Position', [40 80 1000 800]);
end

% filter the data
% ---------------
if(plot_it == 1)
    % Plot original signal
    subplot(5,1,[1 2]); hold on;
    padded_smoothdata = 2*real(ifft(squeeze(data)));
    smoothdata = padded_smoothdata(1:n_samples);
    plot(linspace(0,n_samples/s_rate,n_samples),smoothdata,'b');
    clear padded_smoothdata;
    clear smoothdata;
    xlim([0 n_samples/s_rate]);
    title('Original = blue, Filtered = red');
end

if(plot_it == 1)
    % Plot FFT signal and filter limits
    subplot(5,1,[4 5]);
    plot(freq_vector, abs(data), 'b');
    %set(gca, 'YScale', 'log');
    xlabel('Frequency [Hz]');
    ylabel('Signal power [dB]');
    hold on;
    ymin = min(abs(data));
    ymax = max(abs(data));
    plot([freq_vector(idxl_1) freq_vector(idxl_1)], [ymin ymin+(ymax-ymin)/2], 'k', 'LineWidth', 2);
    plot([freq_vector(idxh_1) freq_vector(idxh_1)], [ymin ymin+(ymax-ymin)/2], 'k', 'LineWidth', 2);
    plot([freq_vector(idxl_2) freq_vector(idxl_2)], [ymin+(ymax-ymin)/2 ymax], 'g', 'LineWidth', 2);
    plot([freq_vector(idxh_2) freq_vector(idxh_2)], [ymin+(ymax-ymin)/2 ymax], 'g', 'LineWidth', 2);
    xlim([plot_freqs(1) plot_freqs(2)]);
end

% Generate gaussian mask
mask_gaussianity = 1.0;
%freq_indices = 1:length(freq_vector);
center_index = round((idxh_1+idxl_1)/2);
width = idxh_1-center_index;
g_1 = ngaussian(1:length(freq_vector),center_index,width,mask_gaussianity);
center_index = round((idxh_2+idxl_2)/2);
width = idxh_2-center_index;
g_2 = ngaussian(1:length(freq_vector),center_index,width,mask_gaussianity);
gaussian_mask = g_1 + g_2;

if(plot_it == 1)
    % Plot FFT mask
    subplot(5,1,3);
    plot(freq_vector, gaussian_mask);
    xlim([plot_freqs(1) plot_freqs(2)]);
end

% Mask the signal
%  If the signal has an odd length, the mask might need to be adjusted
if(length(data) == length(gaussian_mask))
  masked_fft_signal = data .* gaussian_mask;
elseif(length(data) > length(gaussian_mask))
  masked_fft_signal = data .* [gaussian_mask 0];
elseif(length(data) < length(gaussian_mask))
  masked_fft_signal = data .* gaussian_mask(1:(end-1));
end

if(plot_it == 1)
    % Plot masked FFT signal
    y_min = min(abs(data));
    y_max = max(abs(data));
    subplot(5,1,[4 5]);
    plot(freq_vector, abs(masked_fft_signal), 'r');
    xlim([plot_freqs(1) plot_freqs(2)]);
    ylim([y_min y_max]);
    set(gca, 'YScale', 'log');
end

padded_smoothdata = 2*real(ifft(masked_fft_signal));
smoothdata = padded_smoothdata(1:n_samples);

if(plot_it == 1)
    % Plot filtered signal
    subplot(5,1,[1 2]); hold on;
    plot(linspace(0,n_samples/s_rate,n_samples),2*real(ifft(masked_fft_signal)),'r');
    xlim([0 n_samples/s_rate]);
    pause;
    close all;
end

function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Gaussian when n=1, becomes more rectangular as n increases.
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
g = exp(-((x-pos)./(0.6006.*wid)) .^(2*round(n)));







