for channel_n = 1:size(data_preproc.trial{1},1)
    [fft_data,freqs] = pwelch(data_preproc.trial{1}(channel_n,:),2048,0,2048,data_preproc.fsample);
    [fft_data2,freqs2] = pwelch(data_compare.trial{1}(channel_n,:),2048,0,2048,data_compare.fsample);
    loglog(freqs,fft_data,'k'); hold on;
    loglog(freqs2,fft_data2,'r');
    xlim([1 350]);
    ax = gca;
    ax.XTick = [25 30 60 100 120 180 200 240 300 360 420];
    title(['Channel ' num2str(channel_n) ' [' data_preproc.label{channel_n} ']']);
    legend('preproc','compare');
    pause;
end


%%
data_notch_reref = [];
data_notch_reref.fsample = data_ft_notch.fsample;
data_notch_reref.time{1} = data_ft_notch.time{1};
data_notch_reref.trial{1} = NaN([size(data_ft_notch.trial{1},1)-1,size(data_ft_notch.trial{1},2)]);
reref_labels = {};
for ix = 1:size(data_notch_reref.trial{1},1)
    % Medial-Lateral
    data_notch_reref.trial{1}(ix,:) = data_ft_notch.trial{1}(ix,:)-data_ft_notch.trial{1}(ix+1,:);
    data_notch_reref.label{ix} = strcat(data_ft_notch.label{ix},'-',data_ft_notch.label{ix+1});
end

%% next
figure;subplot(2,3,1);
sig_ix = 4;
sig = [];
sig.fsample = raw.header_ecog.sample_rate;
sig.label = raw.header_ecog.channel_labels(sig_ix);
sig.time{1} = linspace(0, size(raw.data_ecog,2), size(raw.data_ecog,2));
sig.trial{1} = raw.data_ecog(sig_ix,:);
plot(sig.trial{1});
title('raw signal');
subplot(2,3,4);
fn_plot_PSD_compare(sig.trial{1},[sig.label],sig.fsample);

subplot(2,3,2);
sig_rs = ft_resampledata(cfg_resamp,sig);
plot(sig_rs.trial{1});
title('sig_resamp');
subplot(2,3,5);
fn_plot_PSD_compare(sig_rs.trial{1},[sig_rs.label],sig_rs.fsample);

subplot(2,3,3);
cfg_notch300 = cfg_notch; cfg_notch300.dftfreq = [300];
sig_rs_n300 = ft_preprocessing(cfg_notch300, sig_rs);
sig_rs_n = ft_preprocessing(cfg_notch, sig_rs);
plot(sig_rs_n300.trial{1},'r'); hold on;
plot(sig_rs_n.trial{1},'k');
legend('n300','n');
subplot(2,3,6);
fn_plot_PSD_compare(vertcat(sig_rs_n300.trial{1},sig_rs_n.trial{1}),[sig_rs_n300.label sig_rs_n.label],sig_rs.fsample);
title('sig_rs_n300 vs. sig_rs_nfull');

simple = raw.data_ecog(sig_ix,:);
[p,q] = rat(1000/raw.header_ecog.sample_rate);
simple_rs = resample(simple, p, q);
simple_rs_notch = cleanline2(simple_rs, 1000, ...
    'LineFrequencies',  line_noise_freqs, ...   % default is 60, 120
    'ScanForLines',     1, ...                  % finds exact line freq around given value
    'LineAlpha',        0.01, ...               % default = 0.01
    'Bandwidth',        3, ...                  % default = 1
    'ChanCompIndices',  1:1,...%header_ecog.n_channels, ...
    'SlidingWinLength', 4.0, ...                % in sec, default = 4
    'SlidingWinStep',   3.0, ...                % in sec, default = 4 (no overlap)
    'SmoothingFactor',  100, ...                % default=100; 1=linear, Inf=no smooth
    'PaddingFactor',    1, ...                  % default = 2; CWH good with this
    'ComputeSpectralPower', 0, ...              % default = 1/True; might be nice
    'PlotFigures', 1, ...                       % maybe alternative to CompSpecPow
    'VerboseOutput',    0);
simple_notch = cleanline2(simple, header_ecog.sample_rate, ...
    'LineFrequencies',  line_noise_freqs, ...   % default is 60, 120
    'ScanForLines',     1, ...                  % finds exact line freq around given value
    'LineAlpha',        0.01, ...               % default = 0.01
    'Bandwidth',        3, ...                  % default = 1
    'ChanCompIndices',  1:1,...%header_ecog.n_channels, ...
    'SlidingWinLength', 4.0, ...                % in sec, default = 4
    'SlidingWinStep',   3.0, ...                % in sec, default = 4 (no overlap)
    'SmoothingFactor',  100, ...                % default=100; 1=linear, Inf=no smooth
    'PaddingFactor',    1, ...                  % default = 2; CWH good with this
    'ComputeSpectralPower', 0, ...              % default = 1/True; might be nice
    'PlotFigures', 1, ...                       % maybe alternative to CompSpecPow
    'VerboseOutput',    0);