%% SBJ01_import_data.m
% Extract data with fieldtrip and save out by data type
function SBJ01_import_data(SBJ,resamp_freq)

addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/'));
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

if isempty(resamp_freq)
    resamp_yn = 'no';
else
    if ischar(resamp_freq)  % for SGE input
        resamp_freq = str2num(resamp_freq);
    end
    resamp_yn   = 'yes';
end

%% Load and preprocess the data
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Process channel labels
bad_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.bad);
eeg_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eeg);
photod_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.photod);
mic_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.mic);

if strcmp(SBJ_vars.raw_file(end-2:end),'mat')
    orig = load(SBJ_vars.dirs.raw_filename);
    % Load Neural Data
    cfg = [];
    cfg.channel = {'all',bad_ch_neg{:},eeg_ch_neg{:},photod_ch_neg{:},mic_ch_neg{:}};
    data = ft_selectdata(cfg,orig.data);
    % EEG data
    if ~isempty(SBJ_vars.ch_lab.eeg)
        cfg.channel = SBJ_vars.ch_lab.eeg;
        eeg = ft_selectdata(cfg,orig.data);
    end
    % Load event data
    cfg.channel = {SBJ_vars.ch_lab.photod{:},SBJ_vars.ch_lab.mic{:}};
    evnt = ft_selectdata(cfg,orig.data);
else
    % Load Neural Data
    cfg            = [];
    cfg.dataset    = SBJ_vars.dirs.raw_filename;
    cfg.continuous = 'yes';
    cfg.channel    = {'all',bad_ch_neg{:},eeg_ch_neg{:},photod_ch_neg{:},mic_ch_neg{:}};
    data = ft_preprocessing(cfg);
    
    % Load EEG data
    if ~isempty(SBJ_vars.ch_lab.eeg)
        cfg.channel = SBJ_vars.ch_lab.eeg;
        eeg = ft_preprocessing(cfg);
    end
    
    % Load event data
    cfg.channel = {SBJ_vars.ch_lab.photod{:},SBJ_vars.ch_lab.mic{:}};
    evnt = ft_preprocessing(cfg);
end

%% Cut to analysis_time
if ~isempty(SBJ_vars.analysis_time)
    for epoch_ix = 1:length(SBJ_vars.analysis_time)
        epoch_len{epoch_ix} = diff(SBJ_vars.analysis_time{epoch_ix});
        cfg_trim = [];
        cfg_trim.latency = SBJ_vars.analysis_time{epoch_ix};
        
        data_pieces{epoch_ix} = ft_selectdata(cfg_trim, data);
        evnt_pieces{epoch_ix} = ft_selectdata(cfg_trim, evnt);
        if ~isempty(SBJ_vars.ch_lab.eeg)
            eeg_pieces{epoch_ix}  = ft_selectdata(cfg_trim, eeg);
        end
    end
    % Stitch them back together
    data = data_pieces{1};
    data.time{1} = data.time{1}-SBJ_vars.analysis_time{1}(1);
    evnt = evnt_pieces{1};
    evnt.time{1} = evnt.time{1}-SBJ_vars.analysis_time{1}(1);
    if ~isempty(SBJ_vars.ch_lab.eeg)
        eeg = eeg_pieces{1};
        eeg.time{1} = eeg.time{1}-SBJ_vars.analysis_time{1}(1);
    end
    if length(SBJ_vars.analysis_time)>1
        for epoch_ix = 2:length(SBJ_vars.analysis_time)
            data.trial{1} = horzcat(data.trial{1},data_pieces{epoch_ix}.trial{1});
            data.time{1} = horzcat(data.time{1},data_pieces{epoch_ix}.time{1}-...
                SBJ_vars.analysis_time{epoch_ix}(1)+data.time{1}(end)+data.time{1}(2));
            
            evnt.trial{1} = horzcat(evnt.trial{1},evnt_pieces{epoch_ix}.trial{1});
            evnt.time{1} = horzcat(evnt.time{1},evnt_pieces{epoch_ix}.time{1}-...
                SBJ_vars.analysis_time{epoch_ix}(1)+evnt.time{1}(end)+evnt.time{1}(2));
            
            if ~isempty(SBJ_vars.ch_lab.eeg)
                eeg.trial{1} = horzcat(eeg.trial{1},eeg_pieces{epoch_ix}.trial{1});
                eeg.time{1} = horzcat(eeg.time{1},eeg_pieces{epoch_ix}.time{1}-...
                    SBJ_vars.analysis_time{epoch_ix}(1)+eeg.time{1}(end)+eeg.time{1}(2));
            end
        end
    end
end

%% Process Microphone data
cfg = [];
cfg.channel = SBJ_vars.ch_lab.mic;
mic_data = ft_selectdata(cfg,evnt);
mic_data = mic_data.trial{1};
%rescale to prevent clipping, add 0.05 fudge factor
mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);

%% Resample data
if strcmp(resamp_yn,'yes') && (data.fsample ~= resamp_freq)
    cfg = [];
    cfg.resamplefs = resamp_freq;
    cfg.detrend = 'no';
    data = ft_resampledata(cfg, data);
    if ~isempty(SBJ_vars.ch_lab.eeg)
        eeg = ft_resampledata(cfg, eeg);
    end
end

%% Save data
nrl_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_',num2str(resamp_freq),'hz.mat');
save(nrl_out_filename, '-v7.3', 'data');

if ~isempty(SBJ_vars.ch_lab.eeg)
    eeg_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_eeg_',num2str(resamp_freq),'hz.mat');
    save(eeg_out_filename, '-v7.3', 'eeg');
end

evnt_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_evnt.mat');
save(evnt_out_filename, '-v7.3', 'evnt');

% Save Microphone data separately as .wav for listening
mic_data_filename = strcat(SBJ_vars.dirs.import,SBJ,'_mic_recording.wav');
audiowrite(mic_data_filename,mic_data_rescale,evnt.fsample);

end