%% SBJ01_import_data.m
% Extract data with fieldtrip and save out by data type
function SU00_extract_data(SBJ,pipeline_id)
error('didnt write pretty much any of this! short cutting for the CPPC2018 abstract, will come back (hoepfully)');

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';apps_dir=[root_dir 'Apps/'];
else root_dir='/Volumes/hoycw_clust/';apps_dir='/Users/colinhoy/Code/Apps/';end
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(genpath([apps_dir 'wave_clus/']));
addpath([apps_dir 'fieldtrip/']);
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m'];
eval(proc_vars_cmd);

%% Process channel labels
if isfield(SBJ_vars.ch_lab,'prefix')
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
        SBJ_vars.ch_lab.bad{bad_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.bad{bad_ix}];
    end
    for eeg_ix = 1:numel(SBJ_vars.ch_lab.eeg)
        SBJ_vars.ch_lab.eeg{eeg_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.eeg{eeg_ix}];
    end
    for eog_ix = 1:numel(SBJ_vars.ch_lab.eog)
        SBJ_vars.ch_lab.eog{eog_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.eog{eog_ix}];
    end
    SBJ_vars.ch_lab.photod = {[SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.photod{1}]};
    SBJ_vars.ch_lab.mic    = {[SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.mic{1}]};
end
if isfield(SBJ_vars.ch_lab,'suffix')
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
        SBJ_vars.ch_lab.bad{bad_ix} = [SBJ_vars.ch_lab.bad{bad_ix} SBJ_vars.ch_lab.suffix];
    end
    for eeg_ix = 1:numel(SBJ_vars.ch_lab.eeg)
        SBJ_vars.ch_lab.eeg{eeg_ix} = [SBJ_vars.ch_lab.eeg{eeg_ix} SBJ_vars.ch_lab.suffix];
    end
    for eog_ix = 1:numel(SBJ_vars.ch_lab.eog)
        SBJ_vars.ch_lab.eog{eog_ix} = [SBJ_vars.ch_lab.eog{eog_ix} SBJ_vars.ch_lab.suffix];
    end
    SBJ_vars.ch_lab.photod = {[SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.suffix]};
    SBJ_vars.ch_lab.mic    = {[SBJ_vars.ch_lab.mic{1} SBJ_vars.ch_lab.suffix]};
end
bad_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.bad);
eeg_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eeg);
eog_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eog);
photod_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.photod);
mic_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.mic);

for b_ix = 1:numel(SBJ_vars.block_name)
    %% Load Data
    if numel(SBJ_vars.raw_file)>1
        block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
    else
        block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
    end
!!! start here    if strcmp(SBJ_vars.raw_file{b_ix}(end-2:end),'mat')
        orig = load(SBJ_vars.dirs.raw_filename{b_ix});
        % Load Neural Data
        cfg = [];
        cfg.channel = {'all',bad_ch_neg{:},eeg_ch_neg{:},eog_ch_neg{:},photod_ch_neg{:},mic_ch_neg{:}};
        data = ft_selectdata(cfg,orig.data);
        % EEG data
        if ~isempty(SBJ_vars.ch_lab.eeg)
            cfg.channel = SBJ_vars.ch_lab.eeg;
            eeg = ft_selectdata(cfg,orig.data);
        end
        % EOG data
        if ~isempty(SBJ_vars.ch_lab.eog)
            cfg.channel = SBJ_vars.ch_lab.eog;
            eog = ft_selectdata(cfg,orig.data);
        end
        % Load event data
        cfg.channel = {SBJ_vars.ch_lab.photod{:},SBJ_vars.ch_lab.mic{:}};
        evnt = ft_selectdata(cfg,orig.data);
    else
        % Load Neural Data
        cfg            = [];
        cfg.dataset    = SBJ_vars.dirs.raw_filename{b_ix};
        cfg.continuous = 'yes';
        cfg.channel    = {'all',bad_ch_neg{:},eeg_ch_neg{:},eog_ch_neg{:},photod_ch_neg{:},mic_ch_neg{:}};
        data = ft_preprocessing(cfg);
        
        % Load EEG data
        if ~isempty(SBJ_vars.ch_lab.eeg)
            cfg.channel = SBJ_vars.ch_lab.eeg;
            eeg = ft_preprocessing(cfg);
        end
        % Load EOG data
        if ~isempty(SBJ_vars.ch_lab.eog)
            cfg.channel = SBJ_vars.ch_lab.eog;
            eog = ft_preprocessing(cfg);
        end
        
        % Load event data
        cfg.channel = {SBJ_vars.ch_lab.photod{:},SBJ_vars.ch_lab.mic{:}};
        evnt = ft_preprocessing(cfg);
    end
    % Toss 'EDF Annotations' Channel (if it exists)
    if any(strcmp(data.label,'EDF Annotations'))
        cfg_s = [];
        cfg_s.channel = {'all','-EDF Annotations'};
        data = ft_selectdata(cfg_s,data);
    end
    
    %% Cut to analysis_time
    if ~isempty(SBJ_vars.analysis_time{b_ix})
        for epoch_ix = 1:length(SBJ_vars.analysis_time{b_ix})
            epoch_len{epoch_ix} = diff(SBJ_vars.analysis_time{b_ix}{epoch_ix});
            cfg_trim = [];
            cfg_trim.latency = SBJ_vars.analysis_time{b_ix}{epoch_ix};
            
            data_pieces{epoch_ix} = ft_selectdata(cfg_trim, data);
            evnt_pieces{epoch_ix} = ft_selectdata(cfg_trim, evnt);
            if ~isempty(SBJ_vars.ch_lab.eeg)
                eeg_pieces{epoch_ix}  = ft_selectdata(cfg_trim, eeg);
            end
            if ~isempty(SBJ_vars.ch_lab.eog)
                eog_pieces{epoch_ix}  = ft_selectdata(cfg_trim, eog);
            end
        end
        % Stitch them back together
        data = data_pieces{1};
        data.time{1} = data.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);
        evnt = evnt_pieces{1};
        evnt.time{1} = evnt.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);
        if ~isempty(SBJ_vars.ch_lab.eeg)
            eeg = eeg_pieces{1};
            eeg.time{1} = eeg.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);
        end
        if ~isempty(SBJ_vars.ch_lab.eog)
            eog = eog_pieces{1};
            eog.time{1} = eog.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);
        end
        if length(SBJ_vars.analysis_time{b_ix})>1
            for epoch_ix = 2:length(SBJ_vars.analysis_time{b_ix})
                data.trial{1} = horzcat(data.trial{1},data_pieces{epoch_ix}.trial{1});
                data.time{1} = horzcat(data.time{1},data_pieces{epoch_ix}.time{1}-...
                    SBJ_vars.analysis_time{b_ix}{epoch_ix}(1)+data.time{1}(end)+data.time{1}(2));
                
                evnt.trial{1} = horzcat(evnt.trial{1},evnt_pieces{epoch_ix}.trial{1});
                evnt.time{1} = horzcat(evnt.time{1},evnt_pieces{epoch_ix}.time{1}-...
                    SBJ_vars.analysis_time{b_ix}{epoch_ix}(1)+evnt.time{1}(end)+evnt.time{1}(2));
                
                if ~isempty(SBJ_vars.ch_lab.eeg)
                    eeg.trial{1} = horzcat(eeg.trial{1},eeg_pieces{epoch_ix}.trial{1});
                    eeg.time{1} = horzcat(eeg.time{1},eeg_pieces{epoch_ix}.time{1}-...
                        SBJ_vars.analysis_time{b_ix}{epoch_ix}(1)+eeg.time{1}(end)+eeg.time{1}(2));
                end
                if ~isempty(SBJ_vars.ch_lab.eog)
                    eog.trial{1} = horzcat(eog.trial{1},eog_pieces{epoch_ix}.trial{1});
                    eog.time{1} = horzcat(eog.time{1},eog_pieces{epoch_ix}.time{1}-...
                        SBJ_vars.analysis_time{b_ix}{epoch_ix}(1)+eog.time{1}(end)+eog.time{1}(2));
                end
            end
        end
    end
    
    %% Resample data
    if strcmp(proc_vars.resample_yn,'yes') && (data.fsample ~= proc_vars.resample_freq)
        cfg = [];
        cfg.resamplefs = proc_vars.resample_freq;
        cfg.detrend = 'no';
        data = ft_resampledata(cfg, data);
        if ~isempty(SBJ_vars.ch_lab.eeg)
            eeg = ft_resampledata(cfg, eeg);
        end
        if ~isempty(SBJ_vars.ch_lab.eog)
            eog = ft_resampledata(cfg, eog);
        end
    end
    
    %% Process Microphone data
    cfg = [];
    cfg.channel = SBJ_vars.ch_lab.mic;
    mic_data = ft_selectdata(cfg,evnt);
    mic_data = mic_data.trial{1};
    %rescale to prevent clipping, add 0.05 fudge factor
    mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);
    
    %% Final Channel Label Corrections
    % Strip Pre/Suffix if Necessary
    for ch_ix = 1:numel(data.label)
        if isfield(SBJ_vars.ch_lab,'prefix')
            data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.prefix,'');
        end
        if isfield(SBJ_vars.ch_lab,'suffix')
            data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.suffix,'');
        end
    end
    if ~isempty(SBJ_vars.ch_lab.eeg)
        for eeg_ix = 1:numel(eeg.label)
            if isfield(SBJ_vars.ch_lab,'prefix')
                eeg.label{eeg_ix} = strrep(eeg.label{eeg_ix},SBJ_vars.ch_lab.prefix,'');
            end
            if isfield(SBJ_vars.ch_lab,'suffix')
                eeg.label{eeg_ix} = strrep(eeg.label{eeg_ix},SBJ_vars.ch_lab.suffix,'');
            end
        end
    end
    if ~isempty(SBJ_vars.ch_lab.eog)
        for eog_ix = 1:numel(eog.label)
            if isfield(SBJ_vars.ch_lab,'prefix')
                eog.label{eog_ix} = strrep(eog.label{eog_ix},SBJ_vars.ch_lab.prefix,'');
            end
            if isfield(SBJ_vars.ch_lab,'suffix')
                eog.label{eog_ix} = strrep(eog.label{eog_ix},SBJ_vars.ch_lab.suffix,'');
            end
        end
    end
    if isfield(SBJ_vars.ch_lab,'prefix')
        evnt.label{1} = strrep(evnt.label{1},SBJ_vars.ch_lab.prefix,'');
        evnt.label{2} = strrep(evnt.label{2},SBJ_vars.ch_lab.prefix,'');
    end
    if isfield(SBJ_vars.ch_lab,'suffix')
        evnt.label{1} = strrep(evnt.label{1},SBJ_vars.ch_lab.suffix,'');
        evnt.label{2} = strrep(evnt.label{2},SBJ_vars.ch_lab.suffix,'');
    end
    
    % Fix any mislabeled channels
    if isfield(SBJ_vars.ch_lab,'mislabel')
        for ch_ix = 1:numel(SBJ_vars.ch_lab.mislabel)
            % Future edit: search for the bad label across data, eeg, evnt
            data.label(strcmp(data.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1))) = SBJ_vars.ch_lab.mislabel{ch_ix}(2);
        end
    end
    
    %% Save data
    nrl_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_',num2str(proc_vars.resample_freq),'hz',block_suffix,'.mat');
    save(nrl_out_filename, '-v7.3', 'data');
    
    if ~isempty(SBJ_vars.ch_lab.eeg)
        eeg_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_eeg_',num2str(proc_vars.resample_freq),'hz',block_suffix,'.mat');
        save(eeg_out_filename, '-v7.3', 'eeg');
    end
    if ~isempty(SBJ_vars.ch_lab.eog)
        eog_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_eog_',num2str(proc_vars.resample_freq),'hz',block_suffix,'.mat');
        save(eog_out_filename, '-v7.3', 'eog');
    end
    
    evnt_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_evnt',block_suffix,'.mat');
    save(evnt_out_filename, '-v7.3', 'evnt');
    
    % Save Microphone data separately as .wav for listening
    mic_data_filename = strcat(SBJ_vars.dirs.import,SBJ,'_mic_recording',block_suffix,'.wav');
    audiowrite(mic_data_filename,mic_data_rescale,evnt.fsample);
    
    clear data evnt eeg eog mic_data mic_data_rescale
end

end
