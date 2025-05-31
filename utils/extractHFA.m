function hfa = extractHFA(roi,bl_start,bl_end,niter,trial_info,alignment) 

%% Extract High-Frequency Activity

% 1. Extract high gamma in 10Hz bins 
% 2. Baseline correction using a bootstrap z-score procedure and average
% across bins
% 3. roi is a fieldtrip data structure 

%% 

ft_defaults;

%% Features

hg_vec   = 75:10:145;      % center frequencies
win      = 5;              % freq window around center

if bl_start ~= 0 && bl_end ~= 0

    %% Cut into Trials
        trial_lim_s_pad = [-2.5 5]; % Add 1.5s buffers at beginning and end of trial
    
        stim_events = trial_info.word_onset;
    
        roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
            round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));
    
    %% Get HFA
    
    % assign memory
    if strcmpi(alignment,'S')
        start_trial = dsearchn(roi_trl.time{1,1}',-1);
        end_trial = dsearchn(roi_trl.time{1,1}',0);

        hgtrace = zeros([length(hg_vec) size(roi_trl.trial,2) length(roi_trl.label) length(roi_trl.time{1,1}(start_trial:end_trial))]);
    else
        start_trial = dsearchn(roi_trl.time{1,1}',-1);
        end_trial = dsearchn(roi_trl.time{1,1}',3.5);
        
        hgtrace = zeros([length(hg_vec) size(roi_trl.trial,2) length(roi_trl.label) length(roi_trl.time{1,1}(start_trial:end_trial))]);

    end
    
    for h = 1:length(hg_vec)

        % filter the data in the range first
        cfg              = [];
        cfg.bpfilter     = 'yes';
        cfg.bpfreq       = [hg_vec(h)-win hg_vec(h)+win];
        cfg.hilbert      = 'abs';
        filt             = ft_preprocessing(cfg, roi_trl);

        if strcmpi(alignment,'S')
            
            % zbaseline correct
            filt             = zbaseline_HFA(filt, bl_start, bl_end, niter);
    
            % Trim buffer periods off
            cfg_trim = [];
            cfg_trim.latency = [-1 0];
            filt = ft_selectdata(cfg_trim,filt);
            t = filt.time{1}';

            hgtrace(h,:,:,:) = permute(cat(3,filt.trial{:}),[3 1 2]);
            clear filt
        else

            % zbaseline correct
            filt             = zbaseline_HFA(filt, bl_start, bl_end, niter);

            % Trim buffer periods off
            cfg_trim = [];
            cfg_trim.latency = [-1 3.5];
            filt = ft_selectdata(cfg_trim,filt);
            t = filt.time{1}';
            
    
            hgtrace(h,:,:,:) = permute(cat(3,filt.trial{:}),[3 1 2]);
            clear filt
        end

    end
    
    if strcmpi(alignment,'S')
        hfa.ztrial = squeeze(mean(hgtrace));
        hfa.time = t;
    else
        hfa.ztrial = squeeze(mean(hgtrace));
        hfa.time = t;
        hfa_r = hfa;
        hfa_r.ztrial = NaN(size(hfa_r.ztrial,1),size(hfa_r.ztrial,2),length(-0.75:1/roi_trl.fsample:0.75));
        for triali = 1:size(hfa_r.ztrial,1)
            rt_idx = nearest(hfa.time,trial_info.response_time(triali));
            hfa_r.ztrial(triali,:,:) = hfa.ztrial(triali,:,rt_idx-375:rt_idx+375);
        end
        hfa_r.time = -0.75:1/roi_trl.fsample:0.75;
        hfa = hfa_r;
        clear hfa_r;
    end
else
    %% Cut into Trials around response for baseline period/ITI
    buffers = [-2.0 2.5]; % Add 1.5s buffers at beginning and end of trial

    stim_events = trial_info.resp_onset;

    roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
        round([buffers(1) buffers(2)]*roi.fsample));

    %% Get HFA
    
    % assign memory
    % hgtrace = zeros([length(hg_vec) size(roi_trl.trial,2) length(roi_trl.label) length(roi_trl.time{1,1}) - (3*roi_trl.fsample)]);
    
    for h = 1:length(hg_vec)

        % filter the data in the range first
        cfg              = [];
        cfg.bpfilter     = 'yes';
        cfg.bpfreq       = [hg_vec(h)-win hg_vec(h)+win];
        cfg.hilbert      = 'abs';
        filt             = ft_preprocessing(cfg, roi_trl);

        cfg.padding = length(roi.trial{1,1})/roi.fsample + 3;
        cfg.padtype = 'mirror';
        cfg.paddir  = 'both';
        filtref     = ft_preprocessing(cfg, roi);

        % Trim buffer periods off
        cfg_trim = [];
        cfg_trim.latency = [-0.5 1.0];
        filt = ft_selectdata(cfg_trim,filt);
        t = filt.time;

        % zbaseline correct
        filt             = zfullrecording_HFA(filt,filtref,niter);
        
        hgtrace(h,:,:,:) = filt.ztrial;
        clear filt
    end
    
    hfa.ztrial = squeeze(mean(hgtrace));
    hfa.time = t;
    
end
    
    
end
