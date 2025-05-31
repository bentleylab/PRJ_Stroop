% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];
proc_id = 'main_ft';


%
cnt = 1;
for subj = 1:numel(files)

    load(['/Users/anaskhan/Desktop/PRJ_Stroop/results/newROIs/Power/HFA/Resp/' files{subj}])
    load(['/Users/anaskhan/Desktop/PRJ_Stroop/results/newROIs/Power/Resp/' files{subj}])

    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])
    
    % dlpfcChans = elec.label(ismember(elec.ROI,'dlPFC'));
    
    % Conflict encoding
    dlpfcChans = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));
    % History encoding
    % dlpfcChans = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,3}.sig_chans));

    if isempty(dlpfcChans)
        continue
    else
       
       
        [lpfcDesign,~] = fn_create_LMM_design(trial_info,1,1,'power');
    
        rowsToremove_power = strcmpi(lpfcDesign.PreviousType,'None');
    
        hfa.ztrial(rowsToremove_power,:,:) = [];
        lpfcDesign(rowsToremove_power,:) = [];

        tempdata = hfa.ztrial(strcmpi(lpfcDesign.CurrentType,'inc'),:,:);

        for itrial = 1:sum(strcmpi(lpfcDesign.CurrentType,'inc'))
            if size(tempdata(itrial,:,:),2) == 1
                hfa_temp.trial{1,itrial} = squeeze(tempdata(itrial,:,:))';
            else
                hfa_temp.trial{1,itrial} = squeeze(tempdata(itrial,:,:));
            end
            hfa_temp.time{1,itrial} = hfa.newtime;
        end

        hfa_temp.label = hfa.label;

        cfgs = [];
        cfgs.channel = dlpfcChans;
        cfgs.latency = [0 0.5];
        hfa_temp = ft_selectdata(cfgs,hfa_temp);

        hfa = hfa_temp;

        clear hfa_temp tempdata

        hfa_data = permute(cat(3,hfa.trial{:}),[3 1 2]);
        hfa_time = hfa.time{1,1};
        %%
        cfgs = [];
        cfgs.channel = dlpfcChans;
        roi = ft_selectdata(cfgs,data);
        origsamp = roi.fsample;
        %% Bring all subjects' fs down to 500
        if origsamp ~= 500
            cfg_resamp = [];
            cfg_resamp.resamplefs = 500;
        
            roi = ft_resampledata(cfg_resamp,roi);
            trial_info.sample_rate = roi.fsample;
                if (origsamp/roi.fsample) == 2 
                    trial_info.word_onset = round(trial_info.word_onset./(origsamp/roi.fsample));
                    trial_info.resp_onset = round(trial_info.resp_onset./(origsamp/roi.fsample));
                    trial_info.sample_rate = roi.fsample;
                    if isfield(trial_info,'run_len')
                        trial_info.run_len = round(trial_info.run_len./(origsamp/roi.fsample));
                    end
                else
                    error('Discrepancies in sampling rate. Check again.')
                end
        end
        %% Cut into Trials
        trial_lim_s_pad = [-1 1.5];
    
        events = trial_info.resp_onset;
    
        roi_trl = fn_ft_cut_trials_equal_len(roi,events,trial_info.condition_n',...
            round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));

        
        for ch = 1:size(hfa_data,2)
            for itrial = 1:size(hfa_data,1)

                [~,locs] = findpeaks(squeeze(hfa_data(itrial,ch,:)));
                temp_peaks = hfa_time(locs);
                temp_peaks(temp_peaks < 0.1) =[];

                peak_times{ch,itrial} = temp_peaks;
            end
        end

        roi_trl.trial(rowsToremove_power) = [];
        roi_trl.time(rowsToremove_power) = [];
        roi_trl.trialinfo(rowsToremove_power) = [];
        roi_trl.sampleinfo(rowsToremove_power,:) = [];

        roi_trl.trial = roi_trl.trial(strcmpi(lpfcDesign.CurrentType,'inc'));
        roi_trl.time = roi_trl.time(strcmpi(lpfcDesign.CurrentType,'inc'));
        roi_trl.trialinfo = roi_trl.trialinfo(strcmpi(lpfcDesign.CurrentType,'inc'));
        roi_trl.sampleinfo = roi_trl.sampleinfo(strcmpi(lpfcDesign.CurrentType,'inc'),:);

        tempdata = permute(cat(3,roi_trl.trial{:}), [3 1 2]);
        half_window = round(0.5*roi.fsample);
        pta_time = roi_trl.time{1,1}(nearest(roi_trl.time{1,1},-0.5):nearest(roi_trl.time{1,1},0.5));

        for ch = 1:size(hfa_data,2)
            for itrial = 1:size(hfa_data,1)
                for time = 1:numel(peak_times{ch,itrial})
                    tidx = nearest(roi_trl.time{1,1},peak_times{ch,itrial}(time));
                    pta(ch,itrial,:) = tempdata(itrial,ch,tidx-half_window:tidx+half_window);
                end
            end
        end
        PTA_data{cnt} = pta;

        % PSD of HFA
        % cfg = [];
        % cfg.foilim     = [1 45];
        % cfg.method     = 'irasa';
        % cfg.keeptrials = 'yes';
        % cfg.pad        = 4;
        % cfg.output     = 'fractal';
        % 
        % fractal = ft_freqanalysis(cfg,hfa);
        % 
        % cfg.output = 'original';
        % 
        % original = ft_freqanalysis(cfg,hfa);
        % 
        % cfg               = [];
        % cfg.parameter     = 'powspctrm';
        % cfg.operation     = 'x2-x1';
        % oscillatory = ft_math(cfg, fractal, original);
        % 
        % GroupFreqs{cnt,1} = oscillatory;
        % GroupFreqs{cnt,2} = original;
        cnt = cnt + 1;
        clear peak_times pta
    end
end

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))

[hl, ~] = boundedline(pta_time,squeeze(mean(mean(PTA_data{1}),2)),squeeze(std(mean(PTA_data{1}),[],2))./sqrt(size(PTA_data{1},2)),'k');
hl(1).LineWidth = 1.25;
ylabel('Amplitude (uV)')
xlabel('Time to Peak (s)')

sine_fit = fit(pta_time',squeeze(mean(mean(pta),2)),'sin1');


addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))

numSubjects = size(GroupFreqs, 1);

% Determine the layout grid size (e.g., 4 columns)
numCols = 4;
numRows = ceil(numSubjects / numCols);

% Create a tiled layout
figure;
t = tiledlayout(numRows, numCols, 'TileSpacing', 'Compact');

for iSubj = 1:numSubjects
    oscillatory = GroupFreqs{iSubj, 2};

    % Calculate the mean and standard error for the PSD
    mean_psd = squeeze(mean(mean(oscillatory.powspctrm, 2)));
    std_psd = squeeze(std(mean(oscillatory.powspctrm, 2), [], 1)) ./ sqrt(size(oscillatory.powspctrm, 1));

    % Plot the PSD for the current subject
    nexttile;
    [hl, ~] = boundedline(oscillatory.freq, mean_psd, std_psd, 'b');
    hl(1).LineWidth = 1.25;
    xlim([1 15]);
    ylabel('Power a.u.');
    xlabel('Frequency (Hz)');
    title(['Subject ' num2str(iSubj)]);
end

% Add a common title if desired
sgtitle('Single-Subject PSDs');