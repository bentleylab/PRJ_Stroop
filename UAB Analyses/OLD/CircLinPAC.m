restoredefaultpath;
root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];

addpath(ft_dir);
ft_defaults;

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);

proc_id = 'main_ft';
SBJ = 'CP24';
%%
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
% eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);

load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
   
%% Select for Medial and Lateral PFC channels
atlas_id = 'Dx';
load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);

MPFC = elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC'));
LPFC = elec.label(ismember(elec.ROI,'dlPFC'));

cfgs.channel = [MPFC;LPFC];
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
trial_lim_s_pad = [-4 4];

stim_events = trial_info.word_onset;

roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
    round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));

% Get theta phases
cfg              = [];
cfg.bpfilter     = 'yes';
cfg.bpfreq       = [4 8];
cfg.hilbert      = 'angle';
theta_phase      = ft_preprocessing(cfg, roi_trl);

% Trim buffer periods off
cfg_trim = [];
cfg_trim.latency = [-2.5 2.5];
theta_phase = ft_selectdata(cfg_trim,theta_phase);

phases = cat(3,theta_phase.trial{:});
phases = squeeze(permute(phases, [3 2 1]));

trials = fn_create_LMM_design(trial_info,1,1,'power');
C = double(strcmpi(trials.CurrentType,'inc')); % conflict trials ('inc') = 1 and nonconflict trials ('con' and 'neu') = 0

phase_time = theta_phase.time{1}(theta_phase.time{1} >= -0.5 & theta_phase.time{1} <= 1.25);

baseline_window = [-0.5 -0.3]; 

base_inds = phase_time >= baseline_window(1) & phase_time < baseline_window(2);
%% Theta Phase Encoding
addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/Circular Statistics Toolbox (Directional Statistics)'))

cnt = 0;
for mpfci = 1:length(MPFC)
    for lpfci = 1:length(LPFC)
        
        t = 1;
        for tidx = nearest(theta_phase.time{1}, phase_time(1)) : ...
                   nearest(theta_phase.time{1}, phase_time(end))
            
            % Phase difference for this pair
            phase_diff = phases(:, :, strcmpi(theta_phase.label, MPFC{mpfci})) ...
                       - phases(:, :, strcmpi(theta_phase.label, LPFC{lpfci}));
            
            % circ_corrcl at this time point
            [rho(t,1), pval(t,1)] = circ_corrcl( phase_diff(:, tidx), C );
            t = t + 1;
        end
        
        % Baseline correction
        baseline_rho = mean(rho(base_inds));
        adjusted_rho = rho - baseline_rho;  

        clusters = bwconncomp(pval < 0.05);
        
        for clusteri = 1:clusters.NumObjects

            cluster_extent = clusters.PixelIdxList{clusteri};

            if length(cluster_extent) < round(0.1 * roi_trl.fsample)
                continue
            elseif all( phase_time(cluster_extent) >= 0 )
                
                cnt = cnt + 1;
                sig_pairs{cnt,1} = MPFC{mpfci};
                sig_pairs{cnt,2} = LPFC{lpfci};
                
                % Time of significant phase encoding
                sig_time(cnt,1) = phase_time(cluster_extent(1));
                
                % Store baseline-adjusted correlation
                phase_encoding(cnt,:) = adjusted_rho; 
                
                break  % only consider first cluster for this pair
            end
        end
        
        clear rho pval adjusted_rho baseline_rho
    end
end
%% PAC
% Get HFA
cfg              = [];
cfg.bpfilter     = 'yes';
cfg.bpfreq       = [70 150];
cfg.hilbert      = 'abs';
hfa      = ft_preprocessing(cfg, roi_trl);

% Trim buffer periods off
cfg_trim = [];
cfg_trim.latency = [-2.5 2.5];
hfa = ft_selectdata(cfg_trim,hfa);

hfas = cat(3,hfa.trial{:});
hfas = squeeze(permute(hfas, [3 2 1]));

PAC_s_time = theta_phase.time{1}(theta_phase.time{1} >= -0.5 & ...
                                  theta_phase.time{1} <= 1.25);

PAC_stim = cell(size(sig_pairs,1),1);

for sp = 1:size(sig_pairs,1)

    t = 1;
    for tidx = nearest(theta_phase.time{1}, PAC_s_time(1)):...
               nearest(theta_phase.time{1}, PAC_s_time(end))
       
        %% Forward direction: MPFC phase -> LPFC amplitude
        MPFCPhase = phases(:, tidx, strcmpi(theta_phase.label, sig_pairs{sp,1}));  
        LPFCHFA   = hfas(:,   tidx, strcmpi(hfa.label,         sig_pairs{sp,2}));  
        
        cTrials = (C == 1);
        [CRho_fwd, ~] = circ_corrcl(MPFCPhase(cTrials), LPFCHFA(cTrials));
        
        ncTrials = (C == 0);
        [NCRho_fwd, ~] = circ_corrcl(MPFCPhase(ncTrials), LPFCHFA(ncTrials));
        
        rho_conflict_fwd(t)    = CRho_fwd;
        rho_nonconflict_fwd(t) = NCRho_fwd;
        
        %% Backward direction: LPFC phase -> MPFC amplitude
        
        LPFCPhase = phases(:, tidx, strcmpi(theta_phase.label, sig_pairs{sp,2}));
        MPFCHFA   = hfas(:,   tidx, strcmpi(hfa.label,        sig_pairs{sp,1}));
        
        [CRho_bwd, ~] = circ_corrcl(LPFCPhase(cTrials), MPFCHFA(cTrials));
        [NCRho_bwd, ~] = circ_corrcl(LPFCPhase(ncTrials), MPFCHFA(ncTrials));
        
        rho_conflict_bwd(t)    = CRho_bwd;
        rho_nonconflict_bwd(t) = NCRho_bwd;
        
        t = t + 1;
    end

    %% Store both directions
    PAC_stim{sp}.pair = {sig_pairs{sp,1}, sig_pairs{sp,2}};
    PAC_stim{sp}.time = PAC_s_time;

    % Forward
    PAC_stim{sp}.rho_conflict_fwd    = rho_conflict_fwd;
    PAC_stim{sp}.rho_nonconflict_fwd = rho_nonconflict_fwd;

    % Backward
    PAC_stim{sp}.rho_conflict_bwd    = rho_conflict_bwd;
    PAC_stim{sp}.rho_nonconflict_bwd = rho_nonconflict_bwd;
    
    clear rho_conflict_fwd rho_nonconflict_fwd ...
          rho_conflict_bwd rho_nonconflict_bwd
end

encoding_timevec = -0.2 : 1/roi_trl.fsample : 0.2;

PAC_encoding = cell(size(sig_pairs,1),1);

for sp = 1:size(sig_pairs,1)
    
    % Time window around this channel pair's encoding onset
    encoding_window = [sig_time(sp) - 0.2, sig_time(sp) + 0.2];
    encoding_idxs = nearest(theta_phase.time{1}, encoding_window(1)) : ...
                  nearest(theta_phase.time{1}, encoding_window(2));
    
    t = 1;
    for tidx = encoding_idxs
        
        %% ------------------------
        %  FORWARD direction: 
        %  MPFC phase -> LPFC amplitude
        %% ------------------------
        MPFCPhase = phases(:, tidx, strcmpi(theta_phase.label, sig_pairs{sp,1}));
        LPFCHFA   = hfas(:,   tidx, strcmpi(hfa.label,         sig_pairs{sp,2}));
        
        cTrials = (C == 1);
        [CRho_fwd, ~] = circ_corrcl(MPFCPhase(cTrials), LPFCHFA(cTrials));
        
        ncTrials = (C == 0);
        [NCRho_fwd, ~] = circ_corrcl(MPFCPhase(ncTrials), LPFCHFA(ncTrials));
        
        rho_conflict_fwd(t)    = CRho_fwd;
        rho_nonconflict_fwd(t) = NCRho_fwd;
        
        %% ------------------------
        %  BACKWARD direction:
        %  LPFC phase -> MPFC amplitude
        %% ------------------------
        LPFCPhase = phases(:, tidx, strcmpi(theta_phase.label, sig_pairs{sp,2}));
        MPFCHFA   = hfas(:,   tidx, strcmpi(hfa.label,         sig_pairs{sp,1}));
        
        [CRho_bwd, ~] = circ_corrcl(LPFCPhase(cTrials), MPFCHFA(cTrials));
        [NCRho_bwd, ~] = circ_corrcl(LPFCPhase(ncTrials), MPFCHFA(ncTrials));
        
        rho_conflict_bwd(t)    = CRho_bwd;
        rho_nonconflict_bwd(t) = NCRho_bwd;
        
        t = t + 1;
    end
    
    %% Store results for this pair
    PAC_encoding{sp}.pair = {sig_pairs{sp,1}, sig_pairs{sp,2}};
    PAC_encoding{sp}.time = encoding_timevec;
    
    % Forward direction
    PAC_encoding{sp}.rho_conflict_fwd    = rho_conflict_fwd;
    PAC_encoding{sp}.rho_nonconflict_fwd = rho_nonconflict_fwd;
    
    % Backward direction
    PAC_encoding{sp}.rho_conflict_bwd    = rho_conflict_bwd;
    PAC_encoding{sp}.rho_nonconflict_bwd = rho_nonconflict_bwd;
    
    clear rho_conflict_fwd rho_nonconflict_fwd ...
          rho_conflict_bwd rho_nonconflict_bwd
end

