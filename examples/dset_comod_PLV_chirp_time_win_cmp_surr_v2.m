function dset_comod_pow_chirp_time_win_cmp_surr_v2(DATASET, event1,...
    win_start1, win_end1, event2, win_start2, win_end2, run_surrogate, varargin)
% win_cmp v2:
%   -changed output directories and filenames to have the directory list
%   the two windows and files to not have that info
%   -NOTE v1 was changed by adding som comments and a little shifting of
%   stuff, most notably by changing the name of the surrogate function and
%   adding documentation to that script
% win_cmp v1 now takes two windows:
%   1) calculates PAC in each using the epoched time method (trial cuts,
%       average across trials for PAC time series)
%   2) averages PAC time series to get final PAC value
%   3) takes the PAC time series from BOTH windows and 
%
% modified from MAT03_dset_comod_PLV_chirp_canolty_circ_surr_v4
% v4 changes:
%   -phase_f_array defaults to 2:30
%   -saves the phase_f_array and amp_f_array with the PAC results for easy plotting
% v3 is SGE compatible, does not do any plotting
%
% As of v2 11/28/15, CWH added his own PAC by PLV functions using the same
% chirp filtered time series, both averaging across time and trials.
% Includes surrogates via circ shift for Canolty PAC and random epochs for CWH PLV.
% Thus, the following are calculated and saved out:
%   PLV across all windows at once + circ stats (canolty)
%   PLV across trials + epoch stats (CWH)
%   PLV across time + epoch stats (CWH)
%   
% Originally modeled after plot_CFC.m, which includes calculations and 
%   plotting for Fig. 1D from Canolty et al. (2006) for two electrodes.
%
%  INPUTS:
%     DATASET           - 'SBJ_task.##' string (ex: 'GP15_EmoGen.01'), split into SBJ, task, and e_num
%     event             - string referring to trial event to lock iwndow to
%                           ('Resp', 'Stim'; future could be 't2')
%     win_start         - string determining when the analysis window starts relative to the event time
%                           (e.g. '-300', which is in ms and gets converted via sampling rate)
%     win_end           - string determining when the analysis window ends relative to the event time
%                           (e.g. '0' would mean end at the event time, in ms and converted via sampling rate)
%     *2                - same as event, win_start, win_end above, but for the second window
%     run_surrogate     - 0 = do not run surrogate analysis. 1 = run surrogate analysis
%
%  OPTIONAL:
%     amp_f_array       - 40 frequencies to use for amplitude
%     phase_f_array     - 19 frequencies to use for phase
%     nboots            - how many iterations for surrogate distribution (default 500)
%
%  EXAMPLE:
% [plv_trials] = plot_CFC('ST22','decision','e1.mat',1,'Events',Events);
%
%     Based upon scripts by Ryan Canolty.  Modified by Dan Bliss and Sara
%     Szczepanski, 2/23/12, Matar Haller 8/14/12
%     Completely reworked by Colin Hoy 11/28/15

%% Set up folders 
PAC_tag  = 'PLV_chirp';
win_tag1  = strcat(event1,'_',win_start1,'.',win_end1);
win_tag2  = strcat(event2,'_',win_start2,'.',win_end2);
run_surrogate = str2num(run_surrogate);

addpath('/home/knight/duration_PLV/Scripts/');
addpath('/home/knight/duration_PLV/Scripts/canolty_CFC/');

PrcsDataDir = '/home/knight/duration_PLV/PrcsData/';
uScore      = strfind(DATASET,'_'); period = strfind(DATASET,'.');
SBJ         = DATASET(1:uScore-1);
task        = DATASET(uScore+1:period-1);
e_num       = str2num(DATASET(period+1:end));
SBJ_dir = strcat(PrcsDataDir,SBJ,'/',task,'/');
% PAC_dir_win1 = strcat(SBJ_dir,'PAC_ePairs_trials/',PAC_tag,'/',win_tag1,'/');
% PAC_dir_win2 = strcat(SBJ_dir,'PAC_ePairs_trials/',PAC_tag,'/',win_tag2,'/');
PAC_dir_diff = strcat(SBJ_dir,'PAC_ePairs_trials/',PAC_tag,'/',win_tag1,'-',win_tag2,'/');

% For the sake of saving space and not saving the same thing as
% dset_comod_PLV_chirp_canolty_circ_surr_v4.m, I'm not saving the PLVs,
% only the difference in PLV between windows
% if ~exist(PAC_dir_win1,'dir') % checks if the appropriate subfolder for this study/block has been created in 'analysis'....
%     mkdir(PAC_dir_win1);
% end
% if ~exist(PAC_dir_win2,'dir') % checks if the appropriate subfolder for this study/block has been created in 'analysis'....
%     mkdir(PAC_dir_win2);
% end
if ~exist(PAC_dir_diff,'dir') % checks if the appropriate subfolder for this study/block has been created in 'analysis'....
    mkdir(PAC_dir_diff);
end
% To make things simpler, I'm saving the pval matrix in the same .mat file
% as the PLV results
if run_surrogate %if we are creating a surrogate distribution
    surr_dir = strcat(SBJ_dir,'PAC_ePairs_trials_surr/',PAC_tag,'/win_cmp/');
    if ~exist(surr_dir,'dir')
        mkdir(surr_dir);
    end
end


%% Set defaults
for n=1:2:length(varargin)-1
    switch lower(varargin{n})
        case 'amp_f_array'
%             amp_f_array = varargin{n+1};
            error('ERROR: Passed amp_f_array in SGE version of script...');
        case 'phase_f_array'
%             phase_f_array = varargin{n+1};
            error('ERROR: Passed phase_f_array in SGE version of script...');
        case 'nboots'
            nboots = str2num(varargin{n+1});
    end
end

if ~exist('amp_f_array', 'var')
    amp_f_array = 5:5:200;
end
if ~exist('phase_f_array', 'var')
    phase_f_array = 2:30;
end
if ~exist('nboots','var')
    nboots = 500;
end

%% Load Data
load(strcat(SBJ_dir,'subj_globals.mat'),'srate');
load(strcat(SBJ_dir,'gdat_notch.mat'));
elec = gdat(e_num,:);

win_start_dp1 = str2num(win_start1)*srate/1000;
win_end_dp1   = str2num(win_end1)*srate/1000;
[starts1,ends1] = func_return_event_wins(SBJ_dir,task,event1,win_start_dp1,win_end_dp1);

win_start_dp2 = str2num(win_start2)*srate/1000;
win_end_dp2   = str2num(win_end2)*srate/1000;
[starts2,ends2] = func_return_event_wins(SBJ_dir,task,event2,win_start_dp2,win_end_dp2);


%% get_signal_parameters, which returns the structure 'sp':
sp = get_signal_parameters('sampling_rate',srate,... 
    'number_points_time_domain',length(elec));

n_amp_freqs = length(amp_f_array); %number of freq for amp
n_phase_freqs = length(phase_f_array); %numer of freq for phase

% plv_time1            = nan(n_amp_freqs, n_phase_freqs, length(elec(starts1(1):ends1(1))));
% plv_time_mean1       = nan(n_amp_freqs, n_phase_freqs);
% plv_time2            = nan(n_amp_freqs, n_phase_freqs, length(elec(starts2(1):ends2(1))));
% plv_time_mean2       = nan(n_amp_freqs, n_phase_freqs);

plv_time_mean_diff   = nan(n_amp_freqs, n_phase_freqs); %difference in mean PLV across windows
plv_time_mean_diff_surr = nan(n_amp_freqs, n_phase_freqs, nboots);
plv_time_mean_diff_pval = nan(n_amp_freqs, n_phase_freqs);


%% loop over phase and amplitude frequencies
e_name =  strcat('e',num2str(e_num,'%03d'));%elec(1:end-4); %string (ex: 'e14')

for i_phase = 1:n_phase_freqs %for each frequency of the phase data
    phase_f = phase_f_array(i_phase); %freq to filter lf signal
    fprintf('\n phase data freq : %i\n', phase_f)
    
    g.center_frequency = phase_f;
    g.fractional_bandwidth = 0.25;
    g.chirp_rate = 0;
    g1 = make_chirplet('chirplet_structure', g, 'signal_parameters', sp);
    % filter raw signal at low frequency, extract phase:
    fs = filter_with_chirplet('raw_signal', elec, ...
        'signal_parameters', sp, ...
        'chirplet', g1);
    lf_phase = angle(fs.time_domain); %phase
    clear g.center_frequency
    
    fprintf('amp data freq : \n')
    for i_amp = 1:n_amp_freqs %for each frequency of the amplitude data
        amp_f = amp_f_array(i_amp); %freq to filter hf data
        fprintf('%i... ',amp_f);
        
        g.center_frequency = amp_f;
        g2 = make_chirplet('chirplet_structure', g, 'signal_parameters', sp);
        % filter raw signal at high frequency, extract amplitude:
        fs = filter_with_chirplet('raw_signal', elec, ...
            'signal_parameters', sp, ...
            'chirplet', g2);
        hf_amp = abs(fs.time_domain); %amp
        % filter high frequency amplitude time-series at low frequency, extract phase:
        fs = filter_with_chirplet('raw_signal', hf_amp, ...
            'signal_parameters', sp, ...
            'chirplet', g1);%filter at low frequency
        hf_phase = angle(fs.time_domain); %extract phase of high frequency amplitude
                
        % Compute cross-frequency phase locking value (PLV).
        plv_ts1 = func_PAC_PLVenvPh_ePair_time(lf_phase, hf_phase, starts1, ends1,'return_complex', 1);
        plv_ts2 = func_PAC_PLVenvPh_ePair_time(lf_phase, hf_phase, starts2, ends2,'return_complex', 1);
        plv1 = mean(abs(plv_ts1));
        plv2 = mean(abs(plv_ts2));
%         plv_time1(i_amp, i_phase, :) = func_PAC_PLVenvPh_ePair_time(lf_phase, hf_phase, starts1, ends1,...
%             'return_complex', 1);
%         plv_time2(i_amp, i_phase, :) = func_PAC_PLVenvPh_ePair_time(lf_phase, hf_phase, starts2, ends2,...
%             'return_complex', 1);
        
        plv_time_mean_diff(i_amp, i_phase) = plv1 - plv2;

        if run_surrogate
            [plv_time_mean_diff_pval(i_amp,i_phase),...
                plv_time_mean_diff_surr(i_amp,i_phase,:)] = func_PLV_time_surr_epoch_diff_rand_PLV_ts(...
                plv_time_mean_diff(i_amp,i_phase),plv_ts1,plv_ts2,nboots);
        end %end if run_surrogate        
    end %end hf amp for loop    
end %end lf phase for loop

% PAC_time_filename1 = strcat(PAC_dir_win1,'PAC_',PAC_tag,'_time_',e_name,'.mat');
% fprintf('saving : %s\n', PAC_time_filename1)
% save(PAC_time_filename1, 'plv_time1', 'phase_f_array', 'amp_f_array');

% PAC_time_filename2 = strcat(PAC_dir_win2,'PAC_',PAC_tag,'_time_',e_name,'.mat');
% fprintf('saving : %s\n', PAC_time_filename2)
% save(PAC_time_filename2, 'plv_time2', 'phase_f_array', 'amp_f_array');

if run_surrogate %in case we need to use these for another analysis later
    PAC_time_diff_filename = strcat(PAC_dir_diff,'PAC_',PAC_tag,'_time_win_cmp_',e_name,'_wPval.mat');
    fprintf('saving : %s\n', PAC_time_diff_filename)
    save(PAC_time_diff_filename, 'plv_time_mean_diff', 'plv_time_mean_diff_pval', 'phase_f_array', 'amp_f_array');

    surr_diff_filename = strcat(surr_dir,'PAC_',PAC_tag,'_time_win_cmp_surr_dist_',e_name,...
        '_nboots',num2str(nboots),'.mat');
    fprintf('saving : %s\n', surr_diff_filename)
    save(surr_diff_filename, 'plv_time_mean_diff_surr');
else
    PAC_time_diff_filename = strcat(PAC_dir_diff,'PAC_',PAC_tag,'_time_win_cmp_',e_name,'_nPval.mat');
    fprintf('saving : %s\n', PAC_time_diff_filename)
    save(PAC_time_diff_filename, 'plv_time_mean_diff', 'phase_f_array', 'amp_f_array');
end

end


