function [filt_r] = fn_realign_filt_s2r(filt_s,rt,trial_lim_sec)
%% Realign a stimulus-locked TFR to the responses indicated by offsets
% INPUTS:
%   filt_s [ft struct] - output of ft_preprocessing to be realigned
%   rt [int array] - n_trialsx1 array of RTs (TIMES in SEC) within that trial
%   trial_lim [int, int] - [start, end] array of TIME (in SEC) indices of epoch around RT
% OUTPUTS:
%   filt_r [ft_struct] - realigned tfr with correct .trial and .time fields

% Check failure cases
if ~isa(filt_s.time{1},'double')
    error('tfr_s.time field is not a double. Trials are probably not the same length!');
end
for t = 1:numel(filt_s.time)
    if sum(filt_s.time{t}==filt_s.time{1})~=numel(filt_s.time{1})
        error('Not all trials are same length!');
    end
end

% Compute output trial/time size
[~, t1_beg_ix] = min(abs(filt_s.time{1}-(rt(1)+trial_lim_sec(1))));
[~, t1_end_ix] = min(abs(filt_s.time{1}-(rt(1)+trial_lim_sec(2))));
[~, t1_rt_ix]  = min(abs(filt_s.time{1}-rt(1)));
new_time = filt_s.time{1}(t1_beg_ix:t1_end_ix)-filt_s.time{1}(t1_rt_ix);

% Create new time field
filt_r = filt_s;
% orig_size = size(filt_s.powspctrm);
% filt_r.powspctrm = NaN([orig_size(1) orig_size(2) orig_size(3) numel(filt_r.time)]);
for t = 1:numel(filt_s.trial)
    filt_r.time{t} = new_time;
    [~, beg_ix] = min(abs(filt_s.time{t}-(rt(t)+trial_lim_sec(1))));
    [~, end_ix] = min(abs(filt_s.time{t}-(rt(t)+trial_lim_sec(2))));
    filt_r.trial{t} = filt_s.trial{t}(:,beg_ix:end_ix);
end

end