%% Fix EEG cleaning on raw data instead of preclean
SBJ = 'IR54';%'IR57';

% load the data
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
load(strcat(SBJ_vars.dirs.events,SBJ,'_eeg_bad_epochs_preclean.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_bob_bad_epochs_preclean.mat'));

cfg = [];
cfg.datafile = [SBJ_vars.dirs.raw SBJ_vars.raw_file];
cfg.channel  = 'DC01';
raw = ft_preprocessing(cfg);

load([SBJ_vars.dirs.preproc SBJ '_preclean.mat']);

%% Fix the epochs
% Find the matches to know which units are correct
bob_match = zeros([size(eeg_bad_epochs,1) 1]);
for ep_ix = 1:size(eeg_bad_epochs,1)
    if ~isempty(find(bad_epochs(:,1)==eeg_bad_epochs(ep_ix,1)))
        bob_match(ep_ix) = find(bad_epochs(:,1)==eeg_bad_epochs(ep_ix,1));
    end
end

% Correct only the bad ones
for ep_ix = 1:size(eeg_bad_epochs,1)
    if bob_match(ep_ix)==0
        eeg_bad_epochs(ep_ix,1) = eeg_bad_epochs(ep_ix,1) * data.fsample / raw.fsample;
        eeg_bad_epochs(ep_ix,2) = eeg_bad_epochs(ep_ix,2) * data.fsample / raw.fsample;
    end
end

% Convert to int for indexing
eeg_bad_epochs = round(eeg_bad_epochs);

% Re-sort to be in order
eeg_bad_epochs = sort(eeg_bad_epochs,1);

%% save fixed one
save(strcat(SBJ_vars.dirs.events,SBJ,'_eeg_bad_epochs_preclean.mat'),'eeg_bad_epochs');
