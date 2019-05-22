pipeline_id = 'srate';
atlas_id    = 'Dx';
roi_id      = 'gROI';

% PSD parameters
wn_len = 4096; %2048
foi    = [1 30];

% FOOOF params
pk_bw_lim = [1 12];%[0.5 12] is default
mn_pk_amp = 0.3;
max_n_pks = 3;
