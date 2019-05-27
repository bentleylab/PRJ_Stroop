proc_id = 'main_ft';%'1000hz'i
atlas_id    = 'Dx';
roi_id      = 'gROI';

% PSD parameters
an.wn_len = 4096; %2048
an.foi    = [1 30];

% FOOOF params
an.pk_bw_lim = [1 12];%[0.5 12] is default
an.mn_pk_amp = 0.3;
an.max_n_pks = 3;
