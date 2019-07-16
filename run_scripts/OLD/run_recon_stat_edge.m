if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Run fn_view_recon_stat_edge
pipeline_id = 'main_ft';
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
view_space = 'pat';
reg_type = '';
show_labels = 1;

%% R-locked
an_id = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';

fn_view_recon_stat_actvEdge('IR21',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   RT: 10 LPFC, 1 MPFC, 1 OFC
%   CNI: 2 INS, 5 LPFC (some overlap with huge row of 6 RT on L)
%   pcon: 2 LPFC/INS overlap with CNI, 1/2 with RT

fn_view_recon_stat_actvEdge('IR27',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   no RT, 1 OFC CNI, 1 LPFC pcon

fn_view_recon_stat_actvEdge('IR31',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   RT: 2 M1; 1 OFC CNI+pcon

% fn_view_recon_stat_actvEdge('IR32',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   FIX IH GRID!!!

fn_view_recon_stat_actvEdge('IR35',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   RT: 2-4 in M1, 6 DLPFC, 2 OFC, 1 MPFC
%   CNI: stron goverlpa with RT, slightly less coverage
%   pcon: 1 M1/PM and 2 DLPFC in common, +1 DLPFC

fn_view_recon_stat_actvEdge('IR39',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   RT: 1 INS?, 8 LPFC
%   CNI: 2LPFC (both overlap with RT)
%   pcon: 4LPFC (2 overlap with CNI)

fn_view_recon_stat_actvEdge('IR41',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   RT: 2 LPFC
%   CNI: 4 INS/LPFC + 3 LPFC (none overlap with RT)
%   pcon: 1 OFC

fn_view_recon_stat_actvEdge('IR48',pipeline_id,stat_id,an_id,'pat','',1,'r')
%   RT: whole row of DLPFC/FPC, some PM/M1/S1, and 3 SMA, 1 SPL
%   CNI: 1preSMA, 3 PM/M1, 1 SPL (not the same)
%   pcon: 3 SMA (1 RT), 3 IPL

fn_view_recon_stat_actvEdge('IR52',pipeline_id,stat_id,an_id,'pat','',1,'b')    %2 INS, 2 MPFC probes on R (+1 LAC)
%   RT: only 1...
%   CNI: 3 insula! 3 LPFC
%   pcon: 2 MPFC, 2 INS?, ~3 DLPFC

fn_view_recon_stat_actvEdge('IR54',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   5 MPFC RT, one is CNI and pcon (trifecta! that's the only pcon and CNI though)

fn_view_recon_stat_actvEdge('IR57',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   all LPFC RT and CNI, only 1 overlap, couple LPFC and OFC pcon

fn_view_recon_stat_actvEdge('IR61',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   couple DLPFC RT and CNI, but none overlap (no pcon)

% fn_view_recon_stat_actvEdge('IR65',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   NO RECON!

fn_view_recon_stat_actvEdge('IR68',pipeline_id,stat_id,an_id,'pat','',1,'l')
%   2 aMCC electrodes with RT+CNI+actv
%   CNI in LPFC and M1/PM

fn_view_recon_stat_actvEdge('IR72',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   RT: 3 DLPFC, 4 MPFC
%   CNI: 2 OFC, 11 LPFC, 8 MPFC
%   pcon: 3 OFC, 5 MPFC, 4 LPFC?

fn_view_recon_stat_actvEdge('IR74',pipeline_id,stat_id,an_id,'pat','',1,'b')
%   RT: 3 M1/S1/PM, 1 DLPFC
%   CNI: 4 M1/S1/PM, 3 DLPFC (all overlap with RT!)
%   pcon: 1 M1, 1 DLPFC (both overlap with CNI+RT), 1 INS?, 1 OFC, 1 PAR


%% Print good labels
SBJ = 'IR74';

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

load([SBJ_vars.dirs.import SBJ '_raw_labels.mat']);

%%
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

%%
good_lab = ft_channelselection({'all',bad_ch_neg{:},eeg_ch_neg{:},eog_ch_neg{:},photod_ch_neg{:},mic_ch_neg{:}},raw_labels);