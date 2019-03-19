function SBJ08a2_HFA_outlier_detection_summary(SBJ,an_id,actv_win,save_fig,fig_vis)%ch_id,steps_out,plot_bad
%% Compare bad_trials across processing steps
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load data
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
if strcmp(SBJ,'IR21')
    ch_id = 'ROI';
    steps_out = {{'bsln',10},{'sm',10},{'cmb',10},{'cmb_ns',10},{'raw',10}};
elseif strcmp(SBJ,'IR54')
    ch_id = 'LOFLAC';
    steps_out = {{'bsln',10},{'sm',10},{'cmb',10}};
elseif strcmp(SBJ,'IR57')
    ch_id = 'LAM456';
    steps_out = {{'bsln',10},{'sm',10},{'cmb',10},{'cmb_ns',10}};
end

data_filename = [SBJ_vars.dirs.proc SBJ '_actv_' ch_id '_' an_id '_mn' actv_win '_' steps_out{1}{1} '.mat'];
fprintf('Loading dataset %s\n',data_filename);
load(data_filename);
n_ch   = numel(hfa.label);
n_freq = size(hfa.powspctrm,3);

%% Aggregate results
bad_trl = cell([numel(steps_out) n_ch]);
bad_cnt = zeros([numel(steps_out) n_ch]);
lgd     = cell([numel(steps_out) n_ch]);
for step_ix = 1:numel(steps_out)
    bad_filename = [SBJ_vars.dirs.proc SBJ '_' an_id...
        '_outliers' num2str(steps_out{step_ix}{2}) '_' steps_out{step_ix}{1} '.mat'];
    tmp = load(bad_filename,'bad_trials');
    for ch_ix = 1:size(hfa.powspctrm,2)
        bad_trl{step_ix,ch_ix} = unique(horzcat(tmp.bad_trials{ch_ix,~cellfun(@isempty,tmp.bad_trials(ch_ix,:))}));
        bad_cnt(step_ix,ch_ix) = numel(bad_trl{step_ix,ch_ix});
    end
    lgd{step_ix} = [steps_out{step_ix}{1} num2str(steps_out{step_ix}{2})];
end

%% Plot Results
fig_dir = ['~/PRJ_Stroop/results/HFA/' SBJ '/actv/outlier_detection/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

figure;
subplot(2,1,1);
imagesc(bad_cnt);

%Plot params
colorbar;
caxis([1 20]);
ax = gca;
ax.YTick      = 1:size(bad_cnt,1);
ax.YTickLabel = lgd;
ax.XTick      = 1:size(bad_cnt,2);
ax.XTickLabel = hfa.label;
ax.Title.String = '# Trials Lost';

subplot(2,1,2);
overlap = zeros([numel(steps_out) n_ch]);
overlap(1,:) = 1;
for ch_ix = 1:n_ch
    for step_ix = 2:numel(steps_out)
        if ~isempty(bad_trl(step_ix,ch_ix)) && ~isempty(bad_trl(step_ix-1,ch_ix))
            tmp = ismember([bad_trl{step_ix,ch_ix}],[bad_trl{step_ix-1,ch_ix}]);
            overlap(step_ix,ch_ix) = sum(tmp)/numel(tmp);
        end
    end
end

imagesc(overlap);

%Plot params
colorbar;
caxis([0 1]);
ax = gca;
ax.YTick      = 1:size(overlap,1);
ax.YTickLabel = lgd;
ax.XTick      = 1:size(overlap,2);
ax.XTickLabel = hfa.label;
ax.Title.String = 'Overlap % Trials Lost';

