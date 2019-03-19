function SBJ08a2_HFA_outlier_detection(SBJ,an_id,actv_win,save_fig,fig_vis)%ch_id,steps_out,plot_bad
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

bad_plot = 'full';
if strcmp(SBJ,'IR21')
    ch_id = 'ROI';
    steps_out = {{'bsln',10},{'sm',10},{'cmb',10},{'cmb_ns',10},{'raw',10}};
elseif strcmp(SBJ,'IR31')
    ch_id = 'LINRACLOF';
    steps_out = {{'bsln',10},{'sm',10},{'cmb',10},{'cmb_ns',10},{'raw',10}};
elseif strcmp(SBJ,'IR54')
    ch_id = 'LOFLAC';
    steps_out = {{'bsln',10},{'sm',10},{'cmb',10}};
elseif strcmp(SBJ,'IR57')
    ch_id = 'LAM456';
    steps_out = {{'bsln',10},{'sm',10},{'cmb',10},{'cmb_ns',10}};
end

%%
if isnumeric(actv_win); actv_win = num2str(actv_win); end
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

for step_ix = 1:numel(steps_out)
    % Load Data
    step = steps_out{step_ix}{1};
    outlier_cut = steps_out{step_ix}{2};
    data_filename = [SBJ_vars.dirs.proc SBJ '_actv_' ch_id '_' an_id '_mn' actv_win '_' step '.mat'];
    fprintf('Loading dataset %s\n',data_filename);
    load(data_filename);
    n_ch   = numel(hfa.label);
    ch_lab = hfa.label;
    n_freq = size(hfa.powspctrm,3);
    fig_dir = ['~/PRJ_Stroop/results/HFA/' SBJ '/actv/outlier_detection/' an_id '/' step '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    bad_trials = cell([n_ch n_freq]);
    n_points = size(hfa.powspctrm,1)*size(hfa.powspctrm,4);
    for ch_ix = 1:n_ch
        figure('Name',hfa.label{ch_ix},'units','normalized',...
            'outerposition',[0 0 1 1],'Visible',fig_vis);
        for f_ix = 1:n_freq
            avg = mean(reshape(hfa.powspctrm(:,ch_ix,f_ix,:),[1 n_points]));
            sd  = std(reshape(hfa.powspctrm(:,ch_ix,f_ix,:),[1 n_points]));
            cut = avg+sd*outlier_cut;
%             unique(vertcat(bad_ix{~cellfun(@isempty,bad_ix)})); %unique can't handle empty cells
            subplot(ceil(sqrt(n_freq)),ceil(sqrt(n_freq)),f_ix);
            hold on;
            % Plot all trials
            for t_ix = 1:size(hfa.powspctrm,1)
                if any(squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,:))>cut)
                    bad_trials{ch_ix,f_ix} = [bad_trials{ch_ix,f_ix} t_ix];
                    if strcmp(bad_plot,'segments')
                        % Plot Individual window segments that were over the limit
                        bad_time_idx = squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,:))>cut;
                        bad_chunks = fn_find_chunks(bad_time_idx);
                        bad_chunks(squeeze(bad_time_idx(bad_chunks(:,1)))==0,:) = [];
                        plot(squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,:)),'k');
                        for chunk_ix = 1:size(bad_chunks,1)
                            win_tmp = bad_chunks(chunk_ix,1):bad_chunks(chunk_ix,2);
                            plot(win_tmp,squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,win_tmp)),'r');
                        end
                    else
                        plot(squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,:)),'r');
                    end
                    line([1 numel(hfa.powspctrm(1,1,1,:))],[cut cut],'Color','r','LineStyle','--','LineWidth',2);
                else
                    if numel(bad_trials{ch_ix,f_ix})<1  % is this necessary? whatever...
                        bad_trials{ch_ix,f_ix} = [];
                    end
                    plot(squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,:)),'k');
                end
            end
            title([hfa.label{ch_ix} ': CF=' num2str(round(hfa.freq(f_ix)),0)...
                '; #bad=' num2str(numel(bad_trials{ch_ix,f_ix}))]);
            set(gca,'XLim',[1 numel(hfa.powspctrm(1,1,1,:))]);
        end
        bad_trl_ch = horzcat(bad_trials{ch_ix,~cellfun(@isempty,bad_trials(ch_ix,:))});
        suptitle([num2str(numel(bad_trl_ch)) ' total bad; ' num2str(numel(unique(bad_trl_ch))) ' unique bad']);
        
        if save_fig
            fig_filename = [fig_dir SBJ '_actv_outliers_std' num2str(outlier_cut) '_' hfa.label{ch_ix} '.png'];
            saveas(gcf,fig_filename);
        end
        
        bad_filename = [SBJ_vars.dirs.proc SBJ '_' an_id '_outliers' num2str(outlier_cut) '_' step '.mat'];
        save(bad_filename,'bad_trials');
    end
    clear hfa bad_trials bad_trl_ch
end

%% Aggregate results
bad_trl = cell([numel(steps_out) n_ch]);
bad_cnt = zeros([numel(steps_out) n_ch]);
lgd     = cell([numel(steps_out) n_ch]);
for step_ix = 1:numel(steps_out)
    bad_filename = [SBJ_vars.dirs.proc SBJ '_' an_id...
        '_outliers' num2str(steps_out{step_ix}{2}) '_' steps_out{step_ix}{1} '.mat'];
    tmp = load(bad_filename,'bad_trials');
    for ch_ix = 1:n_ch
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
ax.XTickLabel = ch_lab;
ax.XTickLabelRotation = 45;
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
ax.XTickLabel = ch_lab;
ax.XTickLabelRotation = 45;
ax.Title.String = 'Overlap % Trials Lost';

if save_fig
    fig_filename = [fig_dir SBJ '_actv_outliers_summary.png'];
    saveas(gcf,fig_filename);
end

end