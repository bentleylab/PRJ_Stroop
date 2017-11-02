function SBJ08a2_HFA_outlier_detection(SBJ,ch_id,an_id,actv_win,steps,outlier_cut
outlier_cut = 12;
% win_len = 10;
% win_step = 10;
% win_lim = fn_sliding_window_lim(squeeze(hfa.powspctrm(1,1,1,:)),win_len,win_step);
bad_trials = cell([size(hfa.powspctrm,2) size(hfa.powspctrm,3)]);
bad_plot = 'full';
% ch_ix = 2;
n_points = size(hfa.powspctrm,1)*size(hfa.powspctrm,4);
for ch_ix = 1:size(hfa.powspctrm,2)
    figure('Name',hfa.label{ch_ix},'units','normalized','outerposition',[0 0 1 1]);
    for f_ix = 1:size(hfa.powspctrm,3)
        avg = mean(reshape(hfa.powspctrm(:,ch_ix,f_ix,:),[1 n_points]));
        sd  = std(reshape(hfa.powspctrm(:,ch_ix,f_ix,:),[1 n_points]));
        cut = avg+sd*outlier_cut;
        unique(vertcat(bad_ix{~cellfun(@isempty,bad_ix)})); %unique can't handle empty cells
        subplot(ceil(sqrt(size(hfa.powspctrm,3))),ceil(sqrt(size(hfa.powspctrm,3))),f_ix);
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
                plot(squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,:)),'k');
            end
        end
        title([hfa.label{ch_ix} ': CF=' num2str(round(hfa.freq(f_ix)),0)...
            '; #bad=' num2str(numel(bad_trials{ch_ix,f_ix}))]);
        set(gca,'XLim',[1 numel(hfa.powspctrm(1,1,1,:))]);
    end
    bad_trl_ch = horzcat(bad_trials{ch_ix,~cellfun(@isempty,bad_trials(ch_ix,:))});
    suptitle([num2str(numel(bad_trl_ch)) ' total bad; ' num2str(numel(unique(bad_trl_ch))) ' unique bad']);
    
    fig_dir = ['~/PRJ_Stroop/results/HFA/' SBJ '/actv/outlier_detection/' an_id '/bsln_sm_cmb/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    fig_filename = [fig_dir SBJ '_actv_outliers_std' num2str(outlier_cut) '_' hfa.label{ch_ix} '.png'];
    saveas(gcf,fig_filename);
    
    bad_filename = [SBJ_vars.dirs.proc SBJ '_' an_id '_outliers' num2str(outlier_cut) '_bsln_sm_cmb.mat'];
    save(bad_filename,'bad_trials');
end

%% Compare bad_trials across processing steps
outlier_std_cuts = [12 12 12];%[10 10 10 8];
steps = {'bsln','bsln_sm','bsln_sm_cmb','bsln_sm_cmb'};
bad = cell([numel(steps) size(hfa.powspctrm,2)]);
bad_cnt = zeros([numel(steps) size(hfa.powspctrm,2)]);
for step_ix = 1:numel(steps)
    bad_filename = [SBJ_vars.dirs.proc SBJ '_' an_id...
        '_outliers' num2str(outlier_std_cuts(step_ix)) '_' steps{step_ix} '.mat'];
    tmp = load(bad_filename,'bad_trials');
    for ch_ix = 1:size(hfa.powspctrm,2)
        bad{step_ix,ch_ix} = unique(horzcat(tmp.bad_trials{ch_ix,~cellfun(@isempty,tmp.bad_trials(ch_ix,:))}));
        bad_cnt(step_ix,ch_ix) = numel(bad{step_ix,ch_ix});
    end
end

figure;
imagesc(bad_cnt);
colorbar;
caxis([1 20]);

overlap_b_bs = zeros([1 16]);
for ch_ix = 1:16
    if ~isempty(bad(2,ch_ix)) && ~isempty(bad(1,ch_ix))
        tmp = ismember([bad{2,ch_ix}],[bad{1,ch_ix}]);
        overlap_b_bs(ch_ix) = sum(tmp)/numel(tmp);
    end
end
