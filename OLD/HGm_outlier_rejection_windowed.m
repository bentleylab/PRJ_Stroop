% try with 10,10 windowing
% compare to voerall distribution, not specific to a given window
outlier_std_lim = 8;
win_len = 10;
win_step = 10;
win_lim = fn_sliding_window_lim(squeeze(hfa.powspctrm(1,1,1,:)),win_len,win_step);
bad_trials = cell([size(hfa.powspctrm,2) size(hfa.powspctrm,3)]);
bad_plot = 'segments';
ch_ix = 2;
% for ch_ix = 1:size(hfa.powspctrm,2)
    figure('Name',hfa.label{ch_ix},'units','normalized','outerposition',[0 0 1 1]);
    for f_ix = 1:size(hfa.powspctrm,3)
        bad_ix = cell([1 size(win_lim,1)]);
        avgs = NaN([size(hfa.powspctrm,1) size(win_lim,1)]);
        avgs_stds = NaN([1 size(win_lim,1)]);
        for win_ix = 1:size(win_lim,1)
            avgs(:,win_ix) = squeeze(mean(hfa.powspctrm(:,ch_ix,f_ix,win_lim(win_ix,1):win_lim(win_ix,2)),4));
            avgs_stds(win_ix) = std(avgs(:,win_ix));
            bad_ix{win_ix} = [bad_ix{win_ix} find(avgs(:,win_ix)>avgs_stds(win_ix)*outlier_std_lim)];
        end
        bad_trials{ch_ix,f_ix} = unique(vertcat(bad_ix{~cellfun(@isempty,bad_ix)})); %unique can't handle empty cells
        subplot(ceil(sqrt(size(hfa.powspctrm,3))),ceil(sqrt(size(hfa.powspctrm,3))),f_ix);
        hold on;
        % Plot all trials
        for t_ix = 1:size(hfa.powspctrm,1)
            plot(squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,:)),'k');
        end
        if strcmp(bad_plot,'segments')
            % Plot Individual window segments that were over the limit
            for win_ix = 1:size(win_lim,1)
                if ~isempty(bad_ix{win_ix})
                    for t_ix = 1:numel(bad_ix{win_ix})
                        plot(win_lim(win_ix,1):win_lim(win_ix,2),squeeze(hfa.powspctrm(bad_ix{win_ix}(t_ix),...
                            ch_ix,f_ix,win_lim(win_ix,1):win_lim(win_ix,2))),'r');%,'LineWidth',2);
                    end
                end
            end
        else
            for t_ix = 1:size(hfa.powspctrm,1)
                if any(bad_trials{ch_ix,f_ix}==t_ix)
                    plot(squeeze(hfa.powspctrm(t_ix,ch_ix,f_ix,:)),'r');%,'LineWidth',2);
                end
            end
        end
        title([hfa.label{ch_ix} ': CF=' num2str(round(hfa.freq(f_ix)),0)...
            '; #bad=' num2str(numel(bad_trials{ch_ix,f_ix}))]);
        set(gca,'XLim',[1 numel(hfa.powspctrm(1,1,1,:))]);
    end
    bad_trl_ch = vertcat(bad_trials{ch_ix,~cellfun(@isempty,bad_trials(ch_ix,:))});
    suptitle([num2str(numel(bad_trl_ch)) ' total bad; ' num2str(numel(unique(bad_trl_ch))) ' unique bad']);
% end