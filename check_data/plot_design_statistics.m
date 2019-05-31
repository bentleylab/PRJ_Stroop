function plot_design_statistics(SBJ)
%% Compute stimulus and design statistics
%   Proportion congruence and stimulus-response associations
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);

%% Load Data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{block});
else
    block_suffix = SBJ_vars.block_name{block};   % should just be ''
end

% Load trial_info
load([SBJ_vars.dirs.events SBJ '_trial_info_final.mat']);
n_trl = numel(trial_info.trial_n);

%% Get Design Information
% Get stimulus and response options
words = unique(trial_info.word);
colors = unique(trial_info.color);
wc_lab_mat = cell([numel(words) numel(colors)]);
for w_ix = 1:numel(words)
    for c_ix = 1:numel(colors)
        wc_lab_mat{w_ix,c_ix} = [words{w_ix} '-' colors{c_ix}];
    end
end

word_idx = zeros([n_trl 1]);
for w_ix = 1:numel(words)
    word_idx(strcmp(trial_info.word,words{w_ix})) = w_ix;
end
color_idx = zeros([n_trl 1]);
for c_ix = 1:numel(colors)
    color_idx(strcmp(trial_info.color,colors{c_ix})) = c_ix;
end

% Get design options
[cni_lab, cni_colors, ~] = fn_condition_label_styles('CNI');
cni_idx = zeros([n_trl 1]);
for cni_ix = 1:numel(cni_lab)
    cni_idx(logical(fn_condition_index(cni_lab{cni_ix},trial_info.condition_n))) = cni_ix;
end
[pc_lab, pc_colors, ~] = fn_condition_label_styles('PC');
pc_idx = zeros([n_trl 1]);
for pc_ix = 1:numel(pc_lab)
    pc_idx(logical(fn_condition_index(pc_lab{pc_ix},trial_info.condition_n))) = pc_ix;
end

%% Compute Design Statistics
% Compute overall word-color frequency
wc_mat = zeros([numel(words) numel(colors)]);
for trl_ix = 1:n_trl
    wc_mat(word_idx(trl_ix),color_idx(trl_ix)) = wc_mat(word_idx(trl_ix),color_idx(trl_ix))+1;
end

% Compute stimulus-response congruence frequency
full_mat = zeros([numel(words) numel(colors) numel(cni_lab) numel(pc_lab)]);
for trl_ix = 1:n_trl
    full_mat(word_idx(trl_ix),color_idx(trl_ix),cni_idx(trl_ix),pc_idx(trl_ix)) = ...
        full_mat(word_idx(trl_ix),color_idx(trl_ix),cni_idx(trl_ix),pc_idx(trl_ix))+1;
end


%% Plot Stimulus-Response Frequency Matrix
% Word-Color frequency
fig_name = 'S-R_frequency';
figure('Name',fig_name);
imagesc(wc_mat./n_trl);
ax = gca; %ax.YDir = 'normal';
ax.XTick = 1:numel(colors);
ax.XTickLabel = colors;
xlabel('Colors');
ax.YTick = 1:numel(words);
ax.YTickLabel = words;
ylabel('Words');
colorbar;

% Word-Color CNI frequency
%   Unnecessary, because C and N columns are already plotted in S-R matrix above
% tmp_mat = sum(full_mat,4);
% fig_name = 'S-R_CNI_frequency';
% figure('Name',fig_name);
% imagesc(reshape(tmp_mat,numel(words)*numel(colors),numel(cni_lab))./n_trl);
% ax = gca; %ax.YDir = 'normal';
% ax.XTick = 1:numel(cni_lab);
% ax.XTickLabel = cni_lab;
% xlabel('Congruence');
% ax.YTick = 1:numel(words)*numel(colors);
% ax.YTickLabel = reshape(wc_lab_mat,numel(words)*numel(colors),1);
% ylabel('Stimulus-Response');
% colorbar;

% Word-Color Incongruent frequency
tmp_mat = squeeze(sum(full_mat(~strcmp(words,'XXXX'),:,strcmp(cni_lab,'I'),:),4));
fig_name = 'S-R_I_frequency';
figure('Name',fig_name);
imagesc(tmp_mat./n_trl);
ax = gca; %ax.YDir = 'normal';
ax.XTick = 1:numel(colors);
ax.XTickLabel = colors;
xlabel('Colors');
ax.YTick = 1:numel(words)-1;
ax.YTickLabel = words(~strcmp(words,'XXXX'));
ylabel('Words');
colorbar;

% CNI across PC
cond_mat = squeeze(sum(sum(full_mat,1),2));
for pc_ix = 1:numel(pc_lab)
    cond_mat(:,pc_ix) = cond_mat(:,pc_ix)./sum(pc_idx==pc_ix);
end
fig_name = 'CNI_PC_frequency';
figure('Name',fig_name);
imagesc(cond_mat);
ax = gca; %ax.YDir = 'normal';
ax.XTick = 1:numel(pc_lab);
ax.XTickLabel = pc_lab;
xlabel('PC');
ax.YTick = 1:numel(cni_lab);
ax.YTickLabel = cni_lab;
ylabel('CNI');
colorbar;

%% Bar graphs for visual inspection
for grp_ix = 1:2
    if grp_ix==1
        fig_name = 'CNI_word-color_design';
        figure('Name',fig_name,'outerposition',[0 0 1 1]);
        cond_lab = cni_lab;
        cond_colors = cni_colors;
        mat = squeeze(sum(full_mat,4));
    else 
        fig_name = 'PC_word-color_design';
        figure('Name',fig_name,'outerposition',[0 0 1 1]);
        cond_lab = pc_lab;
        cond_colors = pc_colors;
        mat = squeeze(sum(full_mat,3));
    end
    
    % Word bar plot
    subplot(3,1,1); hold on;
    tmp_mat = squeeze(sum(mat,2));
    b = cell(size(words));
    bar_offsets = linspace(-0.25,0.25,numel(cond_lab));
    for w_ix = 1:numel(words)
        b{w_ix} = bar(bar_offsets+w_ix,diag(tmp_mat(w_ix,:)),0.9,'stacked');
        for cond_ix = 1:numel(cond_lab)
            set(b{w_ix}(cond_ix),'FaceColor',cond_colors{cond_ix},'EdgeColor','k');
        end
    end
    set(gca,'XTick',1:numel(words));
    set(gca,'XTickLabels',words);
    legend(b{1},cond_lab,'Location','best');
    title('Word-CNI Frequency');
    
    % Color bar plot
    subplot(3,1,2); hold on;
    tmp_mat = squeeze(sum(mat,1));
    b = cell(size(colors));
    bar_offsets = linspace(-0.25,0.25,numel(cond_lab));
    for c_ix = 1:numel(colors)
        b{c_ix} = bar(bar_offsets+c_ix,diag(tmp_mat(c_ix,:)),0.9,'stacked');
        for cond_ix = 1:numel(cond_lab)
            set(b{c_ix}(cond_ix),'FaceColor',cond_colors{cond_ix},'EdgeColor','k');
        end
    end
    set(gca,'XTick',1:numel(colors));
    set(gca,'XTickLabels',colors);
    legend(b{1},cond_lab,'Location','best');
    
    % Word-Color bar plot
    subplot(3,1,3); hold on;
    tmp_mat = reshape(mat,[size(mat,1)*size(mat,2) size(mat,3)]);
    tmp_lab = reshape(wc_lab_mat,[size(wc_lab_mat,1)*size(wc_lab_mat,2) size(wc_lab_mat,3)]);
    b = cell(size(tmp_lab));
    bar_offsets = linspace(-0.25,0.25,numel(cond_lab));
    for wc_ix = 1:numel(tmp_lab)
        b{wc_ix} = bar(bar_offsets+wc_ix,diag(tmp_mat(wc_ix,:)),0.9,'stacked');
        for cond_ix = 1:numel(cond_lab)
            set(b{wc_ix}(cond_ix),'FaceColor',cond_colors{cond_ix},'EdgeColor','k');
        end
    end
    set(gca,'XTick',1:numel(tmp_lab));
    set(gca,'XTickLabels',tmp_lab);
    legend(b{1},cond_lab,'Location','best');
end

end