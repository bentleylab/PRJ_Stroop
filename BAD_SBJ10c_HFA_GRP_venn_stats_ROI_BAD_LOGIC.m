function SBJ10c_HFA_GRP_venn_stats_ROI(SBJs, proc_id, stat_conds, hemi, atlas_id, roi_id,...
                                       plot_out, plt_id, save_fig, fig_ftype)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJs [cell array str] - subject IDs to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_conds [cell array] - {{'stat_id1','an_id1','cond1'},...,{'stat_idn','an_idn','condn'}};
%     stat_id_n [str] - ID of statistical analysis
%     an_id [str] - analysis ID for preprocessing, filtering, etc.
%     cond_n [str] - ID of condition/group to overlap with
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pCNI': ANOVA of previous trial congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT
%   plt_id [str] - ID of the plotting variables (likely 'venn')
%   plot_out [0/1] - include electrodes that don't have an atlas label or in hemi?
%   save_fig [0/1] - save this figure?
%   fig_ftype [str] - file extension for figure saving

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Define default options
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Organize IDs
stat_ids = cell(size(stat_conds)); an_ids = cell(size(stat_conds));
cond_ids = cell(size(stat_conds)); ep_labs = cell(size(stat_conds));
for st_ix = 1:numel(stat_conds)
    stat_ids{st_ix} = stat_conds{st_ix}{1};
    an_ids{st_ix}   = stat_conds{st_ix}{2};
    cond_ids{st_ix} = stat_conds{st_ix}{3};
end
if numel(stat_conds) < 2 || numel(stat_conds) > 3; error('why venn?'); end
for st_ix = 1:numel(stat_conds)
    if numel(stat_conds{st_ix})~=3; error('need stat_id, an_id, cond in each stat_cond'); end
end

% Venn colors
if strcmp(an_ids{1},'HGm_S2t151_zbtA_sm0_wn100') && ...
        strcmp(an_ids{2},'HGm_S2t251_zbtA_sm0_wn100') && ...
        strcmp(an_ids{3},'HGm_R5t101_zbtA_sm0_wn100') && ...
        strcmp(cond_ids{1},cond_ids{2}) && strcmp(cond_ids{2},cond_ids{3})
    color_id = 'SDR';
elseif strcmp(cond_ids{1},'CNI') && strcmp(cond_ids{2},'pCNI') && strcmp(cond_ids{3},'PC')
    color_id = 'ConCSPC';
elseif strcmp(cond_ids{1},'pCNI') && strcmp(cond_ids{2},'PC')
    color_id = 'CSPC';
end
if exist('color_id','var')
    venn_colors = fn_venn_colors(numel(stat_conds), 'cond_id', color_id);
else
    venn_colors = fn_venn_colors(numel(stat_conds));
end
venn_colors = squeeze(venn_colors(logical(eye(size(venn_colors,1)))));
% venn_colors = cell(size(stat_conds));
% venn_colors{1} = plt_vars.cmap(1,:);
% if numel(stat_conds)==3;
%     venn_colors{2} = plt_vars.cmap(50,:);
%     venn_colors{3} = plt_vars.cmap(end,:);
% else
%     venn_colors{2} = plt_vars.cmap(end,:);
% end
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

%% Prep report
% Create output dir
save_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_overlap/' strjoin(stat_ids,'-')...
            '/' strjoin(an_ids,'-') '/' strjoin(cond_ids,'-') '/'];
if ~exist(save_dir,'dir')
    [~] = mkdir(save_dir);
end
% Create report
sig_report_fname = [save_dir 'GRP_' atlas_id '_' roi_id '_' hemi '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');
% Print IDs
title_str = repmat('%-40s',[1 numel(stat_conds)]);
fprintf(sig_report,[title_str '\n'],stat_ids{:});
fprintf(sig_report,[title_str '\n'],an_ids{:});
fprintf(sig_report,[title_str '\n'],cond_ids{:});
fprintf(sig_report,[repmat('=',[1 40*numel(stat_conds)]) '\n']);
% Print result headings
fprintf(sig_report,[repmat('%-20s',[1 numel(stat_conds)+2]) '\n'],'label','ROI',cond_ids{:});
fprintf(sig_report,[repmat('-',[1 40*numel(stat_conds)]) '\n']);
result_str = ['%-20s%-20s' repmat('%-20i',[1 numel(stat_conds)]) '\n'];

%% Load Data
elec_sbj = cell([numel(SBJs) 1]);
good_sbj = true([numel(SBJs) 1]);
sig_mat  = cell([numel(SBJs) 1]);
roi_mat  = cell([numel(SBJs) 1]);
% all_roi_labels = cell([numel(stat_conds) 1]);
% all_roi_colors = cell([numel(SBJs) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% Prepare elec structs
    % Load elec struct
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat'];
    tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    
    % Append SBJ name to labels
    orig_labels = elec_sbj{sbj_ix}.label;
    for e_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{e_ix}];
        % don't need individual elec colors here
    end
    
    % Remove hemi and/or atlas elecs
    if ~plot_out
        % Remove electrodes that aren't in atlas ROIs & hemisphere
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
    else
        % Remove electrodes that aren't in hemisphere
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, [], []);
    end
    cfgs = []; cfgs.channel = roi_elecs;
    elec_sbj{sbj_ix} = fn_select_elec(cfgs,elec_sbj{sbj_ix});
    
    %% Load Stats
    if ~isempty(roi_elecs)
        sig_mat{sbj_ix} = zeros([numel(elec_sbj{sbj_ix}.label) numel(stat_conds)]);
        roi_mat{sbj_ix} = zeros([numel(elec_sbj{sbj_ix}.label) 1]);
        for st_ix = 1:numel(stat_conds)
            if strcmp(cond_ids{st_ix},'actv')
                error('not done actv yet');
                % load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_' stat_id '.mat'],'actv');
            elseif strcmp(cond_ids{st_ix},'RT')
                error('RT not implemented');
            else    % ANOVA
                load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_ids{st_ix} '_' an_ids{st_ix} '.mat']);
                ep_labs{st_ix} = st.ep_lab;
                if iscell(ep_labs{st_ix}); ep_labs{st_ix} = [ep_labs{st_ix}{:}]; end % unpack multi-epoch analyses
                for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
                    stat_ch_ix = strcmp(w2.label,orig_labels{ch_ix});
                    if ~any(stat_ch_ix)
                        fprintf(2,'\tWARNING: elec chan %s not found in w2!\n',elec_sbj{sbj_ix}.label{ch_ix});
                    end
                    % Consolidate to binary sig/non-sig
                    if any(squeeze(w2.qval(strcmp(st.groups,cond_ids{st_ix}),stat_ch_ix,:))<st.alpha)
                        sig_mat{sbj_ix}(ch_ix,st_ix) = 1;
                    end
                end
                clear w2 st actv
            end
        end
        
        % Report on significant electrodes for this SBJ (even if non-sig)
        for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
            fprintf(sig_report,result_str,elec_sbj{sbj_ix}.label{ch_ix},...
                                          elec_sbj{sbj_ix}.(roi_field){ch_ix},sig_mat{sbj_ix}(ch_ix,:));
            % Grab roi_ix
            roi_mat{sbj_ix}(ch_ix) = find(strcmp(roi_list,elec_sbj{sbj_ix}.(roi_field){ch_ix}));
        end
        
%         %% Compile Statistics
%         if any(sig_mat{sbj_ix}(:))
%             print_nums = zeros([2 numel(stat_conds)]);
%             for st_ix = 1:numel(stat_conds)
%                 print_nums(1,st_ix) = sum(sig_mat{sbj_ix}(:,st_ix));
%                 print_nums(2,st_ix) = sum(sig_mat{sbj_ix}(:,st_ix))/size(sig_mat{sbj_ix},1);
%             end
%             fprintf(['%-10s' repmat('%-3i(%.3f)\t',[1 numel(stat_conds)]) '\n'],[SBJ ' sig:'],print_nums(:));
%             
%             % Select color and add to elec
%             all_roi_colors{sbj_ix} = zeros([sum(any(sig_mat{sbj_ix},2)) 3]);
%             sig_elecs = cell([sum(any(sig_mat{sbj_ix},2)) 1]);
%             sig_ix = 0;
%             for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
%                 if any(sig_mat{sbj_ix}(ch_ix,:))
%                     sig_ix = sig_ix + 1;
%                     sig_elecs{sig_ix} = elec_sbj{sbj_ix}.label{ch_ix};
%                     % Select color
%                     cond_ix = find(sig_mat{sbj_ix}(ch_ix,:));
%                     if numel(cond_ix)==numel(stat_conds) && numel(stat_conds) == 3
%                         all_roi_colors{sbj_ix}(sig_ix,:) = all_color;
%                     elseif numel(cond_ix)==2
%                         all_roi_colors{sbj_ix}(sig_ix,:) = venn_colors{cond_ix(1),cond_ix(2)};
%                     else
%                         all_roi_colors{sbj_ix}(sig_ix,:) = venn_colors{cond_ix,cond_ix};
%                     end
%                 end
%             end
%             
%             % Select sig elecs for plotting
%             cfgs = []; cfgs.channel = sig_elecs;
%             elec_sbj{sbj_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix});
%         else
%             % Print no significant elecs
%             elec_sbj{sbj_ix} = {}; good_sbj(sbj_ix) = false;
%             fprintf(2,'\t%s has no significant channels for any stat\n',SBJ);
%         end
        
    else
        % Print no ROI match
        elec_sbj{sbj_ix} = {}; good_sbj(sbj_ix) = false;
        fprintf(2,'\t%s has no channels in %s hemi %s\n',SBJ,atlas_id,hemi);
    end
    clear SBJ SBJ_vars SBJ_vars_cmd sig_elecs sig_ix cond_ix 
end

% Finish report and save sig_mat
fclose(sig_report);
% sig_elecs_fname = [save_dir 'GRP_' atlas_id '_' roi_id '_sig_elecs.mat'];
% save(sig_elecs_fname,'-v7.3','stat_conds','SBJs','elec_sbj','sig_mat');

%% Combine elec structs
elec = ft_appendsens([],elec_sbj{good_sbj});
% elec.roi_color = vertcat(all_roi_colors{:});    % appendsens strips that field
% elec.(roi_field)       = all_roi_labels{cond_ix};    % appendsens strips that field

sig_all = vertcat(sig_mat{:});
sig_all = sig_all(any(sig_all,2),:);
roi_all = vertcat(roi_mat{:});
roi_all = roi_all(any(sig_all,2));

%% Compute Stat Overlaps
% Find number of comparisons
pairs = nchoosek(1:numel(stat_conds),2);
n_groups = numel(stat_ids)+size(pairs,1);
if numel(stat_conds)==3; n_groups = n_groups+1; end

sig_cnt     = zeros([n_groups 1]);
sig_type    = zeros([n_groups 1]);
sig_grps    = cell([n_groups 1]);
sig_cnt_roi = zeros([n_groups numel(roi_list)]);
venn_legend = cell(size(stat_conds));
% roi_legends = cell([numel(stat_conds) numel(roi_list)]);
% Add main effects
sig_ix = 0;
for main_ix = 1:numel(stat_conds)
    sig_ix = sig_ix + 1;
    sig_type(sig_ix) = 1;
    sig_cnt(sig_ix) = sum(sig_all(:,main_ix));
    sig_grps{sig_ix} = main_ix;
    venn_legend{sig_ix} = [cond_ids{main_ix} '-' ep_labs{main_ix}];
    for roi_ix = 1:numel(roi_list)
        sig_cnt_roi(sig_ix,roi_ix) = sum(sig_all(roi_all==roi_ix,main_ix));
    end
end

% Add pairs
for p_ix = 1:size(pairs,1)
    sig_ix = sig_ix + 1;
    sig_type(sig_ix) = 2;
    sig_grps{sig_ix} = pairs(p_ix,:);
    sig_cnt(sig_ix) = sum(sig_all(:,pairs(p_ix,1))==1 & sig_all(:,pairs(p_ix,2))==1);
%     venn_legend{sig_ix} = [cond_ids{pairs(p_ix,1)} '-' ep_labs{pairs(p_ix,1)} ' + '...
%         cond_ids{pairs(p_ix,2)} '-' ep_labs{pairs(p_ix,2)} '(n=' num2str(sig_cnt(sig_ix)) ')'];
    for roi_ix = 1:numel(roi_list)
        sig_cnt_roi(sig_ix,roi_ix) = sum(sig_all(roi_all==roi_ix,pairs(p_ix,1))==1 & ...
                                  sig_all(roi_all==roi_ix,pairs(p_ix,2))==1);
    end
end

% Add triples
if numel(stat_conds)==3
    sig_cnt(end) = sum(all(sig_all,2));
    sig_type(end) = 3;
    sig_grps{end} = 1:3;
%     venn_legend{end} = ['ALL (n=' num2str(sig_cnt(end)) ')'];
    for roi_ix = 1:numel(roi_list)
        sig_cnt_roi(end,roi_ix) = sum(all(sig_all(roi_all==roi_ix,:),2));
    end
end

%% Plot Overall Venn Diagram
fig_name = ['GRP_venn_' atlas_id '_' roi_id '_' strjoin(stat_ids,'-')];
if numel(stat_conds)==3; pos = plt_vars.fig_pos{2};
else pos = plt_vars.fig_pos{1};
end
f = figure('Name',fig_name,'units','normalized','outerposition',pos);

% f = venn(Z), where:
%   n=2: [z1 z2 z12]
%   n=3: [z1 z2 z3 z12 z13 z23 z123]
[v, v_stats] = venn(sig_cnt,'FaceColor',venn_colors);
for z_ix = 1:numel(sig_cnt)
    text(v_stats.ZoneCentroid(z_ix,1),v_stats.ZoneCentroid(z_ix,2),...
         num2str(sig_cnt(z_ix)),'FontSize',plt_vars.text_sz,'HorizontalAlignment','center');
end
axis(gca, plt_vars.axis_vis);
title([atlas_id '-' roi_id ': # SBJ = ' num2str(numel(SBJs)) ', # Sig Elecs = '...
       num2str(size(sig_all,1)) ' / ' num2str(size(vertcat(sig_mat{:}),1)) ...
       ' (' num2str(size(sig_all,1)/size(vertcat(sig_mat{:}),1),'%.4f') ')'],...
       'FontSize',plt_vars.title_sz);
legend(v,venn_legend,'FontSize',plt_vars.legend_sz);

if save_fig
    fig_fname = [save_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Venn Diagrams per ROI
for roi_ix = 1:numel(roi_list)
    fig_name = ['GRP_venn_' atlas_id '_' roi_id '_' strjoin(stat_ids,'-') '_' roi_list{roi_ix}];
    f = figure('Name',fig_name,'units','normalized','outerposition',pos);
    ax = gca; hold on;
    
    if ~any(sig_cnt_roi(:,roi_ix)==0)
        % f = venn(Z), where:
        %   n=2: [z1 z2 z12]
        %   n=3: [z1 z2 z3 z12 z13 z23 z123]
        [v, v_stats] = venn(sig_cnt_roi(:,roi_ix),'FaceColor',venn_colors);
        text_spacers = zeros(size(v_stats.ZoneCentroid));
    else % Not all circles overlap, venn fails
        % Plot circles for each method
        [~,max_idx] = sort(sig_cnt_roi(:,roi_ix),'descend');
        scat_space = sum(sig_cnt_roi(max_idx(1:2),roi_ix))*cos(pi/4);
        scat_locs  = [-scat_space 0; scat_space 0; 0 scat_space];
        v = gobjects(size(stat_conds));
        v_stats.ZoneCentroid = zeros([n_groups 2]);
        for st_ix = 1:numel(stat_conds)
            v(st_ix) = patch(sig_cnt_roi(st_ix,roi_ix)*cos(plt_vars.ang)+scat_locs(st_ix,1),...
                             sig_cnt_roi(st_ix,roi_ix)*sin(plt_vars.ang)+scat_locs(st_ix,2),...
                             venn_colors{st_ix});
            set(v(st_ix),'FaceAlpha',plt_vars.alpha);
            v_stats.ZoneCentroid(st_ix,:) = scat_locs(st_ix,:);
        end
        ax.XLim = [-scat_space-max(sig_cnt_roi(:,roi_ix))-plt_vars.ax_fudge ...
                    scat_space+max(sig_cnt_roi(:,roi_ix))+plt_vars.ax_fudge];
        ax.YLim = [-scat_space-max(sig_cnt_roi(:,roi_ix))-plt_vars.ax_fudge ...
                    scat_space+max(sig_cnt_roi(:,roi_ix))+plt_vars.ax_fudge];
        
        % Plot lines to show overlap
        grp_idx = find(sig_type==2);
        for p_ix = 1:size(pairs,1)
            grp_ix = grp_idx(p_ix);
            if sig_cnt_roi(grp_ix,roi_ix)~=0
                xs = [scat_locs(sig_grps{grp_ix}(1),1) scat_locs(sig_grps{grp_ix}(2),1)];
                ys = [scat_locs(sig_grps{grp_ix}(1),2) scat_locs(sig_grps{grp_ix}(2),2)];
                line(xs, ys, 'LineWidth',sig_cnt_roi(grp_ix,roi_ix),'Color',plt_vars.line_color);
                v_stats.ZoneCentroid(grp_ix,:) = [mean(xs) mean(ys)];
            end
        end
        if numel(stat_conds)==3
            v_stats.ZoneCentroid(end,:) = mean(scat_locs,1);
            text_spacers = plt_vars.text_spacers{2};
        else
            text_spacers = plt_vars.text_spacers{1};
        end
        if ~any(sig_cnt_roi(grp_idx,roi_ix))
            text_spacers = zeros(size(text_spacers));
        end
    end
    
    % Plot text for numbers of electrodes in each zone
    for z_ix = 1:n_groups
        text(v_stats.ZoneCentroid(z_ix,1)+text_spacers(z_ix,1),...
             v_stats.ZoneCentroid(z_ix,2)+text_spacers(z_ix,2),...
            num2str(sig_cnt_roi(z_ix,roi_ix)),'FontSize',plt_vars.text_sz,...
            'HorizontalAlignment',plt_vars.text_align);
    end
    axis equal; axis(gca, plt_vars.axis_vis);
    title([atlas_id '-' roi_id '-' roi_list{roi_ix} ': # SBJ = ' num2str(numel(SBJs)) ...
        ', # Sig Elecs = ' num2str(sum(roi_all==roi_ix)) ' / ' num2str(sum(vertcat(roi_mat{:})==roi_ix)) ...
       ' (' num2str(sum(roi_all==roi_ix)/sum(vertcat(roi_mat{:})==roi_ix),'%.4f') ')'],...
       'FontSize',plt_vars.title_sz,'Color',roi_colors{roi_ix});
    legend(v,venn_legend,'FontSize',plt_vars.legend_sz);
    
    if save_fig
        fig_fname = [save_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
