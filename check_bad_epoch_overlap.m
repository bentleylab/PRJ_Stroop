%% Compare preproc and preclean bad_epochs
SBJs = {'CP26','IR21','IR26','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};
pipeline_id = 'main_ft';

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
% addpath([root_dir 'PRJ_Stroop/scripts/']);
% addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
% addpath(ft_dir);
% ft_defaults

for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %%
    pc = fn_compile_epochs_full2at(SBJ,pipeline_id);
    pc_ix = cell([size(pc,1) 1]);
    for e = 1:numel(pc_ix)
        pc_ix{e} = pc(e,1):pc(e,2);
    end
    
    load([SBJ_vars.dirs.events SBJ '_colin_bad_epochs_preproc.mat']);
    pp = bad_epochs;
    
    pp_ix = cell([size(pp,1) 1]);
    for e = 1:numel(pp_ix)
        pp_ix{e} = pp(e,1):pp(e,2);
    end
    pp_all = cat(2,pp_ix{:});
    
    %%
    overlap = zeros(size(pc_ix));
    for e = 1:numel(pc_ix)
        if ~isempty(union(pc_ix{e},pp_all))
            overlap(e) = 1;
        end
    end
    if sum(overlap)~=numel(overlap)
        fprintf(2,'%s = %i / %i\n',SBJ,sum(overlap),numel(overlap));
    else
        fprintf('%s OK\n',SBJ);
    end
    clear SBJ SBJ_vars pc pp overlap
end