files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

SBJs = cellfun(@(x) x(1:end-4), files, 'UniformOutput', false);

clear files

roi_opts  = {{'l','MPFC',1},{'l','lat',1}};
custom_ROIs = {'dmPFC','dlPFC'};
proc_id   = 'main_ft';
atlas_id  = 'Dx';
reg_type  = 'v';
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'fig';
roi_id = 'ROI';
%%
for roi_ix = 1:numel(roi_opts)
    test_grp_ROI(SBJs, proc_id, reg_type, show_lab,...
                                roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
                                roi_opts{roi_ix}{3},custom_ROIs{roi_ix},'save_fig', save_fig, 'fig_ftype', fig_ftype);
end