%% Quick and dirty R/L labels for elec files
% NOTE: only corrections where automatic detection based on label didn't
%   work were IR68 (all L side implant with AIN, MIN, PIN). Another excpetion
%   was IR32, but it was all R sided except IHL, so automatic labeling was
%   already correct.
pipeline_id = 'main_ft';
show_labels = 1;
hemi        = 'b';
atlas_ids   = {'Dx','Yeo7'};
roi_ids = {'gROI','Yeo7'};
view_spaces  = {'pat','mni'};
reg_types    = {'','_v'};

%%
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};%ALMOST: 'CP24','IR26',  %NEVER: 'IR27','IR37','IR48',
for s = 1:numel(SBJs)
    for atlas_ix = 1:2
        SBJ = SBJs{s};
        SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
        eval(SBJ_vars_cmd);
        fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_spaces{atlas_ix},...
            reg_types{atlas_ix},'_',atlas_ids{atlas_ix},'.mat'];
%         fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_spaces{atlas_ix},...
%             reg_types{atlas_ix},'.mat'];
        load(fname);
        
        % Now do the manual labor
        elec.hemi = repmat({'r'},size(elec.label));
        for e = 1:numel(elec.label)
            if strfind(elec.label{e},'L')
                elec.hemi{e} = 'l';
            end
        end
        
        % Fix IR68 (only exception)
        if strcmp(SBJ,'IR68')
            elec.hemi = repmat({'l'},size(elec.label));
        end
        
        %% Save
        save(fname,'-v7.3','elec');
        clear SBJ_vars SBJ_vars_cmd
    end
end