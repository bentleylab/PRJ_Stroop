% Full SBJ list (completed and ready to run)
SBJs = {'CP24','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};

for s = 13:numel(SBJs)
    % Convert raw pipeline elec files to my SBJ_vars
    fn_convert_elec_struct(SBJs{s},'main_ft','pat','',0);
%     fn_convert_elec_struct(SBJs{s},'main_ft','pat','',1);
    
    % Match elec to atlas labels + tissue (ONLY orig!)
    % run in SGE: fn_match_elec_atlas_ROI_tiss(SBJs{s},'main_ft','pat','','Dx',0);
    
    % Compile orig elec atlas match to reref
%     fn_compile_elec_atlas(SBJs{s},'main_ft','pat','','Dx');
    
    % Export reref atlas info to CSV for manual adjustments
%     fn_export_elec2csv(SBJs{s},'main_ft','pat','','Dx');
    
    % Manual adjustments
    % fn_view_elec_ROI_assignments();
    
    % Reimport ???
    
    % Copy corrected labels to MNI elec files
    % check this: fn_copy_elec_atlas_pat2mni(SBJs{s},'main_ft','v','Dx');
end

%% MNI Check ???
% SBJ = 'IR21';
% 
% % Compare patient and MNI in ortho
% fn_view_recon(SBJ,'main_ft','ortho','pat','',1,'b');
% fn_view_recon(SBJ,'main_ft','ortho','mni','v',1,'b');
% 
% % Check atlas assignments
% fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','DK','gROI');
% fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','Dx','gROI');
% fn_view_recon_atlas(SBJ,pipeline_id,'mni','v',1,'b','Yeo7','Yeo7');

