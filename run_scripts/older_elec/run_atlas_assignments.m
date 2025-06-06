% Full SBJ list (completed and ready to run)
SBJs = {'CP24','CP26','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};

for s = 1:numel(SBJs)
%     fn_compile_elec_atlas(SBJs{s},'main_ft','pat','','Dx');
    fn_copy_elec_atlas_pat2mni(SBJs{s},'main_ft','v','Dx');
end
%     fn_save_elec_atlas(SBJs{s},'main_ft','pat','','DK');
%     fn_save_elec_atlas(SBJs{s},'main_ft','pat','','Dx');
%     fn_save_elec_tissue(SBJs{s},'main_ft','pat','','Dx');
%     fn_save_elec_tissue(SBJs{s},'main_ft','pat','','DK');
%     fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo7');
%     fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo17');

%% MNI Check
SBJ = 'IR21';

% Compare patient and MNI in ortho
fn_view_recon(SBJ,'main_ft','ortho','pat','',1,'b');
fn_view_recon(SBJ,'main_ft','ortho','mni','v',1,'b');

% Check atlas assignments
fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','DK','gROI');
fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','Dx','gROI');
fn_view_recon_atlas(SBJ,pipeline_id,'mni','v',1,'b','Yeo7','Yeo7');

