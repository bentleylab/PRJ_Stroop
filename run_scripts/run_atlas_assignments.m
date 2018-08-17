SBJs = {'IR21','IR27','IR31','IR35','IR39','IR41','IR48','IR52','IR54','IR57','IR61','IR68','IR72','IR74','IR62'};
% SBJs = {'IR65','IR67','IR69'};

for s = 1:numel(SBJs)
    fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo7');
    fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo17');
%     fn_save_elec_atlas(SBJs{s},'main_ft','pat','','DK');
%     fn_save_elec_atlas(SBJs{s},'main_ft','pat','','Dx');
end
% 
% fn_compile_elec_struct('CP24','main_ft','pat','');
% fn_compile_elec_struct('CP24','main_ft','mni','v');
fn_compile_elec_struct('CP24','main_ft','mni','s');

fn_compile_elec_struct('IR62','main_ft','pat','');
fn_compile_elec_struct('IR62','main_ft','mni','v');

fn_compile_elec_struct('IR67','main_ft','pat','');
fn_compile_elec_struct('IR67','main_ft','mni','v');

fn_compile_elec_struct('IR69','main_ft','pat','');
fn_compile_elec_struct('IR69','main_ft','mni','v');