% previous main list: SBJs = {'IR54','IR57','IR61','IR68','IR72','IR74','IR62'};%'IR52','IR21','IR27','IR31','IR35','IR39','IR41','IR48',
% tissue not 100%: 
% need work: SBJs = {'IR67','IR69'};
%   never! 'IR62',

% SBJs list for CPPC 2018 abstract
% SBJs = {'IR21','IR27','IR31','IR35','IR39','IR41','IR48','IR52','IR54','IR57','IR61','IR68','IR72','IR74'};

% Full SBJ list (completed and ready to run)
SBJs = {'CP24','IR21','IR27','IR31','IR35','IR39','IR41','IR48','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};

SBJs = {'IR32'};
fn_compile_elec_struct(SBJs{1},'main_ft','pat','');
fn_compile_elec_struct(SBJs{1},'main_ft','mni','v');

for s = 1:numel(SBJs)
    fn_save_elec_atlas(SBJs{s},'main_ft','pat','','DK');
    fn_save_elec_atlas(SBJs{s},'main_ft','pat','','Dx');
    fn_save_elec_tissue(SBJs{s},'main_ft','pat','','Dx');
    fn_save_elec_tissue(SBJs{s},'main_ft','pat','','DK');
    fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo7');
    fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo17');
end

fn_compile_elec_struct('IR67','main_ft','pat','');
fn_compile_elec_struct('IR67','main_ft','mni','v');

fn_compile_elec_struct('IR69','main_ft','pat','');
fn_compile_elec_struct('IR69','main_ft','mni','v');