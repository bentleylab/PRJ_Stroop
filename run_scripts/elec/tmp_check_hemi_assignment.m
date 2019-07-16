for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Original elec files
    elec_fname = SBJ_vars.recon.elec_pat;
    slash = strfind(elec_fname,'/'); elec_suffix = elec_fname(slash(end)+numel(SBJ)+2:end-4);
    
    tmp = load(elec_fname);
    elec_var_name = fieldnames(tmp);
    if ~strcmp(elec_var_name,elec_suffix)
        warning(['\t!!!! ' SBJ ' elec names in variable and file names do not match! file=' elec_suffix '; var=' elec_var_name{1}]);
    end
    eval(['elec = tmp.' elec_var_name{1} ';']); clear tmp;
    
    % Print default logic
    elec.hemi = repmat({'r'},size(elec.label));
    for e = 1:numel(elec.label)
        if strfind(elec.label{e},'L')
            elec.hemi{e} = 'l';
        end
    end
    horzcat(elec.label,elec.hemi)
    
    % Print warnings for failure cases
    l_str = ~cellfun(@isempty,strfind(elec.label,'L'));
    r_str = ~cellfun(@isempty,strfind(elec.label,'R'));
    if any(all([l_str r_str],2))
        both_ix = find(all([l_str r_str],2));
        for e = 1:numel(both_ix)
            fprintf(2,'%s both in :%s\n',SBJ,elec.label{both_ix(e)});
        end
    end
    if any(~any([l_str r_str],2))
        neither_ix = find(~any([l_str r_str],2));
        for e = 1:numel(neither_ix)
            fprintf(2,'%s neither in :%s\n',SBJ,elec.label{neither_ix(e)});
        end
    end
    pause;
end