function mediation_dream_combineresults_wani(basedir, mediation_setupfile)

% function mediation_dream_combineresults_wani(model_save_dir, mediation_setupfile)

load(mediation_setupfile);
global j;

for i = 1:size(j,1)
    data{i} = sprintf('res_%s_%d.mat\n', scriptname, i); % j from mediation_SETUP.mat
end

for i = 1:modeln % modeln from mediation_SETUP.mat
    eval(['model' num2str(i) '_dat = model' num2str(i) ';']);
    eval(['if any(~cellfun(@isempty, strfind(fields(model' num2str(i) '), ''L2M''))), model' num2str(i) '_dat = rmfield(model' num2str(i) '_dat, ''L2M''); end;']);
    eval(['if any(~cellfun(@isempty, strfind(fields(model' num2str(i) '), ''L2M''))), model' num2str(i) '_L2M = model' num2str(i) '.L2M; end;']); 
end

check_indices = [];

dat_res = {'p', 'beta', 'ste'};

for i = 1:numel(data)
    load(fullfile(basedir, deblank(data{i})));
    
    last_idx = max(last_idx, size(model1.p,1));
    check_indices(end+1:end+2) = [first_idx last_idx];
    
    for ii = 1:modeln
        for jj = 1:numel(dat_res)
            eval(['model' num2str(ii) '_dat.' dat_res{jj} '(first_idx:last_idx,:) = model' num2str(ii) '.' dat_res{jj} '(first_idx:last_idx,:);']);
            L2Mdat = ['model' num2str(ii) '_L2M.' dat_res{jj} '(first_idx:last_idx,:) = model' num2str(ii) '.L2M.' dat_res{jj} '(first_idx:last_idx,:);'];
            eval(['if any(~cellfun(@isempty, strfind(fields(model' num2str(ii) '), ''L2M''))), eval(L2Mdat); end']);
        end
        
        firstlevel = ['model' num2str(ii) '_dat.firstlevelpaths(first_idx:last_idx,:,:) = model' num2str(ii) '.firstlevelpaths(first_idx:last_idx,:,:);'];
        eval(['if any(~cellfun(@isempty, strfind(fields(model' num2str(ii) '), ''firstlevelpaths''))), eval(firstlevel); end']);
    end
end

%% write result images

for i = 1:modeln
    modeldir{i} = fullfile(basedir, models.name{i}); 
    if ~exist(modeldir{i}), mkdir(modeldir{i}); end
end

pathname_twomed = {'X-M', 'M-Y', 'X-Y_direct', 'X-Y_total', 'X-M-Y'};
pathname_threemed = {'X-M1', 'M1-M2', 'M2-Y', 'X-Y_diect', 'X-Y_total', 'X-M1-M2-Y'};
compname = {'_pvals', '_effect', '_ste'};

data = fmri_data(med_vars.imgs{1}(1,:), mask);

strfind(models.fns{1}, 'mediation_three')

for i = 1:modeln
    cd(modeldir{i});
    twomed = strfind(models.fns{i}, 'mediation(');
    threemed = strfind(models.fns{i}, 'mediation_threepaths(');
    field_name = {'p', 'beta', 'ste'};
    % eval(['field_name = fields(model' num2str(i) '_dat);']);
    eval(['isL2M = any(~cellfun(@isempty, strfind(fields(model' num2str(i) '), ''L2M'')));']);
    eval(['isfirstlevel = any(~cellfun(@isempty, strfind(fields(model' num2str(i) '), ''firstlevelpaths'')));']);
    
    for ii = 1:numel(models.savepaths{i})
        for jj = 1:numel(field_name)
            eval(['data.dat = model' num2str(i) '_dat.' field_name{jj} '(:,' num2str(ii) ');']);
            if twomed == 1
                data.fullpath = fullfile(modeldir{i}, [pathname_twomed{models.savepaths{i}(ii)} compname{jj} '.nii']);
            elseif threemed == 1
                data.fullpath = fullfile(modeldir{i}, [pathname_threemed{models.savepaths{i}(ii)} compname{jj} '.nii']);
            end
            write(data);
            
            if isL2M
                eval(['data.dat = model' num2str(i) '_L2M.' field_name{jj} '(:,' num2str(ii) ');']);
                if twomed == 1
                    data.fullpath = fullfile(modeldir{i}, [pathname_twomed{models.savepaths{i}(ii)} '_L2M' compname{jj} '.nii']);
                elseif threemed == 1
                    data.fullpath = fullfile(modeldir{i}, [pathname_threemed{models.savepaths{i}(ii)} '_L2M' compname{jj} '.nii']);
                end
                write(data);
            end
            
        end
        
        if isfirstlevel
            
            eval(['data.dat = squeeze(model' num2str(i) '_dat.firstlevelpaths(:,' num2str(ii) ',:));']);
            if twomed == 1
                data.fullpath = fullfile(modeldir{i}, ['firstlevel_' pathname_twomed{models.savepaths{i}(ii)} 'betas.nii']);
            elseif threemed == 1
                data.fullpath = fullfile(modeldir{i}, ['firstlevel_' pathname_threemed{models.savepaths{i}(ii)} 'betas.nii']);
            end
            write(data);
            
        end
    end
end

end