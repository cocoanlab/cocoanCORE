function [dat, parcel_list] = parcel_view(varargin)
% show orthviews of parcellation in SPM
% 'Fan_et_al_atlas_r280' or 'HarvardOxford parcels'
%
% example ::  [dat, parcel_list] = parcel_view('fan_280');
%         ::  [dat, parcel_list] = parcel_view('harvardoxford');
%
% Copyright Byeol Kim, 2019
% 
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'fan_280'}
                version = 1;
            case {'harvardoxford'}
                version = 2;
        end
    end
end

idx_list{1} = [0 14 28 40 52 64 68 80 88 102 108 120 124 134 146 154 162 174 188 198 210 214 218 230 246 273 274 275 276 278 280];
idx_list{2} = [0:4 6:8 10 13 16 17 18 20 21 23:28 30:33 35 36 39:48];


parcel{1} = which('Fan_et_al_atlas_r280.nii');
parcel{2} = which('HarvardOxford-cort-maxprob-thr0-1mm.nii.gz');


if ~isempty(parcel{version})
    
    filepath = fileparts(parcel{version});
    
    if version == 1        
        short_names = {'sfg','mfg','ifg','org','prg','pcl','stg','mtg','itg','fug','phg','psts',...
            'spl','ipl','pcun','pog','ins','cg','mvocc','locc','amg','hip','bg','tha',...
            'cb','hyp/tha','brainstem','pag','s_nigra','reg_nucleus'};
        
        prompt = sprintf('\n:::::: Choose one parcellation among ::::::\n  ** %s **  \n  :  ', strjoin(short_names));
        while 1
            in_str = input(prompt,'s');
            if ~sum(ismember(short_names, in_str))
                warning('That is improper input.')
            else
                break
            end
        end
        
        dat = fmri_data(parcel{version});
        
        prc_n = find(ismember(short_names, in_str));
        
        parcel_list = load(fullfile(filepath, 'cluster_Fan_Net_r280.mat'));
        parcel_list = parcel_list.cluster_Fan_Net.names_short;

    elseif version == 2
        
        short_names = {'frontal pole','insular cortex','superior frontal gyrus','middle frontal gyrus','inferior frontal gyrus',...
            'precentral gyrus','temporal pole', 'superior temporal gyrus','middle temporal gyrus','inferior temporal gyrus',...
            'postcentral gyrus', 'superior parietal lobule', 'supramarginal gyrus','angular gyrus','lateral occipital cortex',...
            'intracalcarine cortex','frontal medial cortex', 'supplementary motor cortex','subcallosal cortex', ...
            'paracingulate gyrus','cingulate gyrus','precuneous cortex',...
            'cuneal cortex', 'frontal orbital cortex', 'parahip','lingual gyrus','temporal fusiform', 'occipital fusiform', ...
            'frontal operculum cortex', 'central opercular cortex','parietal operculumn cortex', 'planum polare', 'Heschls gyrus', ...
            'planum temporale', 'supracalcarine cortex', 'occipital pole'};
        
        fprintf('\n:::::: Choose one parcellation *number* among below :::::: \n');
        disp(table(short_names', [1:numel(short_names)]', 'VariableNames',{'parcel','number'}));       
        prompt = sprintf('\n  :  ');
        while 1
            prc_n = input(prompt);
            if isempty(prc_n) || prc_n > numel(short_names) || prc_n < 1
                warning('That is improper input.')
            else
                break
            end
        end
        parcel_list = load(fullfile(filepath, 'regions_name.mat'));
        parcel_list = parcel_list.region_names;
    end
    
    dat = fmri_data(parcel{version}, 'noverbose');
    
    idx = idx_list{version}(prc_n)+1:idx_list{version}(prc_n+1);
    dat.dat(~ismember(dat.dat,idx)) = 0;
    parcel_list = parcel_list(idx);
    
    orthviews(dat);
    
    fprintf('\n ***** parcel "%s" from "%s" is loaded! ***** \n', short_names{prc_n}, varargin{1})
else
    expected_dir{1} = 'github/cocoanlab/cocoanCORE/Canonical_brains/Brainnetome';
    expected_dir{2} = 'github/cocoanlab/cocoanCORE/Canonical_brains/HarvardOxford';
    
    error('Cannot find the parcellation files on saved path. Make sure this directory exists and was added on path. \n   .../%s', expected_dir{version})
end
end
