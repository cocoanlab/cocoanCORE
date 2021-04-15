%%
basedir =  '/sas1/cocoanlab/data/SEMIC/';
addpath(genpath(basedir));
modeldir = fullfile(basedir,'/analysis/imaging/first_level/model03a_SPM_SINGLE_TRIAL_PAIN');
%% LOAD CANLABDATASET
D =  [];
input = 'HPC';
D = SEMIC_load_canlabdata(input);

OvrRating = get_var(D, 'overall_rating', 'conditional', {'EventName','OverallRating'});
OvrType = get_var(D, 'overall_question_type', 'conditional', {'EventName','OverallRating'});
PreRating = get_var(D, 'cont_rating_final', 'conditional', {'EventName','OverallRating'}); 

cond_heat = get_var(D, 'stim_level', 'conditional', {'EventName','HeatStart'});
cond_cue =  get_var(D, 'cue_level', 'conditional', {'EventName', 'HeatStart'});

datCont = D.Sub_Event_Level.data{1}; % Continous prediction rating
datOvr = D.Sub_Event_Level.data{2}; % Overall rating  

secs046 = 0.46/0.02;
secs2_5 = 2.5/0.02;
trial_idx =[];
for run_i=1:6; trial_idx = [trial_idx (1:18)+20*(run_i-1)]; end;
%% ====================================================================== %
%                         FOR SELF AND OTHER RATING                       %
%  ====================================================================== %
clc;
M_self = []; M_other = []; 
Cond_heat_self = [];  Cond_heat_other = [];
Cond_cue_self = [];Cond_cue_other = []; 
FinalPredicRating_self = [];FinalPredicRating_other = [];
for sub_i = [1:16, 18:59]   
    subject_modeldir = fullfile(modeldir, sprintf('sub-semic%03d',sub_i));
    out = [];
    load(fullfile(subject_modeldir,'vifs.mat'),'out');
    
    vifs{sub_i} = out.allvifs(trial_idx);
    betanames{sub_i} = out.name(trial_idx);
    
    sorted_idx =[];
    [~,sorted_idx] = sort(betanames{sub_i});
    % FOR HEAT CONDITION
    heaTName = lower({'LV1','LV2','LV3','LV4','LV5'});
    temp_stim = nan(108,1);
    for heat_cond_i = 1:5
        temp_stim(find(contains(betanames{sub_i}, heaTName{heat_cond_i}))) = heat_cond_i;
    end
    Cond_heat_self{sub_i} = temp_stim;
    Cond_heat_other{sub_i} = temp_stim;
    
    % FOR CUE CONDITION
    cuEName = {'low','no','high'};
    temp_cue = nan(108,1);
    for cue_cond_i = 1:3
        %temp_cue(find(strcmp(semic_info.data.condition.cue{sub_i}, cuEName{cue_cond_i}))) = cue_cond_i -2;
        temp_cue(find(contains(betanames{sub_i}, cuEName{cue_cond_i}))) = cue_cond_i - 2;
    end
    Cond_cue_self{sub_i} = temp_cue;
    Cond_cue_other{sub_i} = temp_cue;
    
    % FOR FINAL PARTS OF CONTINUOUS RATING    
    FinalIntesityRating_self{sub_i} = OvrRating(sub_i,sorted_idx);
    FinalIntesityRating_other{sub_i} = OvrRating(sub_i,sorted_idx);
    
    % Mediator
    temp_M = [];
    temp_M = filenames(fullfile(modeldir,sprintf('sub-semic%03d',sub_i),'beta*'));
    temp_M = temp_M(trial_idx,:);
    M_self{sub_i} = temp_M;
    M_other{sub_i} = temp_M;
    

    % FOR ELIMINATING HIGH VIFS AND OTHER ratings quetions 
    temp_vif_idx = [];
    temp_vif_idx = find(vifs{sub_i} > 2.5);
    
    other_idx = [];
    self_idx = [];
    other_idx = find(contains(OvrType(sub_i,sorted_idx),'Other'));
    self_idx = find(contains(OvrType(sub_i,sorted_idx),'Self'));
    
    final_other_idx = [];
    final_self_idx = [];
    
    %
    final_other_idx = unique([temp_vif_idx other_idx]);
    final_self_idx = unique([temp_vif_idx self_idx]);
    
    if ~isempty(final_self_idx) &  ~isempty(final_other_idx)
        % for self ratings
        M_self{sub_i}(final_other_idx) = [];
        Cond_heat_self{sub_i}(final_other_idx) = [];
        Cond_cue_self{sub_i}(final_other_idx) = [];
        FinalIntesityRating_self{sub_i}(final_other_idx) = [];
        
        
        % for other ratings 
        M_other{sub_i}(final_self_idx) = [];
        Cond_heat_other{sub_i}(final_self_idx) = [];
        Cond_cue_other{sub_i}(final_self_idx) = [];
        FinalIntesityRating_other{sub_i}(final_self_idx) = [];
        
    end
    
end

M_self(17) = []; M_other(17) = [];
Cond_cue_self(17) = []; Cond_cue_other(17) = [];
Cond_heat_self(17) = []; Cond_heat_other(17) = [];
FinalIntesityRating_self(17) = []; FinalIntesityRating_other(17) = [];

FinalIntesityRating_self = cellfun(@transpose, FinalIntesityRating_self,'UniformOutput',false);
FinalIntesityRating_other = cellfun(@transpose, FinalIntesityRating_other,'UniformOutput',false);
%% FOR SELF ratings
clear med_vars;
med_vars.M = M_self;
med_vars.X1 = Cond_heat_self;
med_vars.X2 = Cond_cue_self;
med_vars.Y = FinalIntesityRating_self;
med_vars.model_name = modeldir;
med_vars.scripts = '/sas1/cocoanlab/data/SEMIC_for_HPC/functions_matlab_toolbox/Mediation_dream_test/SEMIC_191023_making_variables_for_2mmWholeBrain_only_self.m';
med_vars.Descrip = {'single-trial whole-brain mediation'};
%% save
savenames = fullfile('/sas1/cocoanlab/data','SEMIC','data','meta_data','HPC_SEMIC_191203_mediation_dream_2mmWholeBrain_only_self.mat');
save(savenames,'med_vars');

%% FOR OTHER ratings
clear med_vars;
med_vars.M = M_other;
med_vars.X1 = Cond_heat_other;
med_vars.X2 = Cond_cue_other;
med_vars.Y = FinalIntesityRating_other;
med_vars.model_name = modeldir;
med_vars.scripts = '/sas1/cocoanlab/data/SEMIC_for_HPC/functions_matlab_toolbox/Mediation_dream_test/SEMIC_191023_making_variables_for_2mmWholeBrain_only_self.m';
med_vars.Descrip = {'single-trial whole-brain mediation'};
%% save
savenames = fullfile('/sas1/cocoanlab/data','SEMIC','data','meta_data','HPC_SEMIC_191203_mediation_dream_2mmWholeBrain_only_other.mat');
save(savenames,'med_vars');


