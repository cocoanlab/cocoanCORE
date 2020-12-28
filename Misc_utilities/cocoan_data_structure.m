function out = cocoan_data_structure(project_name, basedir, varargin)

% This function will create data structure in our workstation computer. 
%
% :Usage:
% ::
%
%    out = cocoan_data_structure(project_name, basedir, varargin)
%
% :Inputs:
%
%   **project_name:**
%        project name (e.g., 'Pleasure')
%
%   **basedir:**
%        base directory (e.g., '/Volumes/sein/dropbox')
%
% :Optional Inputs: Enter keyword followed by variable with values
%
%   **'project_only':**
%        With this option, you can make the directories in the project folder only.
%
%   **'data_only':**
%        With this option, you can make the directories in the data folder only.
%
%   **'create':**
%        With this option, this function will actually create the directories.         
%        Thus, without this optional input, it will just show the names of the
%        directories it will create (i.e., 'dry run' is the default)
%
%   **'dry_run' (default):**
%        This will just show the names of the directories it will create 
%
% :Examples:
% ::
%    % inputs
%    project_name = 'Pleasure';
%    basedir = ''/Volumes/sein/dropbox'';
%    out = cocoan_data_structure(project_name, basedir); % dry_run!!!
%    
%    % After you check all the directories are correct
%    out = cocoan_data_structure(project_name, basedir, 'create');
%
%    % project only
%    out = cocoan_data_structure(project_name, basedir, 'project_only');
%
%    % data only
%    out = cocoan_data_structure(project_name, basedir, 'data_only');
%
% ..
% Copyright (C) 04/07/2018, Choong-Wan Woo, Cocoan Lab
% ..

% Programmers' notes:
%
% updated to version 2.0 (June 14 2019) by Suhwan Gim 

dry_run = true;
make_data = true;
make_project = true;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'create'}
                dry_run = false;
            case {'dry_run'}
                dry_run = true;
            case {'project_only'}
                make_data = false;
            case {'data_only'}
                make_project = false;
        end
    end
end

line = '-----------------------------------------------------------------';

disp(line);
disp('   Directory checking...');
disp(line);

if make_data && ~exist(fullfile(basedir, 'data'), 'dir')
    error('\nThe ''data'' directory %s does not exist!!! Please check!!!', fullfile(basedir, 'data'));
end

if make_project && ~exist(fullfile(basedir, 'projects'), 'dir')
    error('\nThe ''projects'' directory %s does not exist!!! Please check!!!', fullfile(basedir, 'projects'));
end

if make_data
    % ---------------------------------------------------------------------------------
    % Data: store raw to preproceessed data that can be used in projects 
    %  *Keep in mind that your data can be used for other projects
    % ---------------------------------------------------------------------------------
    data_dir{1} = fullfile(basedir, 'data', project_name);
    data_dir{2} = fullfile(basedir, 'data', project_name, 'imaging');
    data_dir{3} = fullfile(basedir, 'data', project_name, 'imaging', 'preproc_script');
    data_dir{4} = fullfile(basedir, 'data', project_name, 'imaging', 'dicom_from_scanner');    
    data_dir{5} = fullfile(basedir, 'data', project_name, 'behavioral');
    data_dir{6} = fullfile(basedir, 'data', project_name, 'behavioral', 'raw');
    data_dir{7} = fullfile(basedir, 'data', project_name, 'behavioral', 'preprocessed');
    data_dir{8} = fullfile(basedir, 'data', project_name, 'physio', 'raw');
    data_dir{9} = fullfile(basedir, 'data', project_name, 'physio', 'preprocessed');
    data_dir{10} = fullfile(basedir, 'data', project_name, 'audiovideo');
    data_dir{11} = fullfile(basedir, 'data', project_name, 'metadata');
    data_dir{12} = fullfile(basedir, 'data', project_name, 'documents');
end

if make_project
    % ---------------------------------------------------------------------------------
    % Projects: store models and analysis that are specifially related to your projects 
    % ---------------------------------------------------------------------------------
    proj_dir{1} = fullfile(basedir, 'projects', project_name);
    proj_dir{2} = fullfile(basedir, 'projects', project_name,'sync');
    proj_dir{3} = fullfile(basedir, 'projects', project_name,'sync','writing');
    proj_dir{4} = fullfile(basedir, 'projects', project_name,'sync','literature');
    proj_dir{5} = fullfile(basedir, 'projects', project_name,'sync','data');
    proj_dir{6} = fullfile(basedir, 'projects', project_name,'sync','scripts');        
    proj_dir{7} = fullfile(basedir, 'projects', project_name,'sync','results');
    proj_dir{8} = fullfile(basedir, 'projects', project_name,'sync','figures');
    % These 'sync' folders will be also loacted in personal labtop and
    % shared with Wani through "dropbox"
    % --------------------------------------------------------------------------------
    proj_dir{9} = fullfile(basedir, 'projects', project_name, 'data');    
    proj_dir{10} = fullfile(basedir, 'projects', project_name, 'notes');
    proj_dir{11} = fullfile(basedir, 'projects', project_name, 'analysis');
    proj_dir{12} = fullfile(basedir, 'projects', project_name, 'analysis', 'behavioral');
    proj_dir{13} = fullfile(basedir, 'projects', project_name, 'analysis', 'imaging');
    proj_dir{14} = fullfile(basedir, 'projects', project_name, 'analysis', 'imaging','first_level');
    proj_dir{15} = fullfile(basedir, 'projects', project_name, 'analysis', 'imaging','first_level','model01');
    proj_dir{16} = fullfile(basedir, 'projects', project_name, 'analysis', 'imaging','second_level');
    proj_dir{17} = fullfile(basedir, 'projects', project_name, 'analysis', 'imaging','second_level','model01');
    proj_dir{18} = fullfile(basedir, 'projects', project_name, 'analysis', 'imaging','mediation');    
    proj_dir{19} = fullfile(basedir, 'projects', project_name, 'analysis', 'imaging','pattern');
    proj_dir{20} = fullfile(basedir, 'projects', project_name, 'analysis', 'imaging','connectivity');
    % You can modify these folder based on your project
end

if make_data
    disp(line);
    if dry_run
        disp('   Data directories it will create');
    else
        disp('   Data directories it is creating...');
    end
    disp(line);
    
    for i = 1:numel(data_dir)
        if dry_run
            str = 'WILL_CREATE: ';
        else
            str = 'CREATING: ';
        end
        disp([str data_dir{i}]);
        if ~dry_run
            mkdir(data_dir{i});
        end
    end
    % OUTPUT
    out.data_dir = data_dir;    
end

disp('');
disp('');
    
if make_project
    disp(line);
    if dry_run
        disp('   Project directories it will create');
    else
        disp('   Project directories it is creating...');
    end
    disp(line);
    
    for i = 1:numel(proj_dir)
        if dry_run
            str = 'WILL_CREATE: ';
        else
            str = 'CREATING: ';
        end
        disp([str proj_dir{i}]);
        if ~dry_run
            mkdir(proj_dir{i});
        end
    end
    % OUTPUT
    out.proj_dir = proj_dir;    
end

disp(line);
disp('Done');




end