function visualize_network_webweb(A, varargin)

% visualize_network(A, varargin)
%
% optional inputs
% case {'webwebdir'}
% case {'node_names'}
% case {'groups'}
% case {'group_names'}
% case {'values'}
% case {'value_names'}


webwebdir = '/Users/clinpsywoo/Dropbox/MATLAB/network_tools/webweb/matlab';

vis_groups = false;
vis_values = false;
vis_nodename = false;
use_group_names = false;
use_value_names = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'webwebdir'}
                webwebdir = varargin{i+1};
            case {'node_names'}
                vis_nodename = true;
                node_names = varargin{i+1};
            case {'groups'}
                vis_groups = true;
                if numel(varargin{i+1}) == 1
                    groups = varargin{i+1};
                else
                    for j = 1:numel(varargin{i+1})
                        groups{j} = varargin{i+1}{j};
                    end
                end
            case {'group_names'}
                use_group_names = true;
                if ~iscell(groups)
                    group_names = varargin{i+1};
                else
                    for j = 1:numel(groups)
                        group_names{j} = varargin{i+1}{j};
                    end
                end
            case {'values'}
                vis_values = true;
                if ~iscell(varargin{i+1})
                    values = varargin{i+1};
                else
                    for j = 1:numel(varargin{i+1})
                        values{j} = varargin{i+1}{j};
                    end
                end
            case {'value_names'}
                use_value_names = true;
                if ~iscell(values)
                    value_names = varargin{i+1};
                else
                    for j = 1:numel(values)
                        value_names{j} = varargin{i+1}{j};
                    end
                end
        end
    end
end


try
    cd(webwebdir);
catch
    disp('Please specify the directory that includes the Webweb tool in the code.');
    return
end

dis.w = 300;
dis.h = 300;

% Increase the charge and the gravity
dis.c = 40;
dis.g = 0.05;
dis.l = 13;
dis.r = 4;

% Give the file a name
dis.name = 'visualization_network';
% Name the nodes

if ~vis_nodename
    for i=1:size(A,1)
        dis.nodeNames{i} = num2str(i);
    end
else
    for i=1:size(A,1)
        dis.nodeNames{i} = node_names{i};
    end
end

nets.net.adj = A;

if ~vis_groups && ~vis_values
    webweb(dis, nets);
else
    if vis_groups
        if ~use_group_names
            if ~iscell(groups)
                nets.net.labels.group.type = 'categorical';
                nets.net.labels.group.values = groups;
            else
                for j = 1:numel(groups)
                    eval(['nets.net.labels.group' num2str(j) '.type = ''categorical'';']);
                    eval(['nets.net.labels.group' num2str(j) '.values = groups{' num2str(j) '};']);
                end
            end
        else
            if ~iscell(groups)
                eval(['nets.net.labels.' group_names '.type = ''categorical'';']);
                eval(['nets.net.labels.' group_names '.values = groups;']);
            else
                for j = 1:numel(groups)
                    eval(['nets.net.labels.' group_names{j} '.type = ''categorical'';']);
                    eval(['nets.net.labels.' group_names{j} '.values = groups{' num2str(j) '};']);
                end
            end
            
        end
    end
    if vis_values
        if ~use_value_names
            if ~iscell(values)
                nets.net.labels.group.type = 'scalar';
                nets.net.labels.group.values = values;
            else
                for j = 1:numel(values)
                    eval(['nets.net.labels.value' num2str(j) '.type = ''scalar'';']);
                    eval(['nets.net.labels.value' num2str(j) '.values = values{' num2str(j) '};']);
                end
            end
        else
            if ~iscell(values)
                eval(['nets.net.labels.' value_names '.type = ''scalar'';']);
                eval(['nets.net.labels.' value_names '.values = values;']);
            else
                for j = 1:numel(groups)
                    eval(['nets.net.labels.' value_names{j} '.type = ''scalar'';']);
                    eval(['nets.net.labels.group' value_names{j} '.values = values{' num2str(j) '};']);
                end
            end
            
        end
    end
    
    webweb(dis,nets);
end

end
