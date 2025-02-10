function out = runsort(filelist)
% 
% :: Usage
%       out = runsort(filelist)
%
% :: Input
%       - filelist: A cell array or string array of filenames
%                   containing 'run-XX' patterns.
%
% :: Output
%       - out: The input filelist sorted based on the numerical 
%              value of 'run-XX' in each filename.
%
% :: Description
%       This function extracts 'run-XX' patterns from the input filenames, 
%       sorts them numerically (not lexicographically), and returns 
%       the sorted list.
%
% :: Examples 
%
% filelist = ["experiment_run-02_data.txt", ...
%             "experiment_run-01_data.txt", ...
%             "experiment_run-10_data.txt", ...
%             "experiment_run-03_data.txt"];
%
% sorted_files = runsort(filelist);
% disp(sorted_files)
%
% Expected output:
%     "experiment_run-01_data.txt"
%     "experiment_run-02_data.txt"
%     "experiment_run-03_data.txt"
%     "experiment_run-10_data.txt"
%
% See also, REGEXP, SORT, CELLFUN
% 
% Author: Eui-Jin Jung (your.email@example.com)
% Date: YYYY-MM-DD
%

% Extract 'run-XX' patterns from filenames
runname = regexp(filelist, 'run-[0-9][0-9]', 'match');

% Flatten the nested cell array (if using cell array input)
for i = 1:numel(runname)
    runname{i} = runname{i}{1};
end

% % Extract numerical values from 'run-XX' strings
% run_numbers = cellfun(@(x) str2double(extractAfter(x, "run-")), runname);

% Sort numerically
[~, runsortidx] = sort(runname);

% Reorder the filelist based on the sorted indices
out = filelist(runsortidx);

end