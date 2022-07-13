%% INFO:
%
% Script for converting Tom's ape data for use at the Gottingen Conference
% Toolbox demonstration.
%
% This script is only for demo purposes and has not been validated.
%
%--------------------------------------------------------------------------


close all; clear all; clc;

% Loop through data:
file_data_array = dir('.\data\');
file_data_array(~endsWith({file_data_array.name}, '.tsv')) = [];

for file_data = file_data_array(:)'
    
    fn = [file_data.folder filesep file_data.name];
    
    raw_tpl_data_cell = readcell(fn, 'FileType', 'text');
    
    % Convert to table:
    data_header_row = 1;
    raw_tpl_data_full = cell2table(raw_tpl_data_cell((data_header_row+1):end, :) ...
        ,  'VariableNames' , raw_tpl_data_cell(data_header_row,:));
    
    % Make timestamp:
    raw_tpl_data_full.TimeStamp = raw_tpl_data_full.("Computer timestamp");
    
    % Split per session:
    [sessions, IA, IC] = unique(raw_tpl_data_full.("Timeline name"), 'stable');
    assert(issorted(IC), 'Non contiguous sessions.');
    
    for session_num = 1:numel(sessions)
        
        % Get session data:
        raw_tpl_data = raw_tpl_data_full( ...
            strcmp(raw_tpl_data_full.("Timeline name"), sessions{session_num}), :);
        assert(strcmp(raw_tpl_data.Event(1), 'RecordingStart') ...
            , 'Error slicing table into sessions.');
        
        % Make session name:
        [~, session_name, ~] = fileparts(fn);
        session_name = [session_name ' (' sessions{session_num} ')']; %#ok<AGROW>
        
        % Fix non-monotonically increasing time, and zero the time:
        s_incr = makeStrictlyIncrease(raw_tpl_data.TimeStamp);
        raw_tpl_data = raw_tpl_data(s_incr.new_indx, :);
        zeroTime_ms = raw_tpl_data.TimeStamp(1);
        raw_tpl_data.TimeStamp   = raw_tpl_data.TimeStamp - zeroTime_ms;
        
        % Assemble the diameter data.
        diam_data = struct('t_ms', raw_tpl_data.TimeStamp ...
            , 'L', str2double(raw_tpl_data.("Pupil diameter left")) ...
            , 'R', str2double(raw_tpl_data.("Pupil diameter right")) ...
            );
        noData = isnan(diam_data.L) & isnan(diam_data.R);
        diam_data.t_ms(noData) = [];
        diam_data.L(noData)    = [];
        diam_data.R(noData)    = [];
        diam_data.L            = diam_data.L ./ 1000;
        diam_data.R            = diam_data.R ./ 1000;
        
        % Make events:
        listIn = raw_tpl_data.("Presented Stimulus name");
        [labels, hitSections] = analyzeList(raw_tpl_data.TimeStamp ...
            , listIn);
        
        % Make AOIs:
        labelsCell = raw_tpl_data.("FlangeSize");
        labelsCell(cellfun(@(c) isa(c, 'missing'), labelsCell)) = {''};
        BB = cellfun(@num2str, labelsCell, 'UniformOutput', false);
        blanks = {'_' ' ' '.'};
        labelsCell(ismember(labelsCell,blanks)) = {''};
        [~, eventSections] = analyzeList(raw_tpl_data.TimeStamp, labelsCell);
        gazeData = struct(...
            'snapshotName', 'EXP'...
            ,'snapshotWidth', 800 ...
            ,'snapshotHeight', 600);
        gazeData.RoI   = eventSections;
        
        % Make data:
        pdtData                           = struct();
        pdtData.data.eyeTracking.diameter = diam_data;
        pdtData.data.eyeTracking.labels   = labels;
        pdtData.data.eyeTracking.name     = session_name;
        pdtData.data.eyeTracking.gazeData = gazeData;
        pdtData.data.eyeTracking.eventSections = [];
        pdtData.data.eyeTracking.raw_t_ms_max = raw_tpl_data.TimeStamp(end);
        pdtData.data.eyeTracking.diameterUnit = 'mm';
        
        
        % Before saving, add some metadata to the file (optional):
        pdtData.physioDataInfo.rawDataSource       = fn;
        pdtData.physioDataInfo.pdtFileCreationDate = datestr(now);
        pdtData.physioDataInfo.pdtFileCreationUser = 'Elio';
        
        % Save contents of the data struct (pdtData) to a physioData file:
        save(['.\data\' session_name '.physioData'] ...
            , '-struct', 'pdtData');
        
        
    end
    
end

