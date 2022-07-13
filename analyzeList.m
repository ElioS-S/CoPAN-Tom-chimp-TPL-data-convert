function [labels, hitSections] = analyzeList(t_ms, listIn, fs ...
    , gap_n_dt, isDebug)
%analyzeList Analyzes contiguous sections in timestamped categorical data.
%
%    This function finds contiguous sections in a timstamped list, and
%    returns a struct array containing they're names, and their start and
%    end times.
%
%    If not passed as an input, the function estimates the fs by taking the
%    median of the diff of t_ms. Gaps, defined as delta t errors larger
%    than k the expected delta t given the (estimated) fs, are then
%    identified. Sections spanning such a gap are not considered
%    contiguous.
%
%    Each hit starts at the timestamp of the first value of the contigious
%    section of that hit, and ends at the timestamp of the row after the
%    last value in said section, unless there is a gap following the last
%    section, in which case the end timestamp is defined as the timestamp
%    of the last value plus the dt (calculated from the estimated or passed
%    fs).
%
%--------------------------------------------------------------------------
%
%   Elio Sjak-Shie, August 2019.
%    https://physiodatatoolbox.leidenuniv.nl/
%
%--------------------------------------------------------------------------



arguments
    
    % The time vector and corresponding cellstr array:
    t_ms   = [];
    listIn = {};
    
    % The fs (will be estimted from t_ms if not specified):
    fs     = [];

    % The gap threshold ()in multiples of the fs:
    gap_n_dt = [];
    
    % The debug flag:
    isDebug = false;
    
    
end

% Set default gap value:
if isempty(gap_n_dt)
    gap_n_dt = 2;
end
assert(gap_n_dt >= 0, 'Incorrect Gap N parameter.');

% Set default values, and return if nothing was passed in:
labels = struct(...
    'name',{{}}...
    ,'t_ms', [] ...
    );
hitSections = struct(...
    'name',{}...
    ,'hit',{}...
    ,'desc',{}...
    ,'fs',{}...
    );
if isempty(listIn) || isempty(t_ms)
    return
end

assert(size(listIn,2)==1,'listIn must be a 1 column array.');
assert(size(listIn,1)==size(t_ms,1)...
    ,'listIn and t_ms must have the same number of columns.');

% In no fs was passed, estimate it by taking the median diff of t_ms. Note
% that since there might be gaps in the data, simply taking the mean diff
% of t_ms won't work:
if isempty(fs)
    
    [fs, s] = estimateFs(t_ms);
    dt_ms = s.dt_ms;
    
else
    
    % Calculate the delta t form the fs:
    dt_ms = 1000 / fs;
    
    % If an fs was passed, run checks on the t_ms that estimateFs would
    % have otherwise run:
    assert(size(t_ms,2)==1,'t_ms must be a 1 column vector.');
    assert(all(diff(t_ms)>0)...
        ,'t_ms must be strictly monotonically increasing.');
    
end

% Calculate the error:
t_ms_2           = t_ms;
t_ms_2(2:end)    = t_ms(1:end-1) + dt_ms;
dt_extra_ms      = t_ms - t_ms_2;


% Replace missing with empty string:
listIn(cellfun(@(c) isa(c, 'missing'), listIn)) = {''};
listIn = cellfun(@num2str, listIn, 'UniformOutput', false);

% TODO (HIGH): Add test with missing and numbers.


% Transform data in into a categorical array (making it easier to analyze:
sectionCode = categorical(listIn);

% Find gaps in the data (so that hits are not merged over gaps). A gap is
% defined as an error of N_DT_GAP_THRESH times the dt. Use the
% insertNanAtBool function to insert <undefined> samples in the gaps:
isGapEnd  = dt_extra_ms > (gap_n_dt * dt_ms);
sectionCode_SplitAtGaps ...
    = insertNanAtBool(sectionCode,[isGapEnd(2:end);false]);
t_ms_SplitAtGaps...
    = insertNanAtBool(t_ms,[isGapEnd(2:end);false]);

clearvars sectionCode

% Find the first values of each contiguous section, defined as: A value
% that is: (1) different from the value preceding it or the first value in
% the array; (2) not <undefined>:
firstValueIndx  = [true;sectionCode_SplitAtGaps(2:end) ...
    ~= sectionCode_SplitAtGaps(1:end-1)];
firstValueIndx(ismissing(sectionCode_SplitAtGaps)) = false;
onset_t_ms      = t_ms_SplitAtGaps(firstValueIndx);

% Track the categories corresponding to the sections:
hitCat          = sectionCode_SplitAtGaps(firstValueIndx);

% Similarly, find the last values (values that are different than the next
% or are the last value; and are not undefined):
lastValueIndx = [sectionCode_SplitAtGaps(1:end-1) ...
    ~= sectionCode_SplitAtGaps(2:end); true];
lastValueIndx(ismissing(sectionCode_SplitAtGaps)) = false;
assert(size(firstValueIndx,1) == size(lastValueIndx,1)...
    ,'Onset/offset detection error.');

% Set the section-ending timestamps equal to the timestamps of the rows
% following the last samples (or the last):
nextRowTimeIndx   = min(find(lastValueIndx) + 1, numel(t_ms_SplitAtGaps));
offset_t_ms       = t_ms_SplitAtGaps(nextRowTimeIndx);

% However, if the next timestamp is a gap-end, the timestamps won't work
% because they have been split; i.e., they are NaNs. As such, detect these
% NaNs, and replace them with the local timestamps, plus dt:
offsetNaNs              = isnan(offset_t_ms);
offset_t_ms(offsetNaNs) = t_ms_SplitAtGaps(...
    nextRowTimeIndx(offsetNaNs) -  1) + dt_ms;

% Additionally, if the last section ends at the last sample, add dt to it.
% This corresponds to how it is handled throughout the list, and prevents a
% crash when he last section is one sample long (section onset would be the
% same as the section offset):
if nextRowTimeIndx(end) == numel(t_ms_SplitAtGaps)
    offset_t_ms(end) = offset_t_ms(end) + dt_ms;
end

% Do a final sanity check:
assert(isequal(sum(firstValueIndx),numel(offset_t_ms),numel(onset_t_ms))...
    ,'Error with vectors and such.');

% Make labels (note that the concatenation, transposition and linear
% indexing assures that the labels are sorted):
startEventNames = strcat({'Start_'},cellstr(hitCat));
endEventNames   = strcat({'End_'},cellstr(hitCat));
eventNames      = [startEventNames endEventNames]';
labels.name     = eventNames(:);
StartEndTimes   = [onset_t_ms offset_t_ms]';
labels.t_ms     = StartEndTimes(:);

% Generate the sections:
uniqueSegementNames = categories(sectionCode_SplitAtGaps);
for curEvtIndx = 1 : numel(uniqueSegementNames)
    curName  = uniqueSegementNames{curEvtIndx};
    curRows  = hitCat == curName;
    
    % Enter data into current index:
    hitSections(curEvtIndx).name  = curName;    
    hitSections(curEvtIndx).hit   = [onset_t_ms(curRows) ...
        offset_t_ms(curRows)];
    hitSections(curEvtIndx).desc = ['Section demarked by "' curName '".']; 
    
end

% Save the fs:
[hitSections.fs] = deal(fs);


%% DEBUG CODE:

if isDebug
    % Prints the findings to the console:
    
    maxLen = max(cellfun(@numel,listIn));
    fprintf('\n')
    firstValueIndxDISP = firstValueIndx(~isnan(t_ms_SplitAtGaps));
    lastValueIndxDISP = lastValueIndx(~isnan(t_ms_SplitAtGaps));
    elemDesc = cell(size(t_ms));
    elemDesc(firstValueIndxDISP) ...
        = arrayfun(@(s,e) {sprintf('DUR: %.1f (%.1f -> %.1f)',e-s,s,e)}...
        ,onset_t_ms,offset_t_ms);
    arrayfun(@(ii) fprintf(...
        '\n%7.1f: %s | %s %s %s | %s'...
        ,t_ms(ii)...
        ,pad(listIn{ii},maxLen + 1)...
        ,ifFcn(firstValueIndxDISP(ii),'<-START','       ')...
        ,ifFcn(lastValueIndxDISP(ii),'<-LAST','      ')...
        ,ifFcn(isGapEnd(ii),'<-GAP','     ')...
        ,ifFcn(firstValueIndxDISP(ii),elemDesc{ii},'')...
        ) ...
        ,1:numel(firstValueIndxDISP))
    fprintf('\n\n')
    
end

end



%% Function Tests:

function testMe() %#ok<DEFNU>


% Make string sections:
secCodeAndTimes = {...
    0  '';...  Nothing.
    1  'A';... The first section (A1) starts here.
    2  'A';... 
    3  'A';... Last value for A1 (eventhough the next valid value is A).
    4  '';...  This should become the ending timestamp for A1.
    5  '';...
    6  'A';... A2 starts here.
    7  'A';... This is the last A2 value.
    8  'B';... First B1 value, and end timestamp of A1.
    9  'B';...
    10 '';...  End time of B1.
    11 'C';... Start C1.
    12 'C';...
    13 'C';...
    14 'C';... Last value C1
    15 'A';... Start A3, last value A3, and end time C1.
    ...        Note that since this is a gap, A3 should end at 16, which is
    ...        the timestamp of its last value (15), plus the estimated fs,
    20 '';...  which is 1.
    21 '';...
    22 '';...
    23 'A';... Start and last value of A4.
    ...
    ...
    30 'A';... Start and last value of A5.
    ...
    ...
    40 '';... etc.
    41 'C';...
    42 '';...
    43 '';...
    44 'A';...
    45 'A';...
    46 'A';...    
    49 'A';...
    50 'A';...
    53 '';...
    };
secCode  = secCodeAndTimes(:,2);
sec_t_ms = [secCodeAndTimes{:,1}]';
[~, hitSections] = analyzeList(sec_t_ms,secCode,[],[],true);
disp(' eventSections:')
arrayfun(@(s) disp([{[s.name ':'] ''}; num2cell(s.hit)]),hitSections)

end





