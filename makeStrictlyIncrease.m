function s_out = makeStrictlyIncrease(t_in)
% Fixes a timevector if it is not strickly monotonically increasing.
%
%  The time vector is first sorted, then repeating entries are deleted (the
%  first one is kept). NaNs are also removed.
%
%--------------------------------------------------------------------------


%% 

% Check:
assert(isvector(t_in), 't_in must be a vector.');
t_in = t_in(:);

% Check if it was already sorted:
was_sorted = issorted(t_in);

% Sort time and get the sort index:
[t_sort, t_sort_indx] = sort(t_in);

% Check for duplicates (don't count NaN dupes, or diffs to NaN, which are
% also NaN, and thus don't satisfy >0):
t_sort_diff     = diff(t_sort);
t_dupes         = [~((t_sort_diff > 0) | isnan(t_sort_diff)); false];
had_dupes       = any(t_dupes);

% Make a new time vector that is sorted, and w/o dupes:
t_new      = t_sort(~t_dupes);

% Log the index to transform a vector from the t_in space, to the t_new
% space:
new_indx   = t_sort_indx(~t_dupes);

% Check if t had NaNs:
new_NaNs  = isnan(t_new);
had_NaNs  = any(new_NaNs);

% Remove NaNs:
t_new     = t_new(~new_NaNs);
new_indx = new_indx(~new_NaNs);

% Assemble output:
s_out.was_sorted = was_sorted;
s_out.had_dupes  = had_dupes;
s_out.had_NaNs   = had_NaNs;
s_out.t_new      = t_new;
s_out.new_indx   = new_indx;


end