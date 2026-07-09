function sub_id = rdcm_extract_subid(filename, pattern)
%% Preamble
%{
Extracts a subject identifier from a filename using a regex pattern.

Examples:
  rdcm_extract_subid('sub-CMH001_task-rest_bold.tsv', 'sub-CMH[0-9]+')
    -> 'sub-CMH001'
  rdcm_extract_subid('sub-01_ses-01_bold.tsv', 'sub-[0-9]+')
    -> 'sub-01'

Set params.sub_id_pattern to control matching, e.g.:
  params.sub_id_pattern = 'sub-CMH[0-9]+';   % SUD dataset
  params.sub_id_pattern = 'sub-[A-Za-z0-9]+'; % generic BIDS (default)

If no match is found, sub_id falls back to the filename stem (everything
before the first underscore or dot) with a warning.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
filename:=  char — filename string (basename or full path)
pattern:=   char — MATLAB regexp pattern stored in params.sub_id_pattern
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
sub_id:=    char — matched subject identifier string
%}

[~, basename, ext] = fileparts(filename);
basename_full      = [basename ext];

tokens = regexp(basename_full, pattern, 'match');

if isempty(tokens)
    parts  = strsplit(basename_full, {'_', '.'});
    sub_id = parts{1};
    warning('rdcm_extract_subid:noMatch', ...
        'Pattern "%s" did not match filename "%s". Falling back to "%s".', ...
        pattern, basename_full, sub_id);
else
    sub_id = tokens{1};
end

end