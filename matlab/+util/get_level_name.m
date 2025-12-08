function levelName = get_level_name(D, context)
% GET_LEVEL_NAME  Extract unique level name from a design table.
%
% levelName = util.get_level_name(D, context)
%
% D       : table with a 'level' column (string/char)
% context : optional string for clearer error messages

if nargin < 2
    context = '';
end

if ~istable(D) || ~ismember('level', D.Properties.VariableNames)
    error('get_level_name:missingLevel', ...
        'Table is missing ''level'' column%s.', message_suffix(context));
end

u = unique(string(D.level));
if numel(u) ~= 1
    error('get_level_name:ambiguousLevel', ...
        'Ambiguous level labels (%s) in design table%s.', ...
        strjoin(cellstr(u), ', '), message_suffix(context));
end

levelName = char(u);
end

function s = message_suffix(context)
if ~isempty(context)
    s = sprintf(' (context: %s)', context);
else
    s = '';
end
end
