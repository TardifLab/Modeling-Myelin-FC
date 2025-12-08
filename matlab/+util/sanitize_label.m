function s = sanitize_label(label)
% SANITIZE_LABEL  Make a label safe for folder/filenames.
%
% - Converts to char
% - Trims spaces
% - Replaces internal spaces with underscores
% - Replaces file separators with dashes (just in case)

s = char(label);
s = strtrim(s);
s = strrep(s, ' ', '_');
s = strrep(s, filesep, '-');
end
