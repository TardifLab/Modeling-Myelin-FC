function outi = cellstrfind(targ,str)
% Searches cell arrays for a given string occurence. Returns indices of
% cells in which target string found.
%
% Input:
%   targ    : Cell array to search for string
%   str     : target string
%
% Output:
%   outi    : indices of string
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
opt1='ignorecase';
opt2='UniformOutput';

if ischar(str)                                                              % search input is string vector
    tmp                 = cellfun(@(c) regexp(c,str,opt1),targ,opt2,0);
    outi                = find(cellfun(@(c) ~isempty(c),tmp));
elseif iscell(str)                                                          % search input may be a cell array of strings
    outi                = cell(size(str));
    for ii = 1:numel(str)
        if ischar(str{ii})
            tmp         = cellfun(@(c) regexp(c,str{ii},opt1),targ,opt2,0);
            outi{ii}    = find(cellfun(@(c) ~isempty(c),tmp));
        end
    end
else
    error('Unknown input data type: CELLSTRFIND accepts strings or cell arrays of strings')
end
        
%% @MCN BUILD IN AMBIGUITY HANDLING e.g. (or see conn_selectData.m)
% % % Handle potential ambiguity
% % if any(cellfun(@(c) numel(c),itrg)>1)                                       % if multiple indices found
% %     for kk = 1 : numel(itrg)
% %         keyboard
% %         itrg{kk}    = find(cellfun(@(c) strcmpi(Wtrg{kk},c),W8s));          % more rigid-search criteria
% %     end
% % end
% % if any(cellfun(@(c) numel(c),itrg)>1)                                       % if still ambiguous
% %     error('Unable to resolve ambiguity in your input to WTRG')
% % else
% %     itrg            = cell2mat(itrg);
% % end


%--------------------------------------------------------------------------
end