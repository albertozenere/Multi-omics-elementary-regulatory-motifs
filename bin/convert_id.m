%To convert 'list' from inID to outID identifier, it is required to
%specify both inID and outID because some identifiers are shared by UniProtKBSwissProtID
%and UniProtKBGeneNameID.
%Returns a cell array containing the output identifiers corresponding each input
%identifier (there can be several)

function out = convert_id(list,inID,outID)
martexport = importdata('martexport020718.mat');

%%
key = {'GenestableID', 'ensg';
       'ProteinstableID', 'ensp';
       'Genename', 'name';
       'UniProtKBSwissProtID', 'accp';
       'UniProtKBGeneNameID', 'accg';
       'NCBIgeneID', 'entrez'}; %abbreviations, to be updated if you download new categories from Ensembl/biomart

%% check format
assert(ismember(inID, key(:,2)) & ismember(outID, key(:,2)),'the input & output identifiers must be on the left column of key above');   

if isa(list,'char')
    N = 1;
else
    assert(min(size(list))==1, 'input must be one column'); %check it is a vector
    N = max(size(list));
    
    if size(unique(list),1)~=N
        list_full = list;
        [list, ~,list_ib] = unique(list_full);
    end
end

%%
idx_in = ismember(key(:,2), inID);
idx_out = ismember(key(:,2), outID); 

if isa(list, 'string')
    list = arrayfun(@(x) char(list(x)),1:numel(list),'uni',false)';
end
switch class(list)
    case 'cell'
        out = cell(N,1);
        subset_mart = martexport(ismember(martexport{:,idx_in},list) & ~ismember(martexport{:,idx_out},''),:);
        if ~isempty(subset_mart)
            [u, ia, ib] = unique(subset_mart{:,idx_in});
            
            [idx, ord] = ismember(list, u);
            
            group_u = arrayfun(@(n) unique(subset_mart{ib==n,idx_out}), 1:length(ia), 'UniformOutput', false)';
            out(idx) = group_u(ord(idx))';
        end
    case 'char'
        out = unique(martexport{ismember(martexport{:,idx_in},list), idx_out});
end

if exist('list_full','var')
    out_full = out(list_ib);
    out = out_full;
end
