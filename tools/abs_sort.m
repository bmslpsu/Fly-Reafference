function [sortData] = abs_sort(data)
%% abs_sort: sort a vector by seperating the sorted positive & negative components and concatenating the results
%   INPUTS:
%       data        :   data to sort
%   OUTPUTS:
%       sortData	:   sorted data
%---------------------------------------------------------------------------------------------------------------------------------
if isscalar(data)
    dim = 1;
	sortDim = dim;
else
    dim = size(data);
    dIdx = dim==1;
    if sum(dIdx)>2 ||sum(dIdx)==0
        error('Error: input must be vector')
    end
    sortDim = find(dIdx==0);
end

pIdx    = data>=0;
nIdx    = data<0;
nanIdx  = isnan(data);

if any(nanIdx)
    nanSort = nan(1,dim(2));
else
	nanSort = [];
end

pSort = sort(data(pIdx),'ascend');
nSort = sort(data(nIdx),'descend');
sortData = cat(sortDim,pSort,nSort,nanSort);
end