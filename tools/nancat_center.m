function [centData,cIdx,cVal,R,dR,new_length,n_vector] = nancat_center(data,center,dim,Cent,Even)
%% nancat_center: centers input data around given value and pads with Nan's so all vectors are the same size
%   INPUTS:
%       data       	:   input data (matrix or cell of vectors)
%       center     	:   data center value (set to empty when using Cent)
%       dim         :   dimension to center around
%       Cent     	:   center index , overrides center value
%       Even     	:   if true,equal sizing on each side of "center"
%   OUTPUTS:
%       centData    :   centered data
%       cIdx        :   center indicies
%       cVal        :   computed center values
%       R           :   indicies from center to ends of each column
%       maxR       	:   length of both sides of "center"
%       new_length 	:   new length of Nan padded vectors
%       n_vector 	:   total # of vectors in all cells
%

if nargin < 3
    dim = 1; % defaults to columns
    if nargin < 2
        center = 0; % default center value is 0
    end
elseif nargin == 4
    if ~isempty(center)
        error('Error: the input "center" must be empty when using the input "Cent"')
    end
    center = 0;
end

if iscell(data) % for cell array
elseif ismatrix(data) % if input is a matrix, convert matrix to single cell
    data = {data};
else
    Error('Error: input must be matrix or cell array') 
end

if dim~=1 && dim~=2 % only works on first two dimesions
    error('Error: dimension must be 1 or 2')
else
    data = data(1:end); % put all cells in column
    if dim==2 % transform to column vectors if row vectors are specified
        data = cellfun(@(x) transpose(x), data, 'UniformOutput', false);
    end
end
[n_data, ~] = size(data); % # of cells

[n_center, ~] = cellfun(@(x) size(x), data, 'UniformOutput', true); % size of dimension to center around in each cell

if n_center <= 2
    error('Error: need at least 3 values in array')
end

[cVal, ~] = cellfun(@(x) min(abs(x - center)), data, ... % get value closest to center value for each cell
    'UniformOutput', false);
% cVal = cell2mat(cVal);

R = cell(n_data,1); % store the # of indicies above and below the center value
cIdx = nan(size(data));
n_array = nan(n_data,1);
n_length = nan(n_data,1);
for jj = 1:n_data % for each cell
    [n_length(jj) , n_array(jj)] = size(data{jj});
    R{jj} = nan(2,n_array(jj));
    for kk = 1:n_array(jj)
        if ~isnan(cVal{jj}(kk))
            [cIdx(jj,kk),~] = find( data{jj}(:,kk)==cVal{jj}(1,kk) ); % find location of center value
        else
            cIdx(jj,kk) = 1;
        end
        
        if nargin == 4
            if isempty(Cent)
                % no effect
            elseif isnan(Cent)
                cIdx(jj,kk) = n_center(jj); % for last index
            else
                cIdx(jj,kk) = Cent;
            end
        end
        R{jj}(1,kk) = length(1:cIdx(jj,kk)-1); % length of data below center
        R{jj}(2,kk) = length((cIdx(jj,kk)+1):n_center(jj)); % length of data above center
    end
end

maxR = cellfun(@(x) max(x,[],2)', R, 'UniformOutput', false);
maxR = cell2mat(maxR);
maxR_ALL = max(maxR,[],1)';

if nargin == 5
    if Even % even padding around center
        maxR_ALL = max(maxR_ALL)*ones(size(maxR_ALL));
    end
end

n_vector = sum(n_array); % total # of vectors from cells
new_length = sum(maxR_ALL) + 1; % new length with Nan's

dR = cell(n_data,1); % offset
centData = cell(n_data,1); % centered data
for jj = 1:n_data % for each cell
    for kk = 1:n_array(jj)
        dR{jj}(:,kk) = abs(maxR_ALL - R{jj}(:,kk));
        centData{jj}(:,kk) = cat_pad(data{jj}(:,kk),dR{jj}(:,kk),nan);
    end
end
centData = cat(2,centData{:});

end