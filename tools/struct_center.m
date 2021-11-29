function [comb_struct] = struct_center(struct_array, center, even, dim, varargin)
%% struct_center: centers properties in structure arround one value of a field
%
%   Variable inputs can be used to determine which field is the
%   normalizaiton field, and which are centered. Default in 1st field for
%   normalization & all other fields to be centered.
%
%   INPUT:
%       struct_array 	:   input structure or object array (1st index is
%                          	the normalization field)
%       center          :   center value for normalization parameter
%       even            :   (boolean) for even padding around center
%       dim             :   dimension to operate on
%       min_sz          :   minimium size for nan padding
%   
%   OUTPUT:
%       comb_struct  	:   combined structure array
%
%   Usage:
%       comb_struct = struct_center(struct_array, 0, true, 1, [], {'field_1',field_2',... ,'field_n'})
%

% Get field names in stucture
fnames = string(fieldnames(struct_array));

% Check for user inputs
if nargin == 5 % only specified fields
    norm_fname = string(varargin{1}{1});
    cent_fname = string(varargin{1}(2:end))';
elseif nargin < 5 % minsz not specified, use default
    norm_fname = fnames(1);
    cent_fname = fnames(2:end);
    if nargin < 4
        dim = 1; % default
        if nargin < 3
            even = false; % default
            if nargin < 2
                center = 0; % default
            end
        end
    end
else
    error('Too many input arguments')
end
n_cent = length(cent_fname);

if isempty(center)
    center = 0; % default
end

if isempty(even)
    even = false; % default
end

if isempty(dim)
    dim = 1; % default
end

% Get normalization field
idx_norm = strcmp(norm_fname,fnames);
if ~any(idx_norm)
    error('%s does not exist', norm_fname)
end

% Get centering fields
idx_cent = false(length(fnames),n_cent);
for kk = 1:n_cent
    idx_cent(:,kk) = logical(strcmp(cent_fname(kk),fnames));
    if ~any(idx_cent(:,kk))
        error('%s does not exist', cent_fname(kk))
    end
end

% Normalize & center normalization field
field_norm = {struct_array.(norm_fname)}';

% Check for all nan's
n_cell = length(field_norm);
nan_cond = false(n_cell,1);
for jj = 1:n_cell
    nan_cond(jj) = all(all(isnan(field_norm{jj}))) && isscalar(field_norm{jj});
end
nan_cond = all(nan_cond);

if ~nan_cond
    if ~isempty(field_norm{1})
        [comb_struct.(norm_fname),~,~,~,dR] = nancat_center(field_norm, center, dim, [], even);
        
        % Center other fields
        for kk = 1:n_cent
            field_cent = {struct_array.(cent_fname(kk))}';
            temp = cellfun(@(x,y) padmat(x,y,nan,dim), field_cent, dR, 'UniformOutput', false);
            comb_struct.(cent_fname(kk)) = cat(2,temp{:});
        end
    else
        for n = 1:length(fnames)
            comb_struct.(fnames(n)) = nan;
        end
        % comb_struct = [];
    end
else % return nan's if nan's or empty in norm field
    comb_struct.(norm_fname) = nan;
    for kk = 1:n_cent
        comb_struct.(cent_fname(kk)) = nan;
    end  
end

end