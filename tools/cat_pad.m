function [nanData] = cat_pad(data,L,val)
%% nancat_append: pad a vector with specified value on both sides
%   INPUTS:
%       data     	:   input data (vector)
%       L           :   pad length on each side
%       val         :   pad value
%   OUTPUTS:
%       nanData     :   Nan padded data
%

data = data(:);
L = L(:);
nanData = [val*ones(L(1),1) ; data ; val*ones(L(2),1)];
end