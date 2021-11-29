function [paddata] = padmat(matdata,L,val,dim)
%% padmat: pad a matrix with specified value on both sides
%
%   INPUT:
%       matdata   	:   input data (matrix)
%       L           :   pad length on each side
%       val         :   pad value
%       dim         :   dimension
%   OUTPUT:
%       paddata     :   padded data
%

sz = size(L);

if all(sz==1)
    L(2) = 0;
end

if dim==1
    paddata = [val*ones(L(1),sz(2)) ; matdata ; val*ones(L(2),sz(2))];
elseif dim==2
    paddata = [val*ones(sz(1),L(1)) , matdata , val*ones(sz(1),L(2))];
else
    error('Must be along 1st or 2nd dimension')
end

end