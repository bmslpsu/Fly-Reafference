function [TF,Den,Num] = tf2sym(func, var, delay_flag)
%% tf2sym: convert transfer function object to symbolic transfer function
%   INPUTS:
%       tf    	:   transfer function object
%       var    	:   transform variable (character)
%   OUTPUTS:
%       TF      :   symbolic transfer function
%       Num    	:   transfer function numerator
%       Den    	:   transfer function denominator

a = func.Numerator;
b = func.Denominator;
dim = size(a); % # of tf's in matrix

if nargin < 3
    delay_flag = false;
end

if nargin == 1
    var = "s";
elseif nargin == 2
    if isa(var,'sym') || ischar(var)
        var = string(var);
    elseif isstring(var)
        % pass
    else
        error('Error: variable must be a character, string, or symbolic variable')
    end
end

TF  = sym(1);
Num = cell(dim);
Den = Num;

eval(sprintf('syms %s',var))

for jj = 1:dim(1)
    for kk = 1:dim(2)
        num = '';
        pp = length(a{jj,kk});
        for ii = 1:length(a{jj,kk})
            num = sprintf([num ' + a{jj,kk}(%i)*%s^{%i}'],ii,var,pp-1);
            pp = pp - 1;
        end

        den = '';
        pp = length(b{jj,kk});
        for ii = 1:length(b{jj,kk})
            den = sprintf([den ' + b{jj,kk}(%i)*%s^{%i}'],ii,var,pp-1);
            pp = pp - 1;
        end

        Num{jj,kk} = eval(num);
        Den{jj,kk} = eval(den);
        TF(jj,kk)  = eval(num)/eval(den);
        
        if delay_flag
            TF(jj,kk) = TF(jj,kk)*exp(-s*func.IODelay);
        end
    end
end

end