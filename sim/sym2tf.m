function [TF,Num,Den] = sym2tf(func,var,Ts)
%% sym2tf: convert symbolic transfer function to transfer function object
%   INPUTS:
%       func  	:   symbolic transfer function
%   OUTPUTS:
%       TF      :   transfer function object
%       Num    	:   transfer function numerator coefficents
%       Dem    	:   transfer function denominator coefficents
%

% Let user specify transfer variable
if nargin==1
    var = symvar(func);
elseif nargin==2
    var = sym(var);
elseif nargin==3
    if isempty(var)
        var = symvar(func);
    else
        var = sym(var);
    end
end
tfvar = char(var);

dim = size(func); % dimensions of symbolic matrix
TF  = tf; % transfer function matrix
Num = cell(dim); % numerator coefficents
Den = cell(dim); % denominator coefficents

var = unique(var);
if isempty(var) % purely numeric transfer function, no computations needed
    Num = double(func);
    Den = 1;
 	TF = tf(Num);
    return
elseif length(var)>1 % make sure there is only one variable
   error('TF may contain up to one symbolic variable')    
end

% For all symbolic expresions in the matrix
for jj = 1:dim(1)
    for kk = 1:dim(2)
        % Get numerator and denominator polynomial coefficents, store in cells
        if isnumeric(func(jj,kk)) % purely numeric transfer function, no computations needed
        	Num{jj,kk}  = func(jj,kk);
            Den{jj,kk}  = 1;
        else
            [num,den]   = numden(func(jj,kk));
            Num{jj,kk}  = double(sym2poly(num));
            Den{jj,kk}  = double(sym2poly(den));
        end
        
        % Construct transfer function object
        if nargin==3
            TF(jj,kk) = tf(Num{jj,kk},Den{jj,kk},Ts);
        else
            TF(jj,kk) = tf(Num{jj,kk},Den{jj,kk});
        end
        TF(jj,kk).Variable = tfvar;
    end
end

% If symbolic matrix is 1x1, return numerator and denominator coefficents as arrays (not in cells)
if length(TF)==1
    Num = Num{1};
    Den = Den{1};
end

end