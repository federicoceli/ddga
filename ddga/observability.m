function O_y = observability(sys,T)
% OBSERVABILITY	Generates the observability matrix of sys = (A,B,C).
%   O_y = OBSERVABILITY(sys,T)
%        _          _
%       |   C        |
%       |   CA       |
% O_y = |   CA^2     |
%       |   ....     |
%       |   CA^(T-1) |
%        ¯          ¯

% F. Celi and F. Pasqualetti 2022

A = sys.A;
B = sys.B;
C = sys.C;

O_y = zeros(size(C,1)*T,size(C,2));

for r = 1:T
	O_y(1+(r-1)*size(C,1):size(C,1)+(r-1)*size(C,1),:) = C*A^(r-1);
end