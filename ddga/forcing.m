function F_y = forcing(sys,T)
% FORCING Generates the output forcing matrix of sys = (A,B,C).
%   F_y = FORCING(A,B,C,D,T)
%        _                          _
%       |   C          ...   0     0 |
%       |   CB         ...   0     0 |
% F_y = |   CAB        ...   0     0 |
%       |   ....       ...   ...   0 | 
%       |   CA^(T-2)B  ...   CB    0 |
%        ¯                          ¯

% F. Celi and F. Pasqualetti 2022

A = sys.A;
B = sys.B;
C = sys.C;

F_y = zeros(size(C,1)*T,size(B,2)*T);

for r=2:T
	for c=1:r-1
    	F_y((r-1)*size(C,1)+1:r*size(C,1),(c-1)*size(B,2)+1:c*size(B,2))=C*(A^(r-1-c))*B;
	end
end

