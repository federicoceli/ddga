function [H,is_full_rank] = hankel(signal,order)
% HANKEL	Builds the Hankel matrix of a given order the a (n x T) signal.
%   is_full_rank returns if H is full row-rank. 

% F. Celi and F. Pasqualetti 2022

u = signal;
L = order;

T = size(u,2);
m = size(u,1);

H = zeros(m*L, T - L + 1);

for i = 1 : T - L + 1
    H(:,i) = reshape(u(:,i:i + L - 1),m*L,1);
end

is_full_rank = ( size(H,2) == rank(H) );