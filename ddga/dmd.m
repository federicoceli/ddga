function [A,B,C,D] = dmd(X,U,Y)
% DMD	Dynamic Mode Decomposition (with control)
%   [A,B,C,D] = DMD(X,U,Y)
%   X : state trajectory  (length T+1)
%   U : input trajectory  (length T)
%   Y : output trajectory (length T)

% From Dynamic Mode Decomposition with Control
% by Joshua L. Proctor, Steven L. Brunton, and J. Nathan Kutz
% SIAM J. APPLIED DYNAMICAL SYSTEMS 2016

% F. Celi and F. Pasqualetti 2022

X1 = X(1:end-1,:);
X2 = X(2:end,:);

AB = X2'*pinv([X1'; U']);
CD = Y' *pinv([X1'; U']);

n = size(X,2);

A = AB(:,1:n);
B = AB(:,n+1:end);

C = CD(:,1:n);
D = CD(:,n+1:end);