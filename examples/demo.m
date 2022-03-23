%%
close all
clear all
addpath('../ddga/')
tol = 1e-9;
set_graphics;

%% System setup

A = [1 -.5 -.5 -.5;
	 0  .5   1  -2;
	 0   0   0   1;
	 0   0 -.5 1.5];
B = [1 0  0 0;
	 0 0 -1 0]';
C = [1 -1 0 -3];

n = 4;
m = 2;
p = 3;

sys = ss(A,B,C,[]);

if rank(ctrb(A,B)) == n
	disp('The system is controllable.');
else
	disp('The system is not controllable.');
end

if rank(obsv(A,C)) == n
	disp('The system is observable.');
else
	disp('The system is not observable.');
end

T = n;
N = n + m*T;

%% Model-based geometric control invariants
[Vstar_mb, F_mb] = vstar(A,B,C);
Sstar_mb = sstar(A,B,C);
Rstar_mb = ints(Vstar_mb,Sstar_mb);

%% Data-driven geometric control invariants
% Generate data as multiple (N) experiments of length T
data = generateData(sys,T,[],'multiple');

% Theorem 3.1 and Theorem 3.6 (note that F is not unique)
[Vstar_dd, F_dd] = vstardd(data,tol);
if (subspace(Vstar_dd, Vstar_mb) < tol && size(Vstar_mb,2) == size(Vstar_dd,2))
	disp('Succesfully computed Vstar!');
end

% Theorem 3.4
Sstar_dd = sstardd(data,tol);
if (subspace(Sstar_dd, Sstar_mb) < tol && size(Sstar_mb,2) == size(Sstar_dd,2))
	disp('Succesfully computed Sstar!');
end

% Remark 2
Rstar_dd = ints(Vstar_dd,Sstar_dd);
if (subspace(Rstar_dd, Rstar_mb) < tol && size(Rstar_mb,2) == size(Rstar_dd,2))
	disp('Succesfully computed Rstar!');
end