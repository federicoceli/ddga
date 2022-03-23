%%
close all
clear all
addpath('../ddga/')
tol = 1e-9;
set_graphics;

n = 11;
I = eye(n);
T = 19;
h = 0.2;
tol = 1e-4;

A = [-1 1 0 0 0 0 0 0 0 0 0;
	1 -3 1 0 0 0 0 0 0 0 0;
	0 1 -2 0 0 0 0 0 0 0 0;
	0 0 0 -2 1 0 0 0 0 0 0; 
	0 0 0 1 -3 1 0 0 0 0 0;
	0 0 0 0 1 -2 0 1 0 0 0;
	0 0 0 0 0 0 -1 1 0 0 0;
	0 0 0 0 0 1 1 -4 0 0 1;
	0 0 0 0 0 0 0 0 -2 1 0;
	0 0 0 0 0 0 0 0 1 -2 0;
	0 0 0 0 0 0 0 1 0 0 -1];

C = I([4 11],:);

B = [0 0 0; 1 0 0; 0 0 1; 1 0 0; 1 0 0; 0 0 0; 0 0 0;
	0 1 0; 0 1 0; 0 1 0; 0 0 0];

% inputs
m = size(B,2);
% outputs
p = size(C,1);

% Forward Euler
A = I + h*A;
B = h*B;
C = C;

sys = ss(A,B,C,[],h);

%% Model-based geometric control invariants
[Vstar, F] = vstar(A,B,C);
Sstar = sstar(A,B,C);
Rstar = ints(Vstar,Sstar);

%% Data-driven geometric control invariants
% Generate data as multiple (N) experiments of length T
data = generateData(sys,T,[],'multiple');

[Vstar_dd, F_dd] = vstardd(data);
Sstar_dd = sstardd(data,tol);

for i = 1 : size(B,2)
	B_bar(:,i) = zeros(n,1);
	if subspace(B(:,i), Vstar) < 1e-7
		B_bar(:,i) = B(:,i);
	end
end

x1 = randn(n,1); 
x2 = x1;

attack_time = ceil(T/2);

est_ima_B = ima(data.X(n+1:2*n,:)*null(data.X0));

B_bar = ints(Vstar_dd,est_ima_B);

for t = 1 : T
	u(:,t) = [-2;2;4];
	a(:,t) = zeros(size(B_bar,2),1);
	if t > attack_time
		a(:,t) = 10*rand(size(B_bar,2),1);
	end

	x1(:,t+1) = A*x1(:,t) + B*(u(:,t)) + B_bar*a(:,t);
	x2(:,t+1) = A*x2(:,t) + B*(u(:,t));

	y1(:,t) = C*x1(:,t);
	y2(:,t) = C*x2(:,t);
end

%%
figure;
subplot(3,1,1:2)
hold on
plot(x1')
line([attack_time, attack_time], [min(min(x1)), max(max(x1))], 'Color','red','LineStyle','--')
title(['Undetectable attack (attack from time t = ',num2str(attack_time),'s)'])
ylim([min(min(x1)),max(max(x1))])
hold off
ylabel('state (x)')
subplot(3,1,3)
hold on
plot(y1')
line([attack_time, attack_time], [min(min(y1)), max(max(y1))], 'Color','red','LineStyle','--')
hold off
%title(['Output (attack from time t = ',num2str(attack_time),'s)'])
ylim([min(min(y1)),max(max(y1))])
xlabel('time [s]')
ylabel('output (y)')
set(gcf,'position',[400 250 650 250])
% exportgraphics(gcf,'attack.pdf','ContentType','vector')

%%
%{
figure;
subplot(2,1,1)
plot(x2')
title('State (no attack)')
ylabel('x')
subplot(2,1,2)
plot(y2')
title('Output (no attack)')
xlabel('time [s]')
ylabel('y')

%%
figure;
subplot(2,1,1)
plot(abs(x1'-x2'))
title(['Distance from equilibrium (attack from time t = ',num2str(attack_time),'s)'])
ylabel('$||x-x_{eq}||$','interpreter','latex')
subplot(2,1,2)
plot(abs(y1'-y2'))
%title(['Output distance from eq (attack from time t = ',num2str(attack_time),'s)'])
xlabel('time [s]')
ylabel('$||y-y_{eq}||$','interpreter','latex')
% exportgraphics(gcf,'attack_distance.pdf','ContentType','vector')
%}