clear;
clc;

addpath('methods');

% productivity matrix A1
n = 3;
a1 = tril(ones(n));
A1 = [a1,zeros(n);zeros(n),a1'];
A1 = A1/sum(sum(A1))*n;

% background rate mu1
mu1 = 0.1*ones(1,n);

% time interval [0,T]
T = 1000;

% spatio region
X = 5;
Y = 5;

% parameter in time kernel: omega1
omega1 = 1;

% parameter in spatio kernel: sig1
sig1 = 1;

% parameter in recursive function: alpha1
alpha1 = 2;
h1 = @(x)alpha1*exp(-alpha1*x);
% h1 = @(x)ones(size(x));

% simulation spatiotemporal point process
y=simu_stMRHP(X,Y,T,mu1,A1,omega1,sig1,h1);
H = [y.type, y.t, y.lon, y.lat];

[A,B,omega,Lam,h,sig,tau,p,pb,mu,kk] = stestim_recur(H);
% Nonparametric estimation
%res = nphawkes(H,X,Y);

h00 = h;

alpha = -log(h(1));
A = A/alpha;
h = @(x)alpha*h(x);

M = max(max(Lam));
m = min(min(Lam));

figure(1);
imagesc(A1);
colorbar;
ylabel('Parents','FontSize', 24);
xlabel('Offspring','FontSize', 24);
% title('K matrix: ground truth');
set(gca,'XAxisLocation','top');

figure(2);
imagesc(A);
colorbar;
ylabel('Parents','FontSize', 24);
xlabel('Offspring','FontSize', 24);
% title(strcat('K matrix: ',num2str(length(H(:,2)))));
set(gca,'XAxisLocation','top');
error_K = norm(A-A1)/norm(A1);

% figure(3);
% error_mat = abs(A-A1);
% imagesc(error_mat);
% colorbar;
% ylabel('Parents','FontSize', 24);
% xlabel('Offspring','FontSize', 24);
% title(strcat('Error of K: ',num2str(size(H(:,2)))));
% set(gca,'XAxisLocation','top');

figure(4);
plot(linspace(0,M,1000),h1(linspace(0,M,1000)),'-','linewidth',2,'color',0*c);
hold on;
plot(linspace(0,M,1000),h(linspace(0,M,1000)),'--','linewidth',2,'color',0*c);
legend('ground truth','fitted results','Fontsize',24);
title('Recursive function h','Fontsize',32);
set(gca,'fontsize',18);
error_alpha = abs(h1(0)-h(0))/abs(h1(0));


figure(5);
X = linspace(0,20,1000);
plot(X,omega1*exp(-omega1*X),'-','linewidth',2,'color',0*c);
hold on
plot(X,w*exp(-w*X),'--','linewidth',2,'color',0*c);
legend('ground truth','fitted results','Fontsize',24);
title('Time kernel g','Fontsize',32);
set(gca,'fontsize',18);
error_g = abs(w-omega1)/abs(omega1);

figure(6);
plot(1:size(A1,1), mu1,'.','markersize',20, 'color',0*c);
hold on;
plot(1:size(A1,1), u,'*','markersize',20, 'color',0*c);
axis([0,size(A1,1),0,0.20]);
title('Background rate \mu','Fontsize',32);
legend('ground truth','fitted results','Fontsize',24);
set(gca,'fontsize',18);
error_mu = norm(u-mu1)/norm(mu1);

error_k = zeros(size(Lam));
for i = 1:length(Lam)
    error_k(i) = norm(A1*h1(Lam(i))-A*h(Lam(i)))/norm(A1*h1(Lam(i)));
end
error_prod = mean(error_k);

format long
error_g
error_alpha
error_K
error_mu
error_prod
Num_sampled = length(H(:,2))
