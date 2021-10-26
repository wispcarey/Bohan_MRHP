clear;
clc;

addpath('methods');
addpath('exp_syn_data');

type = 'par';

if strcmp(type,'par')
    load('data_par.mat');
    alpha = -log(h(1));
    h = @(x)alpha*h(x);
    A = A/alpha;
else
    load('data_non_par.mat');
end

M = max(max(Lam));
m = min(min(Lam));

c = [1,1,1];
c_bar = (linspace(0,1,20)').^(0.7)*[1,1,1];

figure(1);
imagesc(A1);
colormap(c_bar);
colorbar;
ylabel('Parents','FontSize', 24);
xlabel('Offspring','FontSize', 24);
%     title(strcat('K matrix: ',num2str(length(H(:,2)))));
set(gca,'XAxisLocation','top');


figure(2);
imagesc(A);
colormap(c_bar);
colorbar;
ylabel('Parents','FontSize', 24);
xlabel('Offspring','FontSize', 24);
%     title(strcat('K matrix: ',num2str(length(H(:,2)))));
set(gca,'XAxisLocation','top');

error_K = norm(A-A1)/norm(A1);

% figure(3);
% error_mat = abs(A-A1);
% imagesc(error_mat);
% colorbar;
% xlabel('Parents');
% ylabel('Offspring');
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
plot(1:size(A1,1), 10*mu1,'.','markersize',20, 'color',0*c);
hold on;
plot(1:size(A1,1), u,'*','markersize',20, 'color',0*c);
axis([0,size(A1,1),0,0.20]);
title('Background rate \mu','Fontsize',32);
legend('ground truth','fitted results','Fontsize',24);
set(gca,'fontsize',18);
error_mu = norm(u-mu1*10)/norm(mu1*10);

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