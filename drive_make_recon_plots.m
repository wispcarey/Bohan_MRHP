clear;
clc;

Omega1 = {1,1,1,1.5,1,1,1.5,1.5};
Alpha1 = {2,2,2,1.5,2,2,1.5,1.5};
K = {'s','r','s','s','r','r','r','s'};
Mu = {'s','r','s','s','r','r','s','r'};
n_list = {10,10,6,6,6,6,6,6};

c = [1,1,1];
c_bar = (linspace(0,1,20)').^(0.7)*[1,1,1];

% suffix can be chosen as 'ori','floor','rand','recon'
suffix = 'rand';

for i = 1
    i
    n = n_list{i};
    if strcmp(K{i},'s')
        if n_list{i} == 10
            load data_10_K
        else
            load data_6_K
        end
    else
        A1 = rand(n);
        A1 = A1/sum(sum(A1))*n/2;
    end
    if strcmp(Mu{i},'s')
        mu1= 0.1*ones(size(A1,1),1);
    else
        mu1 = rand(n,1);
        mu1 = mu1/sum(mu1)*0.1*n;
    end
    omega1 = Omega1{i};
    alpha1 = Alpha1{i};
    h1 = @(x)alpha1*exp(-alpha1*x);
    
    filename = strcat(num2str(omega1),'_',num2str(alpha1),'_',K{i},Mu{i},'_',num2str(n));
    
    load(strcat('exp_syn_data\',filename,'_',suffix,'.mat'));
    
    alpha = -log(h(1));
    A = A/alpha;
    h = @(x)alpha*h(x);
    
    M = max(max(Lam));
    m = min(min(Lam));
    %modifier = 1-exp(-0.5*M);
    
%     figure(1);
%     imagesc(A1);
%     colorbar;
%     xlabel('Parents');
%     ylabel('Offspring');
%     title('K matrix: ground truth');
%     set(gca,'XAxisLocation','top');
    
    figure(2);
    imagesc(A);
    colormap(c_bar);
    colorbar;
    ylabel('Parents','FontSize', 24);
    xlabel('Offspring','FontSize', 24);
%     title(strcat('K matrix: ',num2str(length(H(:,2)))));
    set(gca,'XAxisLocation','top');

%     figure(3);
%     error_mat = abs(A-A1);
%     imagesc(error_mat);
%     colorbar;
%     xlabel('Parents');
%     ylabel('Offspring');
%     title(strcat('Error of K: ',num2str(size(H(:,2)))));
    
    figure(4);
    plot(linspace(0,M,1000),h1(linspace(0,M,1000)),'-','linewidth',2,'color',0*c);
    hold on;
    plot(linspace(0,M,1000),h(linspace(0,M,1000)),'--','linewidth',2,'color',0*c);
    legend('ground truth','fitted results','Fontsize',24);
    title('Recursive function h','Fontsize',32);
    set(gca,'fontsize',18);
    
    figure(5);
    X = linspace(0,20,1000);
    plot(X,omega1*exp(-omega1*X),'-','linewidth',2,'color',0*c);
    hold on
    plot(X,omega*exp(-omega*X),'--','linewidth',2,'color',0*c);
    legend('ground truth','fitted results','Fontsize',24);
    title('Time kernel g','Fontsize',32);
    set(gca,'fontsize',18);
    
    figure(6);
    plot(1:size(A1,1), mu1,'.','markersize',20, 'color',0*c);
    hold on;
    plot(1:size(A1,1), mu,'*','markersize',20, 'color',0*c);
    axis([0,size(A1,1),0,0.20]);
    title('Background rate \mu','Fontsize',32);
    legend('ground truth','fitted results','Fontsize',24);
    set(gca,'fontsize',18);
end