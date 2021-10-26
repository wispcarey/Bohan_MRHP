clear;
clc;

addpath('methods');
addpath('exp_syn_data');

type = {'t','t','t','t','t','t','st','st','st','st','st','st'};
Omega1 = {1,1.5,1.5,1.5,1.5,1, 1,1,1.5,1.5,1.5,1.5};
Alpha1 = {2,1.5,1.5,1.5,1.5,2, 2,1,1.5,1.5,1.5,1.5};
Sig1 = {0,0,0,0,0,0, 2,1,1.5,1.5,1.5,1.5};
K = {'s','s','r','s','r','r','s','s','s','s','r','r'};
Mu = {'s','s','s','r','r','r','s','s','s','r','s','r'};

error_g = zeros(size(type));
error_alpha = zeros(size(type));
error_sig = zeros(size(type));
error_K = zeros(size(type));
error_mu = zeros(size(type));
error_prod = zeros(size(type));
num_sample = zeros(size(type));
num_iter = zeros(size(type));
max_eigen = zeros(size(type));

for i = 1:length(type)
    if strcmp(type{i},'t')
        if strcmp(K{i},'s')
            load data_10_K
        else
            A1 = rand(10);
            A1 = A1/sum(sum(A1))*5;
        end
        if strcmp(Mu{i},'s')
            mu1= 0.1*ones(size(A1,1),1);
        else
            mu1 = rand(10,1);
            mu1 = mu1/sum(mu1);
        end
        omega1 = Omega1{i};
        sig1 = Sig1{i};
        alpha1 = Alpha1{i};
        h1 = @(x)alpha1*exp(-alpha1*x);
        
        filename = strcat(num2str(omega1),'_',num2str(alpha1),'_',num2str(sig1),'_',K{i},Mu{i},'_',type{i});
        
        load(strcat('exp_syn_data\',filename,'.mat'));
        alpha = -log(h(1));
        A = A/alpha;
        h = @(x)alpha*h(x);
        
        error_K(i) = norm(A-A1)/norm(A1);
        error_alpha(i) = abs(h1(0)-h(0))/abs(h1(0));
        error_g(i) = abs(omega-omega1)/abs(omega1);
        error_mu(i) = norm(mu'-mu1)/norm(mu1);
        max_eigen(i) = max(abs(eig(A1)));
        error_k = zeros(size(Lam));
        for j = 1:length(Lam)
            error_k(j) = norm(A1*h1(Lam(j))-A*h(Lam(j)))/norm(A1*h1(Lam(j)));
        end
        error_prod(i) = mean(error_k);
        num_sample(i) = size(H,1);
        num_iter(i) = kk;
    else
        if strcmp(K{i},'s')
            load data_6_K
        else
            A1 = rand(6);
            A1 = A1/sum(sum(A1))*3;
        end
        if strcmp(Mu{i},'s')
            mu1= 0.1*ones(size(A1,1),1);
        else
            mu1 = rand(6,1);
            mu1 = mu1/sum(mu1)*0.6;
        end
        omega1 = Omega1{i};
        sig1 = Sig1{i};
        alpha1 = Alpha1{i};
        h1 = @(x)alpha1*exp(-alpha1*x);
        
        filename = strcat(num2str(omega1),'_',num2str(alpha1),'_',num2str(sig1),'_',K{i},Mu{i},'_',type{i});
        
        load(strcat('exp_syn_data\',filename,'.mat'));
        alpha = -log(h(1));
        A = A/alpha;
        h = @(x)alpha*h(x);
        
        error_K(i) = norm(A-A1)/norm(A1);
        error_alpha(i) = abs(h1(0)-h(0))/abs(h1(0));
        error_g(i) = abs(omega-omega1)/abs(omega1);
        error_mu(i) = norm(mu'-mu1)/norm(mu1);
        error_sig(i) = abs(sig-sig1)/abs(sig1);
        max_eigen(i) = max(abs(eig(A1)));
        error_k = zeros(size(Lam));
        for j = 1:length(Lam)
            error_k(j) = norm(A1*h1(Lam(j))-A*h(Lam(j)))/norm(A1*h1(Lam(j)));
        end
        error_prod(i) = mean(error_k);
        num_sample(i) = size(H,1);
        num_iter(i) = kk;
    end
end