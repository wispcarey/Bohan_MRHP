clear;
clc;

Omega1 = {1,1,1,1.5,1,1,1.5,1.5};
Alpha1 = {2,2,2,1.5,2,2,1.5,1.5};
K = {'s','r','s','s','r','r','r','s'};
Mu = {'s','r','s','s','r','r','s','r'};
n_list = {10,10,6,6,6,6,6,6};

addpath('methods');

for i = 1:length(K)
    i
    T = 3000;
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
    y=simu_tMRHP(T,mu1,A1,omega1,h1);
    y.t = y.t;
    H = [y.type, y.t];
    
    [mu,A,omega,Lam,h,lkh,~,kk] = tempestim_recur(H,1000,1e-5);
    
    filename = strcat(num2str(omega1),'_',num2str(alpha1),'_',K{i},Mu{i},'_',num2str(n));
    
    save(strcat('exp_syn_data\',filename,'_ori.mat'),'H','mu','mu1','A','A1'...
        ,'omega','Lam','lkh','h','kk');
    
    H_floor = [y.type, floor(y.t)];
    
    [mu,A,omega,Lam,h,lkh,~,kk] = tempestim_recur(H_floor,1000,1e-5);
    
    save(strcat('exp_syn_data\',filename,'_floor.mat'),'H','mu','mu1','A','A1'...
        ,'omega','Lam','lkh','h','kk');
    
    for j_indx = 1
        [sort_yt, indx_yt] = sort(floor(y.t)+rand(size(y.t)));
        H_rand_recon = [y.type(indx_yt), sort_yt];

        [mu,A,omega,Lam,h,lkh,~,kk] = tempestim_recur(H_rand_recon,1000,1e-5);
        
        alpha = -log(h(1));
        A0 = A/alpha;
        
        if j_indx == 1
            mu_avg = mu;
            A_avg = A0;
            omega_avg = omega;
            alpha_avg = alpha;
            kk_avg = kk;
        else
            mu_avg = mu + (1-1/j_indx)*(mu_avg-mu);
            A_avg = A0 + (1-1/j_indx)*(A_avg-A0);
            omega_avg = omega + (1-1/j_indx)*(omega_avg-omega);
            alpha_avg = alpha + (1-1/j_indx)*(alpha_avg-alpha);
            kk_avg = kk + (1-1/j_indx)*(kk_avg-kk);
        end
    end
            
    
    save(strcat('exp_syn_data\',filename,'_rand.mat'),'H','mu','mu1','A','A1'...
        ,'omega','Lam','lkh','h','kk','mu_avg','A_avg','omega_avg','alpha_avg','kk_avg');
    
    for j_indx = 1
        H_recon = smooth_t_pchip_MHP_rand(H_floor);
    
        [mu,A,omega,Lam,h,lkh,~,kk] = tempestim_recur(H_recon,1000,1e-5);
        alpha = -log(h(1));
        A0 = A/alpha;
        
        if j_indx == 1
            mu_avg = mu;
            A_avg = A0;
            omega_avg = omega;
            alpha_avg = alpha;
            kk_avg = kk;
        else
            mu_avg = mu + (1-1/j_indx)*(mu_avg-mu);
            A_avg = A0 + (1-1/j_indx)*(A_avg-A0);
            omega_avg = omega + (1-1/j_indx)*(omega_avg-omega);
            alpha_avg = alpha + (1-1/j_indx)*(alpha_avg-alpha);
            kk_avg = kk + (1-1/j_indx)*(kk_avg-kk);
        end
    end
    
    save(strcat('exp_syn_data\',filename,'_recon.mat'),'H','mu','mu1','A','A1'...
        ,'omega','Lam','lkh','h','kk','mu_avg','A_avg','omega_avg','alpha_avg','kk_avg');
end
        
