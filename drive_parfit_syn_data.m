clear;
clc;

type = {'t','t','t','t','t','t','st','st','st','st','st','st'};
Omega1 = {1,1.5,1.5,1.5,1.5,1, 1,1,1.5,1.5,1.5,1.5};
Alpha1 = {2,1.5,1.5,1.5,1.5,2, 2,1,1.5,1.5,1.5,1.5};
Sig1 = {0,0,0,0,0,0, 2,1,1.5,1.5,1.5,1.5};
K = {'s','s','r','s','r','r','s','s','s','s','r','r'};
Mu = {'s','s','s','r','r','r','s','s','s','r','s','r'};

addpath('methods');

for i = 1:length(type)
    i
    T = 3000;
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
        y=simu_tMRHP(T,mu1,A1,omega1,h1);
        H = [y.type, y.t];

        [mu,A,omega,Lam,h,lkh,~,kk] = tempestim_recur(H);
        
        filename = strcat(num2str(omega1),'_',num2str(alpha1),'_',num2str(sig1),'_',K{i},Mu{i},'_',type{i});
        
        save(strcat('exp_syn_data\',filename,'.mat'),'H','mu','mu1','A','A1'...
            ,'omega','Lam','lkh','h','kk');
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
        X = 5;
        Y = 5;
        y=simu_stMRHP(X,Y,T,mu1,A1,omega1,sig1,h1);
        H = [y.type, y.t, y.lon, y.lat];
        [A,B,omega,Lam,h,sig,tau,~,~,mu,kk] = stestim_recur(H);
        
        filename = strcat(num2str(omega1),'_',num2str(alpha1),'_',num2str(sig1),'_',K{i},Mu{i},'_',type{i});
        
        save(strcat('exp_syn_data\',filename,'.mat'),'H','mu','mu1','A','A1'...
            ,'B','omega','Lam','h','sig','kk');
    end
end
        
