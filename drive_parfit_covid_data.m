clear;
clc;

%% data preprocess

% D_start_list = [50,50,50,50,50, 100,100,100, 150, 190, 220, 260];
% D_end_list = [200,240,280,320,360, 250,300,350, 350, 290, 340, 380];
% N_list = [100,150,200,240,280, 150,200,250, 250, 200, 250, 280];
D_start_list = [50,50,50,100,  150,200,100,150,50,50];
D_end_list = [250,280,300,300,  350,400,350,400,360,400];
N_list = [120,200,160,160,  240,320,240,320,280,320];

for iter_num = 1:length(N_list)
    clearvars -except D_start_list D_end_list N_list iter_num
    clf;
    
    addpath('methods');
    addpath('covid_data');
    load('original_data.mat');
    
    D_start = D_start_list(iter_num);
    D_end = D_end_list(iter_num);
    N_size = N_list(iter_num);
    
    titlename = strcat('start=',num2str(D_start),',end=',num2str(D_end),',N=',num2str(N_size));
    filename = strcat(num2str(D_start),'_',num2str(D_end),'_',num2str(N_size));
    fileloc = 'figures_new\exp_val\';
    
    Cum_cases = Cum_cases(1:9,:);
    
    Cum_cases_new = Cum_cases(:,(D_start+1):(D_end+2))-Cum_cases(:,D_start+1);
    
    %% smoothing data by 7 days
    daily_reported_cases = Cum_cases_new(:,2:end)-Cum_cases_new(:,1:end-1);
    daily_reported_cases = round(smoothdata(daily_reported_cases, 2, 'movmean', 7));
    
    %% smoothing data by 7 days
    for column_ind = 2:size(Cum_cases_new,2)
        Cum_cases_new(:,column_ind) = sum(daily_reported_cases(:,1:column_ind-1),2);
    end
    Cum_cases = round(Cum_cases_new/N_size);
    
    daily_cases = Cum_cases(:,2:end)-Cum_cases(:,1:end-1);
    
    %% random smoothing
    t = sum(Cum_cases);
    timeline = zeros(t(end),1);
    types = zeros(t(end),1);
    k = 1;
    
    for date_ind = 1:size(daily_cases,2)
        for type_ind = 1:size(daily_cases,1)
            if daily_cases(type_ind,date_ind) ~= 0
                timeline(k:k + daily_cases(type_ind,date_ind) - 1) = date_ind - 1;
                types(k:k + daily_cases(type_ind,date_ind) - 1) = type_ind;
                k = k + daily_cases(type_ind,date_ind);
            end
        end
    end
    
    H = [types, timeline];
    
    [sort_t, indx_t] = sort(timeline+rand(size(timeline)));
    H_rand = [types(indx_t), sort_t];
    
    H_recon = smooth_t_pchip_MHP_rand(H);
    
    %% fitting
%     [u,A,w,Lam,h,lkh,~,~,~,kk] = tempestim_recur(H_rand);
%     
%     save(strcat('exp_val_data\',filename,'_rand_all.mat'),'H_rand',...
%         'A','u','w','Lam','h','kk','lkh');
    
    
    %% fitting
    [u,A,w,Lam,h,lkh,~,kk] = tempestim_recur(H_recon);
    
    save(strcat('exp_val_data\',filename,'_recon9_all.mat'),'H_recon',...
        'A','u','w','Lam','h','kk','lkh');
    
end




