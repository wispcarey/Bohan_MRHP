clear;
clc;

% D_start_list = [50,100,150,200,50,100,150, 50,50,50,50,50,50];
% D_end_list = [250,300,350,400,300,350,400, 200,240,280,320,360,400];
% N_list = [120,160,240,320,160,240,320, 100,130,200,240,280,320];

D_start_list = [50,50,50,100,  150,200,100,150,50,50];
D_end_list = [250,280,300,300,  350,400,350,400,360,400];
N_list = [120,200,160,160,  240,320,240,320,280,320];

indx = 1:length(N_list);

for ii = 1:length(indx)
    iter_num = indx(ii);
    D_start = D_start_list(iter_num);
    D_end = D_end_list(iter_num);
    N_size = N_list(iter_num);
    
    titlename = strcat('start=',num2str(D_start),',end=',num2str(D_end),',N=',num2str(N_size));
    filename = strcat(num2str(D_start),'_',num2str(D_end),'_',num2str(N_size));
    
    addpath('methods');
    addpath('covid_data');
    addpath('exp_val_data');
    load(strcat(filename, '_recon9_all.mat'));
    H = H_recon;
    
    alpha = 0.05;
    Hist = H(H(:,2)<= (max(H(:,2))*alpha),:);
    
    N = 300;
    bin_size = 1;
    v = simulation_all_hist(ceil(max(H(:,2))),u,A,w,h,Hist,N,bin_size);
    
    p_val = 0.05;
    nn = floor(N*p_val);
    v_sort = sort(v);
    v_cut = v_sort(nn+1:end-nn,:);
    
    average_v = mean(v_cut);
    
    h0 = histogram(H(:,2),0:1:ceil(max(H(:,2))));
    v0 = h0.Values;
    
    save(strcat('exp_simu_data\',filename,'_simulation.mat'), 'H', 'Hist', 'v'); 
        
end