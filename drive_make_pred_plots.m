clear;
clc;

D_start_list = [50,50,50,50,50, 100,100,100, 150, 190, 220, 260];
D_end_list = [200,240,280,320,360, 250,300,350, 350, 290, 340, 380];
N_list = [100,150,200,240,280, 150,200,250, 250, 200, 250, 280];

% D_start_list = [50,50,50,100,  150,200,100,150,50,50];
% D_end_list = [250,280,300,300,  350,400,350,400,360,400];
% N_list = [120,200,160,160,  240,320,240,320,280,320];

indx_list = 1:length(N_list);

T_pred = 15;

error_v = zeros(1,length(indx_list));

c = [1,1,1];
c_bar = sqrt(linspace(0,1,20)')*[1,1,1];

for ii = 1:length(indx_list)
    iter_num = indx_list(ii);
    D_start = D_start_list(iter_num);
    D_end = D_end_list(iter_num);
    N_size = N_list(iter_num);
    
    titlename = strcat('start=',num2str(D_start),',end=',num2str(D_end),',N=',num2str(N_size));
    filename = strcat(num2str(D_start),'_',num2str(D_end),'_',num2str(N_size));
    fileloc = 'figures\';
    
    addpath('methods');
    addpath('covid_data');
    addpath('exp_val_data');
    load(strcat(filename,'_recon9_all_pred.mat'));
    load('original_data.mat');
    
    alpha = 0.05;
    N = 300;
    load(strcat('exp_simu_data\',filename,'_simu_pred.mat')); 
    
    D_end = D_end + T_pred;
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
    
    %% interpolation
    t = sum(Cum_cases);
    x1 = linspace(0,(length(t)-1),20000000);
    y_pchip = interp1(0:(length(t)-1),t,x1,'pchip');
    [~,indx] = unique(floor(y_pchip));
    T_list = x1(indx(2:end));
    
    rearrange_cases = zeros(1, t(end));
    rearrange_index = zeros(1, t(end));
    k = 1;
    for i = 1:size(Cum_cases,1)
        y_pchip = interp1(0:(length(t)-1),Cum_cases(i,:),x1,'pchip');
        [~,indx] = unique(floor(y_pchip));
        x_new = x1(indx(2:end));
        int_indx = (abs(x_new - round(x_new))<1e-4);
        x_new(int_indx) = x_new(int_indx) + 1e-9*rand(1,length(x_new(int_indx)));
        rearrange_cases(k:k+length(indx)-2) = x_new;
        rearrange_index(k:k+length(indx)-2) = i;
        k = k+length(indx)-1;
    end
    
    [~,sort_ind] = sort(rearrange_cases);
    types_new = rearrange_index(sort_ind);
    
    H = [types_new', T_list'];
    
    D_end = D_end_list(iter_num);
    p_val = 0.05;
    nn = floor(N*p_val);
    v_sort = sort(v);
    v_cut = v_sort(nn+1:end-nn,:);
    
    average_v = mean(v_cut);
    
    h0 = histogram(H(:,2),0:1:max(H(:,2)));
    v0 = h0.Values;
    
    V0 = v0(ceil(max(H_train(:,2))):ceil(max(H_train(:,2)))+T_pred-1);
    v1 = v_cut(:,ceil(max(H_train(:,2))):ceil(max(H_train(:,2)))+T_pred-1);
    
    F = figure;
    v_error_neg = abs(v_cut(1,:)-average_v);
    v_error_pos = abs(v_cut(end,:)-average_v);
    
    average_v = average_v*N_size;
    v_error_neg = v_error_neg*N_size;
    v_error_pos = v_error_pos*N_size;
    v0 = v0*N_size;
    
    V1 = average_v(ceil(max(H_train(:,2))):ceil(max(H_train(:,2)))+T_pred-1);
    
    p1 = errorbar(D_end_list(iter_num)+1:D_end_list(iter_num)+T_pred,...
        V1,v_error_neg(ceil(max(H_train(:,2))):ceil(max(H_train(:,2)))+T_pred-1),...
        v_error_pos(ceil(max(H_train(:,2))):ceil(max(H_train(:,2)))+T_pred-1),...
        'linewidth',1,'color',0.5*c);
    hold on;
    p2 = plot(D_start:(D_end+T_pred),v0,'linewidth',2,'color',0*c);
    hold on;
    p3 = plot(D_end_list(iter_num)+1:D_end_list(iter_num)+T_pred,...
        V1,'--','linewidth',2,'color',0*c);
    hold on;
    plot((D_end+1)*[1,1],[0,max(v0)+max(v_error_pos)],':','color',0.3*c,'linewidth',3);
    
%     legend([p1,p2,p3],{'Simulation','Ground Truth','Simulation Average'},'location','northwest','FontSize', 18);
    title(titlename,'FontSize', 28);
    ax = gca;
    ax.YAxis.Exponent = 3;
    set(gca,'fontsize',18);
    img = frame2im(getframe(F));
%     imwrite(img, strcat(fileloc,filename,'_pred.bmp'));
    print(F,strcat(fileloc,filename,'_pred_short.eps'),'-depsc')
    
    error_v(ii) = sum((V0 >= v1(1,:)).*(V0 <= v1(end,:)))/length(V0);
end

close all;

error_v