clear;
clc;

% D_start_list = [50,100,150,200,50,100,150, 50,50,50,50,50,50];
% D_end_list = [250,300,350,400,300,350,400, 200,240,280,320,360,400];
% N_list = [120,160,240,320,160,240,320, 100,130,200,240,280,320];

D_start_list = [50,50,50,100,  150,200,100,150,50,50];
D_end_list = [250,280,300,300,  350,400,350,400,360,400];
N_list = [120,200,160,160,  240,320,240,320,280,320];

indx = 1:10;

A_rec = zeros(9,9,length(indx));

error_A = zeros(1,length(indx));
error_v = zeros(1,length(indx));

c = [1,1,1];
c_bar = sqrt(linspace(0,1,20)')*[1,1,1];

for ii = 1:length(indx)
    iter_num = indx(ii);
    D_start = D_start_list(iter_num);
    D_end = D_end_list(iter_num);
    N_size = N_list(iter_num);
    
    titlename = strcat('start=',num2str(D_start),',end=',num2str(D_end),',N=',num2str(N_size));
    filename = strcat(num2str(D_start),'_',num2str(D_end),'_',num2str(N_size));
    fileloc = 'figures\';
    
    addpath('methods');
    addpath('covid_data');
    addpath('exp_val_data');
    load(strcat(filename, '_recon9_all.mat'));
    
    alpha = 0.05;
    N = 300;
    load(strcat('exp_simu_data\',filename,'_simulation.mat')); 
    
    p_val = 0.05;
    nn = floor(N*p_val);
    v_sort = sort(v);
    v_cut = v_sort(nn+1:end-nn,:);
    
    average_v = mean(v_cut);
    
    figure;
    h0 = histogram(H(:,2),0:1:ceil(max(H(:,2))));
    v0 = h0.Values;
        
    F = figure;
    v_error_neg = abs(v_cut(1,:)-average_v);
    v_error_pos = abs(v_cut(end,:)-average_v);
    
    average_v = average_v*N_size;
    v_error_neg = v_error_neg*N_size;
    v_error_pos = v_error_pos*N_size;
    v0 = v0*N_size;
    
    p1 = errorbar(D_start:D_end,average_v,v_error_neg,v_error_pos,'linewidth',1,'color',0.5*c);
    hold on;
    p2 = plot(D_start:D_end,v0,'linewidth',2,'color',0*c);
    hold on;
    p3 = plot(D_start:D_end,average_v,'--','linewidth',2,'color',0*c);
    
%     legend('Simulation','Ground Truth','Simulation Average','location','northwest','FontSize', 14);
    title(titlename,'FontSize', 28);
    ax = gca;
    ax.YAxis.Exponent = 3;
    set(gca,'fontsize',18);
%     img = frame2im(getframe(F));
%     saveas(gcf, strcat(fileloc,filename,'_',num2str(alpha),'recon'), '-depsc');
    print(F,strcat(fileloc,filename,'_',num2str(alpha),'recon.eps'),'-depsc')
    
    K = figure;
    imagesc(A);
    colormap(c_bar);
    colorbar;
    ylabel('Parents','FontSize', 24);
    xlabel('Offspring','FontSize', 24);
    set(gca,'XAxisLocation','top');
%     imwrite(gcf, strcat(fileloc,filename,'_',num2str(alpha),'K_mat'), '-depsc');
    print(K, strcat(fileloc,filename,'_',num2str(alpha),'K_mat.eps'), '-depsc');
    A = A/norm(A,2);
    
    A_rec(:,:,ii) = A;
    
    if ii == 1
        mean_A = A;
    else
        mean_A = mean_A + (A-mean_A)/ii;
    end
end

for ii = 1:length(indx)
    iter_num = indx(ii);    
    D_start = D_start_list(iter_num);
    D_end = D_end_list(iter_num);
    N_size = N_list(iter_num);
    
    titlename = strcat('start=',num2str(D_start),',end=',num2str(D_end),',N=',num2str(N_size));
    filename = strcat(num2str(D_start),'_',num2str(D_end),'_',num2str(N_size));
    fileloc = 'new_figures_8_19\';
    
    addpath('covid_data');
    addpath('exp_val_data');
    load(strcat(filename, '_recon9_all.mat'));
    load(strcat('exp_8_14\',filename,'_simulation.mat')); 
    
    p_val = 0.05;
    nn = floor(N*p_val);
    v_sort = sort(v);
    v_cut = v_sort(nn+1:end-nn,:);
    
    average_v = mean(v_cut);
    
    figure;
    h0 = histogram(H(:,2),0:1:ceil(max(H(:,2))));
    v0 = h0.Values;
    
    A = A/norm(A,2);
    
    error_A(ii) = norm(A-mean_A,2)/norm(mean_A,2);
    
    v_ind = (v0 <= v_cut(end,:)).*(v0 >= v_cut(1,:));
    error_v(ii) = sum(v_ind)/length(v_ind);
end

close all;

error_A
error_v