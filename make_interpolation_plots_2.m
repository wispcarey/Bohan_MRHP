clear;
clc;

addpath('methods');

t = [0,3,5,10,12];
t1 = [0,1,2,5,6];
t2 = t-t1;
x1 = linspace(0,(length(t)-1),10000);
%Shape-preserving piecewise cubic interpolation
y_pchip = interp1(0:(length(t)-1),t,x1,'pchip');
[~,indx] = unique(floor(y_pchip));
T_list = x1(indx(2:end))';
y_recon = interp1(T_list,1:length(T_list),x1,'previous');
y_recon(isnan(y_recon)) = 0;

y_pchip = interp1(0:(length(t)-1),t1,x1,'pchip');
[~,indx] = unique(floor(y_pchip));
T_list_1 = x1(indx(2:end))';
y_recon_1 = interp1(T_list,1:length(T_list),x1,'previous');
y_recon_1(isnan(y_recon_1)) = 0;

y_pchip = interp1(0:(length(t)-1),t2,x1,'pchip');
[~,indx] = unique(floor(y_pchip));
T_list_2 = x1(indx(2:end))';
y_recon_2 = interp1(T_list,1:length(T_list),x1,'previous');
y_recon_2(isnan(y_recon_2)) = 0;

c = [1,1,1];
c_bar = (linspace(0,1,20)').^(0.7)*[1,1,1];

figure(1);
plot(T_list, 3*ones(size(T_list)),'>','markersize',10,'color',0*c);
hold on;
plot([0,length(t)-1],[3,3],'k')
hold on;
plot(T_list_2, 2*ones(size(T_list_2)),'>','markersize',10,'color',0*c);
hold on;
plot([0,length(t)-1],[2,2],'k')
hold on;
plot(T_list_1, 1*ones(size(T_list_1)),'>','markersize',10,'color',0*c);
hold on;
plot([0,length(t)-1],[1,1],'k')
hold on;
yticks([1,2,3]);
yticklabels({'Class 1','Class 2','Whole'});
axis([0,4.5,0,4]);
ytickangle(45)
set(gca,'fontsize',18);

figure(2);
indx_0 = [ones(size(T_list_1));2*ones(size(T_list_2))];
T_list_0 = [T_list_1;T_list_2];
[~,sort_ind] = sort(T_list_0);
sort_indx = indx_0(sort_ind);
plot(T_list(sort_indx == 1), 3*ones(size(T_list(sort_indx == 1))),'>','markersize',10,'color',0*c);
hold on;
plot(T_list(sort_indx == 2), 3*ones(size(T_list(sort_indx == 2))),'>','markersize',10,'color',0*c);
hold on;
plot([0,length(t)-1],[3,3],'k');
hold on;
plot(T_list_2, 2*ones(size(T_list_2)),'>','markersize',10,'color',0*c);
hold on;
plot([0,length(t)-1],[2,2],'k');
hold on;
plot(T_list_1, 1*ones(size(T_list_1)),'>','markersize',10,'color',0*c);
hold on;
plot([0,length(t)-1],[1,1],'k');
hold on;
cate1 = 0;
cate2 = 0;
for i = 1:t(end)
    if sort_indx(i) == 1
        cate1 = cate1+1;
        dir = [T_list(i) - T_list_1(cate1), 2];
        quiver(T_list_1(cate1), 1, dir(1), dir(2),'linewidth',1,'color',0*c);
    else
        cate2 = cate2+1;
        dir = [T_list(i) - T_list_2(cate2), 1];
        quiver(T_list_2(cate2), 2, dir(1), dir(2),'linewidth',1,'color',0*c);
    end    
end
yticks([1,2,3]);
yticklabels({'Class 1','Class 2','Whole'});
axis([0,4.5,0,4]);
ytickangle(45)
set(gca,'fontsize',18);