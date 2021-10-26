clear;
clc;

addpath('methods');

t = [0,3,5,10,12];
x1 = linspace(0,(length(t)-1),10000);
y_step = interp1(0:(length(t)-1),t,x1,'previous');
%Shape-preserving piecewise cubic interpolation
y_pchip = interp1(0:(length(t)-1),t,x1,'pchip');
[~,indx] = unique(floor(y_pchip));
T_list = x1(indx(2:end))';
y_recon = interp1(T_list,1:length(T_list),x1,'previous');
y_recon(isnan(y_recon)) = 0;
% y_cs = interp1(0:(length(t)-1),t,x1,'spline');

fun = @(x)sign(x).*(abs(x)>1)+x.*(abs(x)<=1);
x_test = linspace(-3,3,200);
y_test = fun(x_test);
x0 = -3:3;
y0 = fun(x0);
y1 = interp1(x0,y0,x_test,'spline');
y2 = interp1(x0,y0,x_test,'pchip');

c = [1,1,1];
c_bar = (linspace(0,1,20)').^(0.7)*[1,1,1];

figure(1);
plot(x0,y0,'.','color',0*c,'markersize',30);
hold on;
plot(x_test,y1,'--','linewidth',2,'color',0*c);
hold on;
plot(x_test,y2,'-','linewidth',2,'color',0*c);
legend('Data','Spline','Akima','location','southeast','Fontsize',24);
set(gca,'fontsize',18);
axis([-3,3,-1.5,1.5]);

figure(2);
plot(x1,y_step,'-','linewidth',2,'color',0*c);
hold on;
plot(x1,y_pchip,'--','linewidth',2,'color',0*c);
legend('Rounded data','Interpolation','Location','northwest','Fontsize',24);
set(gca,'fontsize',18);

figure(3);
plot(x1,y_step,'-','linewidth',2,'color',0*c);
hold on;
plot(x1,y_pchip,'--','linewidth',2,'color',0*c);
hold on;
plot(x1,y_recon,':','linewidth',2,'color',0*c);
legend('Rounded data','Interpolation','Reconstruction','Location','northwest','Fontsize',22);
set(gca,'fontsize',18);

figure(4);
p31 = plot(x1,y_pchip,'--','linewidth',2,'color',0*c);
hold on;
p32 = plot(x1,y_recon,'-','linewidth',2,'color',0*c);
hold on;
[~,indx] = unique(y_recon);
t_recon = x1(indx(2:end));
p33 = plot(t_recon,zeros(size(t_recon)),'.','markersize',30,'color',0*c);
hold on;
for i = 1:length(t_recon)
    plot([t_recon(i),t_recon(i)],[0,i],':','linewidth',2,'color',0*c);
end
legend([p31,p32,p33],...
    {'Interpolation','Reconstruction','Time step {T_i}'},...
    'Location','northwest','Fontsize',21);
set(gca,'fontsize',18);

% rand_T_list = zeros(size(t_recon));
% 
% for i = 1:length(t_recon)
%     if i == 1 && t_recon(i)<=1
%         rand_T_list(i) = t_recon(1)*rand();
%     elseif i == 1 && t_recon(i)>1
%         rand_T_list(i) = (t_recon(i)-floor(t_recon(i)))*rand() + floor(t_recon(i));
%     elseif t_recon(i)-t_recon(i-1) <= 1
%         rand_T_list(i) = (t_recon(i)-t_recon(i-1))*rand() + t_recon(i-1);
%     elseif t_recon(i) == floor(t_recon(i))
%         rand_T_list(i) = rand() + t_recon(i) - 1;
%     else
%         rand_T_list(i) = (t_recon(i)-floor(t_recon(i)))*rand() + floor(t_recon(i));
%     end
% end

load rand_T_10_23

figure(5);
% p41 = plot(x1,y_pchip,'r');
% hold on;
% p42 = plot(x1,y_recon,'black');
% hold on;
% y_rand_t = interp1(0:(length(t)-1),t,rand_T_list,'pchip');

p43 = plot(t_recon, ones(size(t_recon)),'.','markersize',30,'color',0*c);
hold on;
p44 = plot(rand_T_list, zeros(size(rand_T_list)),'.','markersize',10,'color',0*c);
for i = 1:length(t_recon)
    hold on;
    plot([rand_T_list(i),rand_T_list(i)],[0,1],'--','linewidth',0.5,'color',0*c);
end
hold on;
plot([0,4],[0,0],'-','linewidth',0.5,'color',0*c);
hold on;
plot([0,4],[1,1],'-','linewidth',0.5,'color',0*c);

axis([0,4,-0.5,2.5]);
yticks([0,1]);
yticklabels({'t_i','T_i'});

legend([p43,p44],...
    {'Time step {T_i}', 'Reconstructed sequence {t_i}'},...
    'Location','northwest','Fontsize',21);

set(gca,'fontsize',18);

% legend([p41,p42,p43,p44],...
%     {'interpolated curve','reconstructed counting process','time sequence {T_i}','reconstructed time sequence {t_i}'},...
%     'Location','northwest');
