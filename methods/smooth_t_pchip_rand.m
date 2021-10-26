function [rand_T_list, T_cum] = smooth_t_pchip_rand(T_list)
%
% reconstruct the time sequence from a given integer time sequence
% details refer to algorithm 4
%
% input: 
% T_list: integer time sequence
% 
% outputs:
% rand_T_list:reconstructed time sequence
% T_cum: cumulative number of cases
%
% Bohan Chen, 8/24/2020

T_list = floor(T_list);
T = max(T_list)+1;
N = length(T_list);

T_cum = zeros(1,T+1);

for i = 1:T
    T_cum(i+1) = sum(T_list < i);
end

x1 = linspace(0,T,20000000);
y1 = interp1(0:T,T_cum,x1,'pchip');
[~,indx] = unique(floor(y1));
T_list = x1(indx(2:end))';
int_indx = abs(T_list - round(T_list))<1e-4;
T_list(int_indx) = round(T_list(int_indx));
rand_T_list = zeros(size(T_list));

for i = 1:length(T_list)
    if i == 1 && T_list(i)<=1
        rand_T_list(i) = T_list(1)*rand();
    elseif i == 1 && T_list(i)>1
        rand_T_list(i) = (T_list(i)-floor(T_list(i)))*rand() + floor(T_list(i));
    elseif T_list(i)-T_list(i-1) <= 1
        rand_T_list(i) = (T_list(i)-T_list(i-1))*rand() + T_list(i-1);
    elseif T_list(i) == floor(T_list(i))
        rand_T_list(i) = rand() + T_list(i) - 1;
    else
        rand_T_list(i) = (T_list(i)-floor(T_list(i)))*rand() + floor(T_list(i));
    end
end