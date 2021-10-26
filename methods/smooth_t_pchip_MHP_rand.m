function H = smooth_t_pchip_MHP_rand(H)

% reconstruct the point process from a integer time point process
% details refer to algorithm 4 and 5 
%
% input: H: N*2 or N*4 matrix, integer time point process 
%
% output: H: matrix with the same size as the input H, reconstructed point
% process
%
% Bohan Chen, 8/24/2020

N=length(H);
tdata=H(:,2);
topics=H(:,1);
U = length(unique(topics));

T_list = zeros(size(tdata));
for u = 1:U
    subtdata = smooth_t_pchip_rand(tdata(topics == u));
    T_list(topics == u) = subtdata;
end

[~,indx] = sort(T_list);
tdata = smooth_t_pchip_rand(tdata);

H = [topics(indx),tdata];