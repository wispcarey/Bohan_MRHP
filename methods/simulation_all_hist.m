function v = simulation_all_hist(T,u,A,w,h,Hist,N,bin_size)

v = zeros(N,length(0:bin_size:T)-1);

for i = 1:N
    i
    y=simu_tMRHP_hist(T,u,A,w,h,Hist);
    h1 = histogram(y.t,0:bin_size:T);
    v(i,:) = h1.Values;
    clf;
end