function Lam = update_Lam(mu,A,omega,H,h)

tdata=H(:,2)';
topics=H(:,1)';
N=length(tdata);
U=length(unique(topics));
T=tdata(end);

Lam = zeros(size(tdata));

for i = 1:N
    if i == 1
        Lam(i) = mu(topics(1));
    else
        Lam(i) = mu(topics(i))+...
            (h(Lam(1:i-1)).*omega.*exp(-omega*(tdata(i)-tdata(1:i-1))))*A(topics(1:i-1),topics(i));
    end
end