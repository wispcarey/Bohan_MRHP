function Lam = update_Lam_st_mu(A,B,omega,sig,tau,H,h)

tdata=H(:,2)';
topics=H(:,1)';
N=length(tdata);
U=length(unique(topics));
T=tdata(end);
LAT=H(:,4)';
LON=H(:,3)';

Lam = zeros(size(tdata));

mu = cell(U,1);
for i = 1:U
    mu{i} = @(x,y) exp(-((x-LON).^2+(y-LAT).^2)/(2*tau^2))*B(topics, i)/(2*pi*tau^2*T);
end


% Dis_t = distance';

for i = 1:N
    if i == 1
        Lam(i) = mu{topics(i)}(LON(i),LAT(i));
%         Lam(i) = 0.1;
    else
        Lam(i) = mu{topics(i)}(LON(i),LAT(i))+...
            (h(Lam(1:i-1)).*omega.*exp(-omega*(tdata(i)-tdata(1:i-1))).*...
            exp((-(LAT(i)-LAT(1:i-1)).^2-(LON(i)-LON(1:i-1)).^2)/(2*sig^2))/(2*pi*sig^2))*A(topics(1:i-1),topics(i));
    end
end