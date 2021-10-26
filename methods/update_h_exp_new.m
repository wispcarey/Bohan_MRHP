function h = update_h_exp_new(p,A,omega,H,Lam,h)

tdata=H(:,2)';
topics=H(:,1)';
T=tdata(end);

al = -log(h(1));

At = A';

B = sum(sum(triu(p,1).*Lam'));
C_fun = @(x)sum(sum(bsxfun(@etimes,(1-exp(-omega*(T-tdata))).*Lam.*exp(-x*Lam),At(:,topics))))-B;
r = rand();
if r<0.05
    al_new = fzero(C_fun,1);
else
    al_new = fzero(C_fun,al);
end

if al_new < 0
    al_new = max(al,0);
end

h = @(x)exp(-al_new*x);