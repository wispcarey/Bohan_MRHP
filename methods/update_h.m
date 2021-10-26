function h = update_h(p,A,omega,H,Lam)

tdata=H(:,2)';
topics=H(:,1)';
N=length(tdata);
U=length(unique(topics));
T=tdata(end);

M = max(Lam);
m = min(Lam);
C = linspace(0,M,30);

% sort_Lam = sort(Lam);
% inds = round(linspace(1,N,10));
% C = sort_Lam(inds);
% C(1) = 0;

C_val = zeros(1,length(C)-1);
At = A';

for i = 1:length(C)-1
    indicator = C(i)<=Lam & Lam<C(i+1);
    nume = sum(sum(triu(p,1).*indicator'));
    deno = sum(sum(bsxfun(@etimes,(1-exp(-omega*(T-tdata))).*indicator,At(:,topics))));
    if deno~=0
        C_val(i) = nume/deno;
    end
end

sum_h = (C(2:end)-C(1:end-1))*C_val';
C_val = C_val/sum_h;

h = @(x)(x>=C(1)).*(x<=C(end)).*...
    reshape(C_val(min(max(sum(x>permute(C,[1,3,2]),3),ones(size(x))),length(C_val)*ones(size(x)))),size(x));