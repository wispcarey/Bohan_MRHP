function y=simu_tMRHP_hist(T,mu,A,omega,h,H)
%simulate multivariate Hawkes process with with exp(t) kernel in time and bivariate Gaussian(x) in space.
%output in chronological order

% Inputs:
% T: scalar time window [0,T] for simulation
% mu: background rates 1*M vector
% A: M-by-M matrix   A(i,j) # of events in j trigger by i
% omega: scalar, time scale
% h: function handle, recursive function
% H: history point process
%    the simulation interval is (max(H(:,2),T]
%
% Outputs 
% simulated point processes data =(u(i), t(i)) i = 1,...,y.n. Here u(i) is the subprocess for event i.
%
% Bohan Chen. 8/24/2020
 
t0 = max(H(:,2));
n0 = size(H,1);
%HP struct
y=struct('n',0,'t',[],'type',[],'lam',[],'father',[],'id',[]);

%step1: background points
bmu=sum(mu);
n_types=length(mu);
n_start=poissrnd(bmu*(T-t0));
bp = struct('n',n_start,'t',[],'type',[],'lam',[]);
bp.t=sort(rand(1,bp.n)*(T-t0))+t0;
bp.type=ones(1,bp.n);
temp=rand(1,bp.n);

for i=1:n_types-1
    bp.type((temp>sum(mu(1:i))/bmu & temp<=sum(mu(1:i+1))/bmu))=i+1;
end

bp.n = bp.n + n0;
bp.type = [H(:,1)',bp.type];
bp.t = [H(:,2)',bp.t];

bp.father=zeros(1,bp.n);
bp.lam = 0;
bp.id = [1:1:bp.n];
w=bp;

%step2: aftershocks
af_types=sum(A,2)';
indx = 1;
while indx <= w.n
    af=struct('n',[],'t',[],'type',[],'id',[],'father',[]);
  
    
    if indx == 1
        w.lam = mu(w.type(1));
    else
        w.lam = mu(w.type(indx)) +...
            h(y.lam)*(A(y.type, w.type(indx)).*omega.*exp(-omega*(w.t(indx)-y.t)));
    end
        
           
    
    af.n=poissrnd(af_types(w.type(indx))*h(w.lam));

    if af.n>0.5
        b1=exprnd(1/omega,1,af.n);
        af.t = b1+w.t(indx);
        less_T = af.t<=T & af.t> t0;
        af.t = af.t(less_T);
        temp=rand(1,af.n);
        temp_type=ones(1,af.n);
        for j=1:n_types-1
            temp_type((temp>sum(A(w.type(indx),1:j))/af_types(w.type(indx)) & temp<=sum(A(w.type(indx),1:j+1))/af_types(w.type(indx))))=j+1;
        end
        af.type = temp_type(less_T);
        af.n = sum(less_T);
        af.father= ones(1,af.n)*indx;
        af.id = w.n + [1:1:af.n];
    end

    %combine
    if af.n>0.5
        w.n = w.n+af.n;
        w.t = [w.t,af.t];
        w.type = [w.type,af.type];
        w.father = [w.father,af.father];
        w.id = [w.id,af.id];
        [w.t,ind]=sort(w.t);
        w.type = w.type(ind);
        w.id = w.id(ind);
        w.father = w.father(ind);        
    end
    
  
    
    y.n = y.n + 1;
    y.type = [y.type; w.type(indx)];
    y.t = [y.t; w.t(indx)];
    y.father = [y.father, w.father(indx)];
    y.id = [y.id, w.id(indx)];
    y.lam = [y.lam, w.lam];
    
    indx = indx+1;
end

end   