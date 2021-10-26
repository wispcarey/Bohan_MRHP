function y=simu_stMRHP(X,Y,T,mu,A,omega,sig,h)
%simulate multivariate recursive Hawkes process with with exp(t) kernel in time and bivariate Gaussian(x) in space.
%output in chronological order
%here sig is var instead of std!

% Inputs:
% X,Y: scalars spatial window [0,X]*[0,Y]
% h: recursive function
% mu: background rates
% A: M-by-M matrix   A(i,j) # of events in j trigger by i
% omega: scalar time scale
% sig: scalar spatial scale
% h: recursive function

% Outputs 
% simulated point processes data =(u(i),t(i),x(i),y(i)) i = 1,...,y.n. Here u(i) is the subprocess for event i. 

%HP struct
y=struct('n',0,'t',[],'type',[],'lam',[],'father',[],'id',[],'lon',[],'lat',[]);

%step1: background points
bmu=sum(mu);
n_types=length(mu);
n_start=poissrnd(bmu*T);
bp = struct('n',n_start,'t',[],'type',[],'lam',[],'lon',[],'lat',[]);
bp.lon=rand(1,bp.n)*X;
bp.lat=rand(1,bp.n)*Y;
bp.t=sort(rand(1,bp.n)*T);
bp.type=ones(1,bp.n);
bp.father=zeros(1,bp.n);
temp=rand(1,bp.n);

for i=1:n_types-1
%     ind=find(temp>mu(i)/bmu & temp<=mu(i+1)/bmu);
    bp.type((temp>sum(mu(1:i))/bmu & temp<=sum(mu(1:i+1))/bmu))=i+1;
end
bp.lam = 0;
bp.id = [1:1:bp.n];
w=bp;

%step2: aftershocks
af_types=sum(A,2)';
indx = 1;
while indx <= w.n
    af=struct('n',[],'t',[],'lon',[],'lat',[],'type',[],'id',[],'father',[]);
  
    
    if indx == 1
        w.lam = mu(w.type(1));
    else
        w.lam = mu(w.type(indx)) +...
            h(y.lam)*(A(y.type, w.type(indx)).*omega.*exp(-omega*(w.t(indx)-y.t)).*...
            exp((-(w.lon(indx)-y.lon).^2-(w.lat(indx)-y.lat).^2)/(2*sig^2))/(2*pi*sig^2));
    end
        
           
    
    af.n=poissrnd(af_types(w.type(indx))*h(w.lam));



    if af.n>0.5
        b1=exprnd(1/omega,1,af.n);
        af.t = b1+w.t(indx);
        less_T = af.t<=T;
        af.t = af.t(less_T);
        xy=mvnrnd([0,0],[sig,0;0,sig],af.n);
        af.lon=xy(:,1)'+w.lon(indx);
        af.lat=xy(:,2)'+w.lat(indx);
        af.lon=af.lon(less_T);
        af.lat=af.lat(less_T);
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
        w.lon = [w.lon, af.lon];
        w.lat = [w.lat, af.lat];
        [w.t,ind]=sort(w.t);
        w.type = w.type(ind);
        w.id = w.id(ind);
        w.father = w.father(ind); 
        w.lon = w.lon(ind);
        w.lat = w.lat(ind);
    end
    
  
    
    y.n = y.n + 1;
    y.type = [y.type; w.type(indx)];
    y.t = [y.t; w.t(indx)];
    y.lon = [y.lon; w.lon(indx)];
    y.lat = [y.lat; w.lat(indx)];
    y.father = [y.father, w.father(indx)];
    y.id = [y.id, w.id(indx)];
    y.lam = [y.lam, w.lam];
    
    indx = indx+1;
end

end   