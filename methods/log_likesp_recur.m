function l=log_likesp_recur(xx,H,deltat,distance,M,h,Lam) 

% Calculate the negative log likelihood function of a spatiotemporal point process
% with given parameters

% inputs:
% xx: record of A,B,omega,sig,tau,
%     A: xx(1:M^2); B: xx(M^2+1:2M^2)
%     omega: xx(2M^2+1); sig=tau=xx(2M^2+2)
% H: N*4 matrix 
%    observed point processes data H(i,:)=[u(i), t(i), x(i), y(i)] i = 1,...,N. Here u(i) is the subprocess for event i.
% deltat: N*N upper trianguar matrix, deltat(i,j) = t_j-t_i, j>= i
% distance: N*N distance matrix, distance(i,j) is the distance between
%           (x(i),y(i)) and (x(j),y(j))
% M: number of subprocesses (A is M*M)
% h: recursive function
% Lam: 1*N vector, conditional intensity
%
% outputs:
% l: log-likelihood function
%
% Bohan Chen, 8/24/2020

A=reshape(xx(1:M^2),M,M);
B=reshape(xx(M^2+1:2*M^2),M,M);
omega=xx(2*M^2+1);
sig=xx(2*M^2+2);
tau=xx(2*M^2+2);


N=length(H);
tdata=H(:,2)';
topics=H(:,1)';
T=tdata(end);

P=omega.*A(topics,topics).*triu(exp(-omega*deltat)).*triu(exp(-distance/(2*sig.^2)))/(2.*pi.*sig.^2).*(h(Lam)'); 
Pb=exp(-distance/(2.*tau.^2))/(2.*pi.*tau.^2.*T);  
P=P-diag(diag(P));
Pb=Pb-diag(diag(Pb));
  
Pb=B(topics,topics).*Pb;
temp=sum(Pb+P);
for i=1:N
   if temp(i)==0
       Pb(i,i)=1;
   end
end
       

l=-sum(sum(A(topics,:),2)'.*(1-exp(-omega.*(T-tdata)).*h(Lam)))-sum(sum(B(topics,:)))+sum(log(sum(P+Pb)));
l=-l; %negative log-likelihood
