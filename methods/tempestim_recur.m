function [u,A,w,Lam,h,lkh,p,kk] = tempestim_recur(H, max_iter, error_bound, fit_mtd)

% EM for temporal multivariate Hawkes process
% based on 'Topic time series analysis of microblogs', IMA Journal of Applied Mathematics, 2016.
% based on Baichuan Yuan's code, https://github.com/ybcmath/MultiSTHP

% Inputs
% H: N*2 matrix, temporal point process  
%    H(i,:) = [u(i),t(i)]
% max_iter: scalar, maximum number of iteration
%           default: 3000
% error_bound: scalar, terminate if the error is less than this bound
%           default: 5e-5
% fit_mtd: 'par' -- fit h parametrically as h(x) = exp(-a*x)
%          'nonpar' -- fit h nonparametrically
%          default: 'par'

% Outputs   -- para
% u: 1*M vector, background rate 
% A: M*M matrix A(i,j) to # of events in j trigger by i, 
%    productivity matrix
% w: time scale, for the time triggering function g(t) = w*exp(-w*t)
% Lam: 1*N vector, conditional intensity
% h: recursive function
% lkh: neg log-likelihood
% p: N-by-N matrix triggering/background prob defined in EM
% 
% Bohan Chen 2020/10/12

if nargin == 1
    max_iter = 3000;
    error_bound = 5e-5;
    fit_mtd = 'par';
end

if nargin == 3
    fit_mtd = 'par';
end

N=length(H);
tdata=H(:,2)';
topics=H(:,1)';
T=H(end,2);
M = length(unique(topics));

%definition: deltat(i,j)=t(j)-t(i)
deltat=triu(bsxfun(@minus,tdata,tdata(tril(ones(N))*ones(N))));

para=rand(1,M^2+1+M);

A=reshape(para(1:M^2),M,M);
u=para(M^2+1:M^2+M);
w=para(M^2+M+1);

lastomega=inf;
lastA=inf;
lastu=inf;
lastl=inf;
  
topic_ind=cell(1,M);
for i=1:M
	topic_ind{i}=(topics==i);
end
inv_t=(T-tdata);

h = @(x)ones(size(x));

for kk =1:max_iter
    tic;    
    % update lambda
    Lam = update_Lam(u,A,w,H,h);
    
    % E steps: estimate P
    [p,lkh] = ExpcstepTemp_recur(u,A,w,H,deltat,h,Lam);
    error=max(max(abs(lastA-A)))+abs(lastomega-w)+max(abs(lastu-u));
    fprintf('iter %d: error = %g, lkh = %g\n', kk, error, lkh);
    aic=2*(M^2+1+M)-2*lkh;
    if  error< error_bound || abs(lkh-lastl)< 1e-6
        kk
        break
    end   
        
    % M steps
    lastomega=w;
    lastA=A;
    lastu=u;
    lastl=lkh;
    At=A';
    diagp=diag(p);
    
    temp_omega1=sum(sum(p.*deltat));
    temp_p=sum(sum(triu(p,1)));
    
    % update omega
    for i = 1:5
        w=temp_p/(temp_omega1+sum(sum(bsxfun(@etimes,h(Lam).*(T-tdata).*exp(-w*(T-tdata)),At(:,topics)))));
    end
    
    % update h, depending on nonpar or par
    if strcmp(fit_mtd, 'nonpar')
        h = update_h(p,A,w,H,Lam);
    elseif strcmp(fit_mtd, 'par')
        h = update_h_exp_new(p,A,w,H,Lam,h);
    else
        error('fit_mtd should be either par or nonpar');
    end
	
    etotimes=exp(-w*inv_t);
    pnodiag=p-diag(diagp);
	
    % update u and A
    u=zeros(1,M);
    for i=1:M
        u(i)=sum(diagp(topic_ind{i}))/T;
        for j=1:M
            A(i,j)=sum(sum(pnodiag(topic_ind{i}, topic_ind{j})))/(sum((1-etotimes(topic_ind{i})).*h(Lam(topic_ind{i}))));
        end
    end
    
    toc;
           
end