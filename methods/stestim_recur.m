function [A,B,omega,Lam,h,sig,tau,p,pb,mu,kk] = stestim_recur(H, max_iter, error_bound, fit_mtd)

% EM parameter fitting algorithm for multivariate recursive spatiotemporal Hawkes processes 
% Fit a exponential kernel in time and a bivariate Gaussian in space; 
% Fit the recursive function h parametrically as an exponential funciton or nonparametrically as a piecewise constant function
%
% based on 'Topic time series analysis of microblogs', IMA Journal of Applied Mathematics, 2016. 
% based on 'Multivariate Spatiotemporal Hawkes Processes and Network Reconstruction', SIAM Journal on Mathematics of Data Science, 2019. https://github.com/ybcmath/MultiSTHP 
%
% Inputs:
% H: N*4 matrix 
% observed point processes data H(i,:)=[u(i), t(i), x(i), y(i)] i = 1,...,N. Here u(i) is the subprocess for event i.
% max_iter: scalar, maximum number of iteration
%           default: 3000
% error_bound: scalar, terminate if the error is less than this bound
%           default: 5e-5
% fit_mtd: 'par' -- fit h parametrically as h(x) = exp(-a*x)
%          'nonpar' -- fit h nonparametrically
%          default: 'par'
%
% Outputs:
% A: M*M matrix   A(i,j) # of events in j trigger by i
% B: M*M matrix   Similar to A, for background
% omega: scalar time scale, for the time triggering function g(t) = w*exp(-w*t)
% Lam: 1*N vector, conditional intensity
% h: recursive function
% sig=tau: scalar spatial scale 
% p, pd: N-by-N matrix triggering/background prob defined in EM
% mu: 1*U
% kk: number of iterations
%
% Bohan Chen 2020/10/12

if nargin == 1
    max_iter = 1000;
    error_bound = 5e-5;
    fit_mtd = 'par';
end

if nargin == 3
    fit_mtd = 'par';
end

LAT=H(:,4);
LON=H(:,3);
N=length(H);
M=length(unique(H(:,1)));

coord = 0; % coord = 1 if use Great-circle distance 
if coord
    rp_lat=repmat(LAT,1,length(LAT));
    rp_lon=repmat(LON,1,length(LON));
    [arclen] = distance(reshape(rp_lat,1,N^2),reshape(rp_lon,1,N^2),...
    reshape(rp_lat',1,N^2),reshape(rp_lon',1,N^2));
    distance1 = reshape(arclen,N,N);
else
    distance1 = dist2([LAT,LON],[LAT,LON]);
    distance1 = distance1 - spdiags(diag(distance1), 0, N, N);
end

tdata=H(:,2)';
topics=H(:,1)';
M = max(topics);
T=tdata(end);

deltat = triu(bsxfun(@minus, tdata, repmat(tdata', [1 N])));

% Default: Using Random initialization
para=rand(1,2*M^2+2+M);
A=reshape(para(1:M^2),M,M);
A=A/sum(sum(A))*5;
B=reshape(para(M^2+1:2*M^2),M,M);
omega=para(2*M^2+1);
sig=para(2*M^2+2);
tau=para(2*M^2+2);
mu = para(2*M^2+3:2*M^2+2+M);
h = @(x)ones(size(x));

lkh=0;

for kk =1:max_iter
  tic;
  lasts=sig;
  lastomega=omega;
  lastA=A;
  lastB=B;
  lastt=tau;
  lastl=lkh;
  
  % update lambda
  Lam =  update_Lam_st_mu(A,B,omega,sig,tau,H,h);
     

  % E steps: estimate P  
  [pb,p] = Expcstepst_recur(A,B,omega,sig,tau,H,N,deltat,distance1,h,Lam);
    
  % M steps 
  At=A';
  
  % update omega
  for s = 1:5
    omega=sum(sum(p))/(sum(sum(p.*deltat))+sum(sum(bsxfun(@etimes,h(Lam).*(T-tdata).*exp(-omega*(T-tdata)),At(:,topics)))));
  end
  
  % update h, depending on nonpar or par
  if strcmp(fit_mtd, 'nonpar')
      h = update_h(p,A,w,H,Lam);
  elseif strcmp(fit_mtd, 'par')
      h = update_h_exp_new(p,A,omega,H,Lam,h);
  else
      error('fit_mtd should be either par or nonpar');
  end

  % update A
  etotimes=exp(-omega*(T-tdata));
  for i=1:length(unique(topics))
    for j=1:length(unique(topics))
      A(i,j)=sum(sum(p(topics==i, topics==j)))/(sum((1-etotimes(topics==i)).*h(Lam(topics==i))));
    end
  end  
  
  % update B
  for i=1:length(unique(topics))
    for j=1:length(unique(topics))
      B(i,j)=sum(sum(pb(topics==i, topics==j)))/(sum(topics==i));
    end
  end   
    
  % update sigma and tau
  sig=sqrt(sum(sum(p.*distance1+pb.*distance1))/(2*sum(sum(p+pb))));
  tau=sig;

  xx=[reshape(A,1,M^2),reshape(B,1,M^2),omega,sig];
  lkh=log_likesp_recur(xx,H,deltat,distance1,M,h,Lam);
  error=max([max(max(abs(lastA-A))),abs(lastomega-omega),abs(lastt-tau),abs(lasts-sig),max(max(abs(lastB-B)))]);
  
  % update mu, 
  % Here mu work as an average estimation only for the output
  % mu is not included in the iteration 
  % accuracy mu is calculated when updating lambda
  mu = zeros(1,M);
  for i = 1:M
      mu(i) = sum(sum(pb(:,topics == i)))/T;
  end  
   
  % terminal conditions
  if  error< error_bound || abs(lastl-lkh)<error_bound
      break
  end
  
  % nan detection
  if isnan(sum(xx))
      fprintf('NAN!!')
      break
  end
  
  fprintf('iter %d: error = %g, lkh = %g\n', kk, error, lkh);
  toc;
    
end

