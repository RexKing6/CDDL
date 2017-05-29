function [opts] = FDDL_INIC (ipts,par)
par.initM       =     'zero';  % 初始化方法
par.nIter       =     200;     % 最大迭代次数
par.isshow      =     false;   %
par.twist       =     true;    % 'true': use twist 这个不懂什么含义
par.citeT       =     1e-5;    %  stop criterion
par.cT          =     1e+10;   %  stop criterion
m    =    size(ipts.D,2);
n    =    size(ipts.X,2);

switch lower(par.initM)
    case {'zero'}
        A    =    zeros(m,n);
    case {'transpose'}
        A  =  ipts.D'*ipts.X;
    case {'pinv'}
        A  =  pinv(ipts.D)*ipts.X;
    case {'last'}
        A    =    ipts.last_coef;
    otherwise
        error('Nonknown method!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 下面开始用初始化的字典来迭代更新系数，但涉及到掩模什么的，看不懂
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D        =    ipts.D;
X        =    ipts.X;
tau      =    par.tau;
nIter    =    par.nIter;
c        =    par.c;
sigma    =    c;
tau1     =    tau/2;
B        =    eye(n)-ones(n,n)/n;
% B为文中的eta,艾塔

At_pref   =    A(:);
At_now    =    A(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TWIST parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for_ever           =         1;
IST_iters          =         0;
TwIST_iters        =         0;
sparse             =         1;
verbose            =         1;
enforceMonotone    =         1;
lam1               =         1e-4;   %default minimal eigenvalues
lamN               =         1;      %default maximal eigenvalues
rho0               =         (1-lam1/lamN)/(1+lam1/lamN); 
alpha              =         2/(1+sqrt(1-rho0^2));        %default,user can set
beta               =         alpha*2/(lam1+lamN);         %default,user can set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xm2       =      At_pref;
xm1       =      At_pref;
   
A     =   reshape(At_pref,[m,n]);
gap1  =   norm((X-D*A),'fro')^2;
gap3  =   sum(abs(A(:)));
prev_f   =   gap1+2*tau1*gap3;
   
for n_it = 2 : nIter;
     A     =   reshape(At_now,[m,n]);
     gap1  =   norm((X-D*A),'fro')^2;
     gap3  =   sum(abs(A(:)));
     ert(n_it-1)   =   gap1+2*tau1*gap3;
    while for_ever
     % IPM estimate
        v1   =    [];
        for i  =   1:n
            A     =   reshape(xm1,[m,n]);
            tem1  =   X(:,i)-D*A(:,i);
            tem2  =   D'*tem1;
            v1    =   [v1;tem2];
        end
        A     =   reshape(xm1,[m,n])';
        v     =  xm1+v1/sigma;
        x_temp  =  soft(v,tau1/sigma);
        
        if (IST_iters >= 2) || ( TwIST_iters ~= 0)
            % set to zero the past when the present is zero
            % suitable for sparse inducing priors
            if sparse
                mask    =   (x_temp ~= 0);
                xm1     =   xm1.* mask;
                xm2     =   xm2.* mask;
            end
            % two-step iteration
            xm2    =   (alpha-beta)*xm1 + (1-alpha)*xm2 + beta*x_temp;
            % compute residual
            A     =   reshape(xm2,[m,n]);
            gap1  =   norm((X-D*A),'fro')^2;
            gap3  =   sum(abs(A(:)));
            f   =   gap1+2*tau1*gap3;
          
            if (f > prev_f) && (enforceMonotone)
                TwIST_iters   =  0;  % do a IST iteration if monotonocity fails
            else
                TwIST_iters =   TwIST_iters+1; % TwIST iterations
                IST_iters   =    0;
                x_temp      =   xm2;
                if mod(TwIST_iters,10000) ==0
                   c = 0.9*c; 
                   sigma= c;
                end
                break;  % break loop while
            end
        else
          A     =   reshape(x_temp,[m,n]);
          gap1  =   norm((X-D*A),'fro')^2;
          gap3  =   sum(abs(A(:)));
          f   =   gap1+2*tau1*gap3;

          if f > prev_f
                % if monotonicity  fails here  is  because
                % max eig (A'A) > 1. Thus, we increase our guess
                % of max_svs
                c         =    2*c; 
                sigma     =    c;
                if verbose
%                     fprintf('Incrementing c=%2.2e\n',c);
                end
                if  c > par.cT
                    break;  % break loop while    
                end
                IST_iters = 0;
                TwIST_iters = 0;
           else
                TwIST_iters = TwIST_iters + 1;
                break;  % break loop while
           end
        end
    end

    citerion      =   abs(f-prev_f)/prev_f;
    if citerion < par.citeT || c > par.cT
%        fprintf('Stop!\n c=%2.2e\n citerion=%2.2e\n',c,citerion);
       break;
    end
    
    xm2           =   xm1;
    xm1           =   x_temp;
    At_pref       =   At_now;
    At_now        =   x_temp;
    prev_f        =   f;
    
end

opts.A     =       reshape(At_now,[m,n]);
opts.ert   =       ert;