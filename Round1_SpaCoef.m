function [opts] = Round1_SpaCoef(ipts,par)

par.nIter    =     200;   % maximal iteration number
par.isshow   =     false;
par.citeT    =     1e-6;  % stop criterion
par.cT       =     1e+10; % stop criterion

fish_tau3    =    1;  % parameter of ||A_i-D_iX_i^i||_F^2 in fidelity term
fish_tau2    =    1;  % parameter of ||D_jX_i^j||_F^2 in fidelity term Eq.(4)
m            =    size(ipts.D,2);
drls         =    par.dls;
D            =    ipts.D;
X            =    ipts.X;
A            =    ipts.A;
tau          =    par.tau;
lambda1      =    par.tau;

lambda2      =    par.lambda2;
lambda3      =    par.lambda2;
lambda4      =    par.lambda2;

trls         =    ipts.trls;
classn       =    length(unique(trls));
nIter        =    par.nIter;
c            =    par.c;
sigma        =    c;
tau1         =    tau/2;
index        =    par.index;

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

%preprocessing
Ai                 =          X;  % the i-th training data
Xa                 =          A;  %
Xi                 =          A(:,trls==index);
Xt_now             =          A(:,trls==index);

newpar.n_d          =   size(Ai,2);             % the sample number of i-th training data
n                   =   size(Xa,2);             % the total sample number of training data

for ci = 1:classn
    t_n_d = sum(trls==ci);
    t_b_line_i     =   ones(t_n_d,newpar.n_d)./t_n_d;
    t_c_j = ones(t_n_d,newpar.n_d)./n;
    t_b_angle_i           =   t_b_line_i-t_c_j;   
    CJ(ci).M  = t_c_j;
    BAI(ci).M = t_b_angle_i;   
end

newpar.B_line_i     =   ones(newpar.n_d,newpar.n_d)./newpar.n_d;
newpar.C_j          =   ones(newpar.n_d,newpar.n_d)./n;
newpar.CjCj         =   (newpar.C_j)*(newpar.C_j)';
newpar.C_line       =   ones(n,newpar.n_d)./n;
B_i                 =   eye(newpar.n_d,newpar.n_d)-newpar.B_line_i;
newpar.BiBi         =   B_i*(B_i)';
B_angle_i           =   newpar.B_line_i-newpar.C_j;
newpar.Bai          =   B_angle_i;
newpar.BaiBai       =   B_angle_i*(B_angle_i)';
Xo                  =   Xa;
Xo(:,trls==index)   =   0;
G_X_i               =   Xo*newpar.C_line;
newpar.BaiGxi       =   B_angle_i*(G_X_i)';
newpar.DD           =   D'*D;
newpar.DAi          =   D'*Ai;
Di0                 =   D;
Di0(:,drls~=index)  =   0;
newpar.Di0Di0       =   (Di0)'*Di0;
newpar.Di0Ai        =   (Di0)'*Ai;

newpar.DoiDoi       =   zeros(size(D,2));
for t_i  =  1:classn
    if t_i ~= index
    Doi                 =   D;
    Doi(:,drls~=t_i)      =   0;
    newpar.DoiDoi       =   newpar.DoiDoi+(Doi)'*Doi;
    end 
end

newpar.m            =   m;                 % the number of dictionary column atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xa(:,trls==index)  =  Xi;
xm2       =      Xi;%A(:,trls==index);
xm1       =      Xi;%A(:,trls==index); % now

[gap] = Round1_Class_Energy(Ai,D,xm1,lambda1);
prev_f   =   gap;
ert(1) = gap;
for n_it = 2 : nIter
   Xa(:,trls==index)  =  Xi;
   while for_ever
        % IPM estimate
        % 梯度，第一轮用的是原本FDDL的函数
        grad = FDDL_Gradient_Comp(xm1,Xa,classn,index,...
        lambda2,lambda3,lambda4,fish_tau2,fish_tau3,trls,drls,newpar,...
        BAI,CJ);
    
        v        =   xm1(:)-grad./(2*sigma);
        tem      =   soft(v,tau1/sigma);
        x_temp   =   reshape(tem,[size(D,2),size(xm1,2)]);
        
        if (IST_iters >= 2) | ( TwIST_iters ~= 0)
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
           [gap] = Round1_Class_Energy(Ai,D,xm1,lambda1);

            f = gap;
          
            if (f > prev_f) & (enforceMonotone)
                TwIST_iters   =  0;  % do a IST iteration if monotonocity fails
            else
                TwIST_iters =   TwIST_iters+1; % TwIST iterations
                IST_iters   =    0;
                x_temp      =   xm2;
                if mod(TwIST_iters,10000) ==0
                   c = 0.9*c; 
                   sigma = c;
                end
                break;  % break loop while
            end
        else
          
        [gap] = Round1_Class_Energy(Ai,D,xm1,lambda1);
    
        f   =   gap;
         
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

    citerion      =   abs(f-prev_f)/abs(prev_f);
    if citerion < par.citeT | c > par.cT
%        fprintf('Stop!\n c=%2.2e\n citerion=%2.2e\n',c,citerion);
       break;
    end
    
    xm2           =   xm1;
    xm1           =   x_temp;
    Xt_now        =   x_temp;
    Xi            =   Xt_now; 
    prev_f        =   f;
    ert(n_it)     =   f;
%     fprintf('Iteration:%f  Total gap:%f\n',n_it,ert(n_it-1));
end  

opts.A     =       Xt_now;
opts.ert   =       ert;