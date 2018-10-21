clear,clc;

% load('continuous_AwA.mat')
% 使用AwA连续属性

load('binary_AwA.mat')
% 使用AwA二进制属性

% load('continuous_CUB.mat')
% 使用CUB连续属性

% load('continuous_aPY.mat')
% 使用aPY连续属性

classna = length(unique(trfl))+length(unique(tefl)); % 总的类个数
classnb = length(unique(trfl)); % 训练类个数

trff=[];
trffl=[];

% p为训练样本每类的个数，总共40个训练类
p=10;

for i=unique(trfl)
   tmp=trf(:,trfl==i);
   if(~isempty(tmp))
       trff=[trff tmp(:,1:p)];
   end
   tmp=trfl(:,trfl==i);
   if(~isempty(tmp))
       trffl=[trffl tmp(:,1:p)];
   end
end

% 下面进行归一化，之后PCA降维
trff = trff./repmat(sqrt(sum(trff.^2)), [size(trff, 1),1]);
trff = trff-mean(mean(trff));

[pc,~,latent,~]=princomp(trff');
% PCA维数降为100维
ccc=trff'*pc(:,1:100);
trffpca=ccc';

teff=[];
teffl=[];
% p为测试样本每类的个叔
p=10;
for i=unique(tefl) 
   tmp=tef(:,tefl==i);
   if(~isempty(tmp))
       teff=[teff tmp(:,1:p)];
   end
   tmp=tefl(:,tefl==i);
   if(~isempty(tmp))
       teffl=[teffl tmp(:,1:p)];
   end
end

% 对测试样本归一化，然后PCA降维
teff = teff./repmat(sqrt(sum(teff.^2)), [size(teff, 1),1]);
teff = teff-mean(mean(teff));
% 用上面降维的矩阵降到同一空间
ccc=teff'*pc(:,1:100);
teffpca=ccc';

% 下面是参数设置
opts.show          =   true; % 是否画图
opts.wayInit       =   'random'; % 字典初始化方式，原本FDDL初始化方式是pca，改成了随机初始化，FDDL_INID里有这两种初始化方式
opts.lambda1       =   0.2125; 
opts.lambda2       =   0; % FDDL原本的f(x)系数，但这里我们第一轮没有f(x)，所以取为0
opts.nIter         =   4; % 迭代次数
opts.nClass        =   classna;
[Dict1,Drls1,Coef1] = InitRound1(trc,trcl,opts);

% 因为第一轮都是用同一个函数，所以这里是用了两次lambda1
temp               =   opts.lambda1;
opts.lambda1       =   0.001;
opts.nClass        =   classnb;
[Dict2,Drls2,Coef2] = InitRound1(trffpca,trffl,opts);

opts.lambdaw       =   0.1;
opts.lambda2       =   opts.lambda1;
opts.lambda1       =   temp;
opts.nIter         =   10;
opts.eta1          =   0.1;
opts.eta2          =   100;
opts.eta3          =   100;

A=0;
B=0;
tem_XA = Coef1;
tem_XB = Coef2;
for i=unique(trffl)
    t_X_icb = tem_XB(:,trffl==i);
    t_X_ica = tem_XA(:,trcl==i);
    A=A+t_X_icb*t_X_icb'./size(t_X_icb,2);
    B=B+repmat(t_X_ica,[1,size(t_X_icb,2)])*t_X_icb';
end
w=B/(opts.lambdaw*classnb/opts.eta2*eye(size(Coef2,1))+A);

[Da,Xa,Db,Xb] = Round2(trc,Coef1,Dict1,trcl,classna,Drls1,trffpca,Coef2,Dict2,trffl,classnb,Drls2,opts);

% 训练结束，下面开始测试

A=0;
B=0;
tem_XA = Xa;
tem_XB = Xb;
for i=unique(trffl)
    t_X_icb = tem_XB(:,trffl==i);
    t_X_ica = tem_XA(:,trcl==i);
    A=A+t_X_icb*t_X_icb'./size(t_X_icb,2);
    B=B+repmat(t_X_ica,[1,size(t_X_icb,2)])*t_X_icb';
end
w=B/(opts.lambdaw*classnb/opts.eta2*eye(size(Coef2,1))+A);

lambda   =   opts.lambda2;
nClass   =    classnb;
td1_ipts.tau1 = lambda;
td1_ipts.D    = Db;
if size(Db,1)>=size(Db,2)
   td1_par.eigenv = eigs(Db'*Db,1);
else
   td1_par.eigenv = eigs(Db*Db',1);  
end

% IPM_SC是FDDL自带的一个测试函数
tt_dat = teffpca;
for indTest = 1:size(tt_dat,2)
    fprintf(['Totalnum:' num2str(size(tt_dat,2)) ' Nowprocess:' num2str(indTest) '\n']);
    td1_ipts.y          =      tt_dat(:,indTest);   
    [opts]              =      IPM_SC(td1_ipts,td1_par);
    s                   =      opts.x;
    
    i=1;
    for indClass  =  unique(teffl)
          error(i)=norm(w*opts.x-Xa(:,indClass),2);
          i=i+1;
    end
    % 计算error，选出最小的那个作为结果
    error
    [~,n]=min(error);
    nn=unique(teffl);
    fanal(indTest)=nn(n);
end  

c=0;
d=teffl;
for i=1:p*(classna-classnb)
    if fanal(i)==d(i)
        c=c+1;
    end
end
c/size(tt_dat,2)
% c为正确的个数