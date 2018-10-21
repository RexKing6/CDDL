clear,clc;

% load('continuous_AwA.mat')
% ʹ��AwA��������

load('binary_AwA.mat')
% ʹ��AwA����������

% load('continuous_CUB.mat')
% ʹ��CUB��������

% load('continuous_aPY.mat')
% ʹ��aPY��������

classna = length(unique(trfl))+length(unique(tefl)); % �ܵ������
classnb = length(unique(trfl)); % ѵ�������

trff=[];
trffl=[];

% pΪѵ������ÿ��ĸ������ܹ�40��ѵ����
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

% ������й�һ����֮��PCA��ά
trff = trff./repmat(sqrt(sum(trff.^2)), [size(trff, 1),1]);
trff = trff-mean(mean(trff));

[pc,~,latent,~]=princomp(trff');
% PCAά����Ϊ100ά
ccc=trff'*pc(:,1:100);
trffpca=ccc';

teff=[];
teffl=[];
% pΪ��������ÿ��ĸ���
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

% �Բ���������һ����Ȼ��PCA��ά
teff = teff./repmat(sqrt(sum(teff.^2)), [size(teff, 1),1]);
teff = teff-mean(mean(teff));
% �����潵ά�ľ��󽵵�ͬһ�ռ�
ccc=teff'*pc(:,1:100);
teffpca=ccc';

% �����ǲ�������
opts.show          =   true; % �Ƿ�ͼ
opts.wayInit       =   'random'; % �ֵ��ʼ����ʽ��ԭ��FDDL��ʼ����ʽ��pca���ĳ��������ʼ����FDDL_INID���������ֳ�ʼ����ʽ
opts.lambda1       =   0.2125; 
opts.lambda2       =   0; % FDDLԭ����f(x)ϵ�������������ǵ�һ��û��f(x)������ȡΪ0
opts.nIter         =   4; % ��������
opts.nClass        =   classna;
[Dict1,Drls1,Coef1] = InitRound1(trc,trcl,opts);

% ��Ϊ��һ�ֶ�����ͬһ��������������������������lambda1
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

% ѵ�����������濪ʼ����

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

% IPM_SC��FDDL�Դ���һ�����Ժ���
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
    % ����error��ѡ����С���Ǹ���Ϊ���
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
% cΪ��ȷ�ĸ���