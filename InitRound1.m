function [Dict,Drls,coef,CMlabel] = InitRound1(TrainDat,TrainLabel,opts)
TrainDat = TrainDat*diag(1./sqrt(sum(TrainDat.*TrainDat)));

% 初始化字典，每一类分别初始化，初始化函数FDDL_INID
Dict_ini = [];
Dlabel_ini = [];
for ci = unique(TrainLabel)
    cdat          =    TrainDat(:,TrainLabel==ci);
    dict          =    FDDL_INID(cdat,size(cdat,2),opts.wayInit);
    Dict_ini      =    [Dict_ini dict];
    Dlabel_ini    =    [Dlabel_ini repmat(ci,[1 size(dict,2)])];
end
% 这里定义了一个字典标签，可能是字典训练后可能存在相似的几列，在后面是用这个来删减字典的列数

ini_par.tau         =     opts.lambda1;
ini_ipts.D          =     Dict_ini;
coef = zeros(size(Dict_ini,2),size(TrainDat,2));

if size(Dict_ini,1)>size(Dict_ini,2)
  ini_par.c        =    1.05*eigs(Dict_ini'*Dict_ini,1);
else
  ini_par.c        =    1.05*eigs(Dict_ini*Dict_ini',1);
end

% 用初始化的字典来初始化系数
for ci =  unique(TrainLabel)
  fprintf(['Initializing Coef:  Class ' num2str(ci) '\n']);
  ini_ipts.X      =    TrainDat(:,TrainLabel==ci);
  [ini_opts]      =    FDDL_INIC (ini_ipts,ini_par);
  coef(:,TrainLabel ==ci) =    ini_opts.A;
end

% 开始对字典和系数分别迭代更新
Fish_par.dls        =     Dlabel_ini; % 字典标签
Fish_ipts.D         =     Dict_ini;   % 字典
Fish_ipts.trls      =     TrainLabel; % 训练样本标签
Fish_par.tau        =     opts.lambda1;
Fish_par.lambda2    =     opts.lambda2;
Fish_nit            =     1;          % 当前迭代次数
drls                =     Dlabel_ini; % 字典标签

while Fish_nit<=opts.nIter
  if size(Fish_ipts.D,1)>size(Fish_ipts.D,2)
    Fish_par.c        =    1.05*eigs(Fish_ipts.D'*Fish_ipts.D,1);
  else
    Fish_par.c        =    1.05*eigs(Fish_ipts.D*Fish_ipts.D',1);
  end

  % 先更新系数，也是每一类分别更新
  CMlabel = [];
  CoefM = [];
  for ci = unique(TrainLabel)
    fprintf(['Updating coefficients, class: ' num2str(ci) '\n'])
    Fish_ipts.X         =  TrainDat(:,TrainLabel==ci);
    Fish_ipts.A         =  coef;
    Fish_par.index      =  ci;
    % 这里更新系数用的是FDDL里的函数，改了下最后的class Energy，其他直接用上了
    [Copts]             =  Round1_SpaCoef (Fish_ipts,Fish_par);
    coef(:,TrainLabel==ci)    =  Copts.A;
    CMlabel(ci)         =  ci;
    CoefM(:,ci)         =  mean(Copts.A,2);
  end
  % 计算目标函数
  [GAP_coding(Fish_nit)]  =  Round1_FDL_Energy(TrainDat,coef,Fish_par,Fish_ipts);
  
  % 再更新字典
  for ci = unique(TrainLabel)
    fprintf(['Updating dictionary, class: ' num2str(ci) '\n']);
    [Fish_ipts.D(:,drls==ci),Delt(ci).delet]= FDDL_UpdateDi (TrainDat,coef,...
      ci,opts.nClass,Fish_ipts,Fish_par);
  end
  % 计算目标函数
  [GAP_dict(Fish_nit)]  =  Round1_FDL_Energy(TrainDat,coef,Fish_par,Fish_ipts);
  
  % 下面这段应该就是用字典标签来删减字典的列数，当训练样本大到一个程度的时候会删减，目前小样本不会
  newD = []; newdrls = []; newcoef = [];
  for ci = unique(TrainLabel)
    delet = Delt(ci).delet;
    if isempty(delet)
      classD = Fish_ipts.D(:,drls==ci);
      newD = [newD classD];
      newdrls = [newdrls repmat(ci,[1 size(classD,2)])];
      newcoef = [newcoef; coef(drls==ci,:)];
    else
      temp = Fish_ipts.D(:,drls==ci);
      temp_coef = coef(drls==ci,:);
      for temp_i = 1:size(temp,2)
        if sum(delet==temp_i)==0
          newD = [newD temp(:,temp_i)];
          newdrls = [newdrls ci];
          newcoef = [newcoef;temp_coef(temp_i,:)];
        end
      end
    end
  end
  Fish_ipts.D  = newD;
  coef         = newcoef;
  drls         = newdrls;
  Fish_par.dls        =     drls;
  
  Fish_nit = Fish_nit +1;
end

Dict = Fish_ipts.D;
Drls = drls;

if opts.show
  subplot(1,2,1);plot(GAP_coding,'-*');title('GAP_coding');
  subplot(1,2,2);plot(GAP_dict,'-o');title('GAP_dict');
end
return;