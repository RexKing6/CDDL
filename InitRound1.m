function [Dict,Drls,coef,CMlabel] = InitRound1(TrainDat,TrainLabel,opts)
TrainDat = TrainDat*diag(1./sqrt(sum(TrainDat.*TrainDat)));

% ��ʼ���ֵ䣬ÿһ��ֱ��ʼ������ʼ������FDDL_INID
Dict_ini = [];
Dlabel_ini = [];
for ci = unique(TrainLabel)
    cdat          =    TrainDat(:,TrainLabel==ci);
    dict          =    FDDL_INID(cdat,size(cdat,2),opts.wayInit);
    Dict_ini      =    [Dict_ini dict];
    Dlabel_ini    =    [Dlabel_ini repmat(ci,[1 size(dict,2)])];
end
% ���ﶨ����һ���ֵ��ǩ���������ֵ�ѵ������ܴ������Ƶļ��У��ں������������ɾ���ֵ������

ini_par.tau         =     opts.lambda1;
ini_ipts.D          =     Dict_ini;
coef = zeros(size(Dict_ini,2),size(TrainDat,2));

if size(Dict_ini,1)>size(Dict_ini,2)
  ini_par.c        =    1.05*eigs(Dict_ini'*Dict_ini,1);
else
  ini_par.c        =    1.05*eigs(Dict_ini*Dict_ini',1);
end

% �ó�ʼ�����ֵ�����ʼ��ϵ��
for ci =  unique(TrainLabel)
  fprintf(['Initializing Coef:  Class ' num2str(ci) '\n']);
  ini_ipts.X      =    TrainDat(:,TrainLabel==ci);
  [ini_opts]      =    FDDL_INIC (ini_ipts,ini_par);
  coef(:,TrainLabel ==ci) =    ini_opts.A;
end

% ��ʼ���ֵ��ϵ���ֱ��������
Fish_par.dls        =     Dlabel_ini; % �ֵ��ǩ
Fish_ipts.D         =     Dict_ini;   % �ֵ�
Fish_ipts.trls      =     TrainLabel; % ѵ��������ǩ
Fish_par.tau        =     opts.lambda1;
Fish_par.lambda2    =     opts.lambda2;
Fish_nit            =     1;          % ��ǰ��������
drls                =     Dlabel_ini; % �ֵ��ǩ

while Fish_nit<=opts.nIter
  if size(Fish_ipts.D,1)>size(Fish_ipts.D,2)
    Fish_par.c        =    1.05*eigs(Fish_ipts.D'*Fish_ipts.D,1);
  else
    Fish_par.c        =    1.05*eigs(Fish_ipts.D*Fish_ipts.D',1);
  end

  % �ȸ���ϵ����Ҳ��ÿһ��ֱ����
  CMlabel = [];
  CoefM = [];
  for ci = unique(TrainLabel)
    fprintf(['Updating coefficients, class: ' num2str(ci) '\n'])
    Fish_ipts.X         =  TrainDat(:,TrainLabel==ci);
    Fish_ipts.A         =  coef;
    Fish_par.index      =  ci;
    % �������ϵ���õ���FDDL��ĺ���������������class Energy������ֱ��������
    [Copts]             =  Round1_SpaCoef (Fish_ipts,Fish_par);
    coef(:,TrainLabel==ci)    =  Copts.A;
    CMlabel(ci)         =  ci;
    CoefM(:,ci)         =  mean(Copts.A,2);
  end
  % ����Ŀ�꺯��
  [GAP_coding(Fish_nit)]  =  Round1_FDL_Energy(TrainDat,coef,Fish_par,Fish_ipts);
  
  % �ٸ����ֵ�
  for ci = unique(TrainLabel)
    fprintf(['Updating dictionary, class: ' num2str(ci) '\n']);
    [Fish_ipts.D(:,drls==ci),Delt(ci).delet]= FDDL_UpdateDi (TrainDat,coef,...
      ci,opts.nClass,Fish_ipts,Fish_par);
  end
  % ����Ŀ�꺯��
  [GAP_dict(Fish_nit)]  =  Round1_FDL_Energy(TrainDat,coef,Fish_par,Fish_ipts);
  
  % �������Ӧ�þ������ֵ��ǩ��ɾ���ֵ����������ѵ��������һ���̶ȵ�ʱ���ɾ����ĿǰС��������
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