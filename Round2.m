function [Da,Xa,Db,Xb] = Round2(Aa,Xa,Da,trlsa,classna,drls1,Ab,Xb,Db,trlsb,classnb,drls2,opts)

Aa = Aa*diag(1./sqrt(sum(Aa.*Aa)));
Ab = Ab*diag(1./sqrt(sum(Ab.*Ab)));
Fish_par.tau        =     opts.lambda1;
Fish_par.tau2       =     opts.lambda2;
Fish_par.lambdaw    =     opts.lambdaw;
Fish_par.eta1       =     opts.eta1;
Fish_par.eta2       =     opts.eta2;
Fish_par.eta3       =     opts.eta3;
Fish_nit            =     1;

ipts.Aa                  =     Aa;
ipts.Xa                  =     Xa;
ipts.Da                  =     Da;
ipts.trlsa               =     trlsa;
drlsa                    =     drls1;
ipts.classna             =     classna;
ipts.Ab                  =     Ab;
ipts.Xb                  =     Xb;
ipts.Db                  =     Db;
ipts.trlsb               =     trlsb;
drlsb                    =     drls2;
ipts.classnb             =     classnb;
ipts.drlsa               =     drlsa;
ipts.drlsb               =     drlsb;

while Fish_nit<=opts.nIter
    if size(Da,1)>size(Da,2)
        ca        =    1.05*eigs(Da'*Da,1);
    else
        ca        =    1.05*eigs(Da*Da',1);
    end
    if size(Db,1)>size(Db,2)
        cb        =    1.05*eigs(Db'*Db,1);
    else
        cb        =    1.05*eigs(Db*Db',1);
    end
  
    obj='a';
    % obj判断是属性还是特征

    % 先更新Z1
    for ci = unique(trlsa)
        fprintf(['Updating coefficientsA, class: ' num2str(ci) '\n'])
        ipts.A              =  Aa(:,trlsa==ci);
        Fish_par.index      =  ci;
        % 每一类更新系数
        [Copts]             =  Round2_SpaCoef(ipts,Fish_par,ca,obj);
        Xa(:,trlsa==ci)     =  Copts.A;
    end
    % 计算目标函数
        [GAP_codinga(Fish_nit),~,AA(Fish_nit*4-3),BB(Fish_nit*4-3),CC(Fish_nit*4-3),DD(Fish_nit*4-3),EE(Fish_nit*4-3),FF(Fish_nit*4-3),GG(Fish_nit*4-3)]  =  Round2_FDL_Energy(Aa,Xa,Da,trlsa,Ab,Xb,Db,trlsb,Fish_par);

  % 再更新D1
   for ci = unique(trlsa)
       fprintf(['Updating dictionaryA, class: ' num2str(ci) '\n']);
       Fish_par.dls = drlsa;
       Fish_ipts.D = Da;
       Fish_ipts.trls = trlsa;
       [Da(:,drlsa==ci),Delta(ci).delet]= FDDL_UpdateDi (Aa,Xa,ci,classna,Fish_ipts,Fish_par);
       % 每一类更新字典
   end
   % 计算目标函数
  [GAP_dicta(Fish_nit),~,AA(Fish_nit*4-2),BB(Fish_nit*4-2),CC(Fish_nit*4-2),DD(Fish_nit*4-2),EE(Fish_nit*4-2),FF(Fish_nit*4-2),GG(Fish_nit*4-2)]  =  Round2_FDL_Energy(Aa,Xa,Da,trlsa,Ab,Xb,Db,trlsb,Fish_par);
  
  % 然后再用字典标签删减字典的列数
  newDa = []; newdrlsa = []; newcoefa = [];
  for ci = unique(trlsa)
    delet = Delta(ci).delet;
    if isempty(delet)
      classDa = Da(:,drlsa==ci);
      newDa = [newDa classDa];
      newdrlsa = [newdrlsa repmat(ci,[1 size(classDa,2)])];
      newcoefa = [newcoefa; Xa(drlsa==ci,:)];
    else
      temp = Da(:,drlsa==ci);
      temp_coef = Xa(drlsa==ci,:);
      for temp_i = 1:size(temp,2)
        if sum(delet==temp_i)==0
          newDa = [newDa temp(:,temp_i)];
          newdrlsa = [newdrlsa ci];
          newcoefa = [newcoefa;temp_coef(temp_i,:)];
        end
      end
    end
  end
  
  Da    = newDa;
  Xa    = newcoefa;
  drlsa = newdrlsa;
  ipts.drlsa = drlsa;
  
  % 再更新Z2
  obj='f';
  for ci = unique(trlsb)
    fprintf(['Updating coefficientsF, class: ' num2str(ci) '\n'])
    ipts.A         =  Ab(:,trlsb==ci);
    Fish_par.index      =  ci;
    [Copts]             =  Round2_SpaCoef(ipts,Fish_par,cb,obj);
    Xb(:,trlsb==ci)    =  Copts.A;
  end
  % 计算目标函数
  [GAP_codingb(Fish_nit),~,AA(Fish_nit*4-1),BB(Fish_nit*4-1),CC(Fish_nit*4-1),DD(Fish_nit*4-1),EE(Fish_nit*4-1),FF(Fish_nit*4-1),GG(Fish_nit*4-1)]  =  Round2_FDL_Energy(Aa,Xa,Da,trlsa,Ab,Xb,Db,trlsb,Fish_par);

  % 再更新D2
   for ci = unique(trlsb)
       fprintf(['Updating dictionaryF, class: ' num2str(ci) '\n']);
       Fish_par.dls = drlsb;
       Fish_ipts.D = Db;
       Fish_ipts.trls = trlsb;
       [Db(:,drlsb==ci),Deltb(ci).delet]= FDDL_UpdateDi (Ab,Xb,ci,classnb,Fish_ipts,Fish_par);
   end
   % 计算目标函数
  [GAP_dictb(Fish_nit),~,AA(Fish_nit*4),BB(Fish_nit*4),CC(Fish_nit*4),DD(Fish_nit*4),EE(Fish_nit*4),FF(Fish_nit*4),GG(Fish_nit*4)]  =  Round2_FDL_Energy(Aa,Xa,Da,trlsa,Ab,Xb,Db,trlsb,Fish_par);
  
  % 最后用字典标签删减字典的列数
  newDb = []; newdrlsb = []; newcoefb = [];
  for ci = unique(trlsb)
    delet = Deltb(ci).delet;
    if isempty(delet)
      classDb = Db(:,drlsb==ci);
      newDb = [newDb classDb];
      newdrlsb = [newdrlsb repmat(ci,[1 size(classDb,2)])];
      newcoefb = [newcoefb; Xb(drlsb==ci,:)];
    else
      temp = Db(:,drlsb==ci);
      temp_coef = Xb(drlsb==ci,:);
      for temp_i = 1:size(temp,2)
        if sum(delet==temp_i)==0
          newDb = [newDb temp(:,temp_i)];
          newdrlsb = [newdrlsb ci];
          newcoefb = [newcoefb;temp_coef(temp_i,:)];
        end
      end
    end
  end
  
  Db    = newDb;
  Xb    = newcoefb;
  drlsb = newdrlsb;
  ipts.drlsb = drlsb;
  
  Fish_nit = Fish_nit +1;
end

if opts.show
  subplot(2,2,1);plot(GAP_codinga,'-*');title('GAP codinga');
  subplot(2,2,2);plot(GAP_dicta,'-o');title('GAP dicta');
  subplot(2,2,3);plot(GAP_codingb,'-*');title('GAP codingb');
  subplot(2,2,4);plot(GAP_dictb,'-o');title('GAP dictb');
%     subplot(2,4,1);plot(AA,'-*');
%     subplot(2,4,2);plot(BB,'-*');
%     subplot(2,4,3);plot(CC,'-*');
%     subplot(2,4,4);plot(DD,'-*');
%     subplot(2,4,5);plot(EE,'-*');
%     subplot(2,4,6);plot(FF,'-*');
%     subplot(2,4,7);plot(GG,'-*');
end
return;