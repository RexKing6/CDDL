function grad = Round2_Gradient_Comp(Xi,ipts,index,obj,eta1,eta2,eta3,lambdaw,newpar,BAI,CJ)
n_d             =      newpar.n_d;                % the sample number of i-th training data
C_j             =      newpar.C_j;
C_line          =      newpar.C_line;   
DD              =      newpar.DD;
DAi             =      newpar.DAi;
BaiBai          =      newpar.BaiBai;
BaiGxi          =      newpar.BaiGxi;
CjCj            =      newpar.CjCj;
m               =      newpar.m;

% 先判断是属性还是特征，再赋予相应的数据计算
if obj=='a'
    classn = ipts.classna;
    X = ipts.Xa;
    trls = ipts.trlsa;
elseif obj=='f'
    classn = ipts.classnb;
    X = ipts.Xb;
    trls = ipts.trlsb;
end

% for k = 1:classn
%     Z(k).Matrix   =   X(:,trls==k)*BAI(k).M-X*C_line+Xi*C_j+X(:,trls==k)*CJ(k).M;
% end

XiT      =   Xi';

tem      =   2*DD*Xi-2*DAi;
grad1    =   tem(:);

grad56   =   0;
X(trls==index) = [];
if obj=='a'
    grad56=-2*eta3/ipts.classna*(1-1/ipts.classna)*((1-1/ipts.classna)*ipts.Xa(trls==index)-mean(X,2));
end
% tem       =  -eta3*(2*BaiBai*XiT-2*BaiGxi);
% grad5     = tem(:);

% grad6 = zeros(size(grad5));
% for k = 1:classn
%     temz  =  Z(k).Matrix;
%     if k~=index
%         tem = -eta3*(2*CjCj*XiT-2*C_j*temz');
%         grad6 = grad6+tem(:);
%     end
% end 
% 
% grad56 = reshape(grad5+grad6,[n_d m])';
% grad56 = grad56(:);

A=0;
B=0;
tem_XA = ipts.Xa;
tem_XB = ipts.Xb;
for i=unique(ipts.trlsb)
    t_X_icb = tem_XB(:,ipts.trlsb==i);
    t_X_ica = tem_XA(:,ipts.trlsa==i);
    A=A+t_X_icb*t_X_icb'./size(t_X_icb,2);
    B=B+repmat(t_X_ica,[1,size(t_X_icb,2)])*t_X_icb';
end
w=B/(lambdaw*length(unique(ipts.classnb))/eta2*eye(size(ipts.Xb,1))+A);

grad2 = 0;
if obj=='a'
    t_X_icb = tem_XB(:,ipts.trlsb==index);
    grad2 = w*t_X_icb-repmat(Xi,[1,size(t_X_icb,2)]);
    grad2 = -2*eta2/size(t_X_icb,2)/ipts.classnb*sum(grad2,2);
elseif obj=='f'
    t_X_ica = tem_XA(:,ipts.trlsa(index));  
    grad2 = w'*(w*Xi-repmat(t_X_ica,[1,size(Xi,2)]));
    grad2 = 2*eta2/size(Xi,2)/ipts.classnb*grad2;
    grad2 = grad2(:);
end
grad = grad1+grad56+grad2;
if obj=='f'
    grad=grad1+grad2./eta1;
end