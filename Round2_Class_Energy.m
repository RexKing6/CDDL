function [gap] = Round2_Class_Energy(Xi,Aa,Xa,Da,trlsa,Ab,Xb,Db,trlsb,index,lambda1,lambda2,lambdaw,eta1,eta2,eta3)

if size(Xa(:,trlsa==index)) == size(Xi)
    Xa(:,trlsa==index)  =  Xi;
    A_index = Aa(:,trlsa==index);
    GAP1 = norm((A_index-Da*Xi),'fro')^2;
    GAP2 = lambda1*sum(abs(Xi(:)));
else
    Xb(:,trlsb==index)  =  Xi;
    A_index = Ab(:,trlsb==index);
    GAP1 = eta1*norm((A_index-Db*Xi),'fro')^2;
    GAP2 = eta1*lambda2*sum(abs(Xi(:)));
end

A=0;
B=0;
tem_XA = Xa;
tem_XB = Xb;
for i=unique(trlsb)
    t_X_icb = tem_XB(:,trlsb==i);
    t_X_ica = tem_XA(:,trlsa==i);
    A=A+t_X_icb*t_X_icb'./size(t_X_icb,2);
    B=B+repmat(t_X_ica,[1,size(t_X_icb,2)])*t_X_icb';
end
w=B/(lambdaw*length(unique(trlsb))/eta2*eye(size(Xb,1))+A);
GAPw = lambdaw*norm(w,'fro')^2;

GAPa3=0;
if size(Xa(:,trlsa==index)) == size(Xi)
    gapa3 = 0;
    gapa4 = 0;
    tem_XA = Xa;
    for i_ca = unique(trlsa)
        t_X_ica = tem_XA(:,trlsa==i_ca);
        gapa4 = gapa4+size(t_X_ica,2)*(mean(t_X_ica,2)-mean(tem_XA,2))'*(mean(t_X_ica,2)-mean(tem_XA,2));
    end
    for i_cb = 1:unique(trlsb)
        t_X_ica = tem_XA(:,trlsa==i_cb);
        t_X_icb = tem_XB(:,trlsb==i_cb);
        gapa3 = gapa3+sum(diag((w*t_X_icb-repmat(t_X_ica,[1,size(t_X_icb,2)]))*(w*t_X_icb-repmat(t_X_ica,[1,size(t_X_icb,2)]))'./size(t_X_icb,2)));
    end
    GAPa3 = eta2*gapa3-eta3*gapa4;
end

gap=GAP1+GAP2+GAPa3+GAPw;