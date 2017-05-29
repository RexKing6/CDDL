function [gap,w,A,B,C,D,E,F,G] = Round2_FDL_Energy(Aa,Xa,Da,trlsa,Ab,Xb,Db,trlsb,Fish_par)
lambda1 = Fish_par.tau;
lambda2 = Fish_par.tau2;
lambdaw = Fish_par.lambdaw;
eta1    = Fish_par.eta1;
eta2    = Fish_par.eta2;
eta3    = Fish_par.eta3;

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
w=B/(lambdaw*length(unique(trlsb))*eye(size(Xb,1))+A);
GAPw = lambdaw*norm(w,'fro')^2;

gapa3 = 0;
gapa4 = 0;
GAPa1 = norm((Aa-Da*Xa),'fro')^2;
GAPa2 = lambda1*sum(abs(Xa(:)));
tem_XA = Xa;
for i_ca = unique(trlsa)
    t_X_ica = tem_XA(:,trlsa==i_ca);
    gapa4 = gapa4+size(t_X_ica,2)*(mean(t_X_ica,2)-mean(tem_XA,2))'*(mean(t_X_ica,2)-mean(tem_XA,2));
end
for i_cb =unique(trlsb)
    t_X_ica = tem_XA(:,trlsa==i_cb);
    t_X_icb = tem_XB(:,trlsb==i_cb);
    gapa3 = gapa3+sum(diag((w*t_X_icb-repmat(t_X_ica,[1,size(t_X_icb,2)]))*(w*t_X_icb-repmat(t_X_ica,[1,size(t_X_icb,2)]))'./size(t_X_icb,2)));
end
GAPa3 = eta2*gapa3-eta3*gapa4;

GAPb1 = norm((Ab-Db*Xb),'fro')^2;
GAPb2 = lambda2*sum(abs(Xb(:)));
A=GAPa1;
B=sum(abs(Xa(:)));
C=GAPb1;
D=sum(abs(Xb(:)));
E=gapa3;
F=gapa4;
G=norm(w,'fro')^2;
gap=GAPa1+GAPa2+GAPa3+eta1*(GAPb1+GAPb2)+GAPw;