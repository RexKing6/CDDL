function [gap] = Round1_FDL_Energy(Aa,Xa,Fish_par,Fish_ipts)
D = Fish_ipts.D;
lambda1 = Fish_par.tau;
GAP1 = norm((Aa-D*Xa),'fro')^2;
GAP2 = lambda1*sum(abs(Xa(:)));
gap = GAP1+GAP2;