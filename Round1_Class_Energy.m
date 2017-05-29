function [gap] = Round1_Class_Energy(Ai,D,Xi,lambda1)
GAP1  =   norm((Ai-D*Xi),'fro')^2;
GAP2  =   lambda1*sum(abs(Xi(:)));
gap = GAP1+GAP2;