function [ESRE_value,EIGs_vector]=f_eval_KPCA(data,A,basis,K_m,options)
%============================================================================
%
% function [ESRE_value,EIGs_vector]=f_eval_KPCA(data,A,basis,K_m,options)
%
% 	Produces ESRE and eigenvalue estimates
%
%	Inputs: data  - input data matrix (samples in lines)
%			A - matrix A 
%			K_m - Gram Matrix of dictionary vectors
%			options.Extraction.kernel_type='RBF' 	- Kernel type  
%			options.Extraction.RBF.kernel_width 	- RBF kernel parameter
%
%	Outputs: ESRE_value 	- estimate of ESRE
%            EIGs_vector - estimate of eigenvalues
%
% Used in FPO-KPCA and FPO-R-KPCA Algorithm (DOI: 10.1109/TSP.2017.275011)	
%
% Authors: Prof. J.B.O Souza Filho, P.S.R.Diniz - PEE/POLI/COPPE/EP/FEDERAL
%                   UNIVERSITY OF RIO DE JANEIRO
%   
% E-mails: jbfilho@poli.ufrj.br & diniz@smt.ufrj.br  
%
% Year: 2017 - Version 1.0
%============================================================================

% - ESRE
kappa_matrix=f_kappa(data,basis,options);

if (strcmp(upper(options.Extraction.remove_mean_FS),'YES'))
    kappa_mean=mean(kappa_matrix,2);
    kappa_centered=kappa_matrix-kappa_mean*ones(1,size(data,1));
else
    kappa_centered=kappa_matrix;
end

Y=A'*kappa_centered;

T=A*Y;

prod_K_ti=K_m*T;

y_sq_mean=mean(sum(Y.*Y,1));
t_sq_mean=mean(sum(T.*prod_K_ti,1));	

ESRE_value=mean(diag(K_m))-2*y_sq_mean+t_sq_mean;

% Includes an extra term in case of centering in FS
if (strcmp(upper(options.Extraction.remove_mean_FS),'YES'))
    K_tst=f_kappa(data,data',options);
    N=size(data,1);
    ESRE_value=ESRE_value-(1/(N^2))*sum(sum(K_tst));
end

% - EIGS_mat
if (strcmp(upper(options.Evaluate.eigs),'YES'))
    EIGs_vector=sum(Y.*Y,2)/size(Y,2);
else
    EIGs_vector=[];
end

