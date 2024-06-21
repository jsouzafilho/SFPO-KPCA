function [ref_par]=f_KPCA_reference_extraction(sample_extraction,options)
%============================================================================
%
% function [ref_par]=f_KPCA_reference_extraction(sample_extraction,options)
%
% 	Produces Standard KPCA components
%
%	Inputs: sample_extraction  - input data matrix (samples in lines)
%			options.Extraction.kernel_type='RBF' 	- Kernel type  
%			options.Extraction.RBF.kernel_width 	- RBF kernel parameter
%			options.Extraction.number_of_components	- Number of kernel components to be extracted
%
%	Outputs: ref_par.A				- Matrix A
%			ref_par.eig_val			- Eigenvalues
%			ref_par.K_uncentered 	- Matrix K uncentered
%			ref_par.K_centered=Kc 	- Matrix K centered in FS 
%			ref_par.mse_ref			- ESRE 
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

num_of_events=size(sample_extraction,1);

K=f_kappa(sample_extraction,sample_extraction',options);

% - Centering kernel matrix in FS if ENABLED!
if (strcmp(upper(options.Extraction.remove_mean_FS),'YES'))
    J = ones(num_of_events,num_of_events)/num_of_events;
    Kc = K - J*K - K*J + J*K*J;
else
    Kc=K;
end

% - Eigendecomposition
option = struct('disp',0);
[eigvector , eigvalue] = eig(Kc);%,options.Extraction.number_of_components,...
%    'la',option);
eigvalue = diag(eigvalue);

% - Sort in decreasing order
[~,sort_ind]=sort(eigvalue,'descend');

eigvalue=eigvalue(sort_ind());
eigvector=eigvector(:,sort_ind); 

% - Reference MSE
mse_ref=(sum(eigvalue((options.Extraction.number_of_components+1):end)...
    /num_of_events));

% - Select eigenpairs of interest
eigvalue=eigvalue(1:options.Extraction.number_of_components);
eigvector=eigvector(:,1:options.Extraction.number_of_components);

% - Cut irrelevant KPCA basis
maxEigValue = max(abs(eigvalue));
eigIdx = find(abs(eigvalue)/maxEigValue < 1e-6);
eigvalue (eigIdx) = [];
eigvector (:,eigIdx) = [];

eigvector=bsxfun(@times,1./sqrt(eigvalue)',eigvector);
eigvalue=eigvalue/num_of_events;

% - Mounting output structure
ref_par.A=eigvector;
ref_par.eig_val=eigvalue;
ref_par.K_uncentered=K;
ref_par.K_centered=Kc;
ref_par.mse_ref=mse_ref;
