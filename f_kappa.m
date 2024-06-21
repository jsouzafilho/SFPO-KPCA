function kappa=f_kappa(data,basis,options)
%============================================================================
%
% function kappa=f_kappa(data,basis,options)
%
% 	Produces empirical kernel mapping
%
%	Inputs: data  - input data matrix (samples in lines)
%			basis - matrix whose columns are basis mapping vectors 
%
%			options.Extraction.kernel_type='RBF' 	- Kernel type  
%			options.Extraction.RBF.kernel_width 	- RBF kernel parameter
%
%	Output: kappa - output matrix (the mapped samples are in columns)
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
if(strcmp(options.Extraction.kernel_type,'RBF'))
    
    aa = sum(basis.*basis,1)';
    bb = sum(data.*data,2);
    ab = basis'*data';

    D = bsxfun(@plus,aa,bb') - 2*ab;
    D(D<0) = 0;

    kappa = exp(-D/(2*options.Extraction.RBF.kernel_width^2));

else
    error('>>> Kernel function is not currently implemented!')
end
