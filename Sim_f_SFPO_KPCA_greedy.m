%==================================================================================================
%
% This code exemplifies the use of the function [f_SFPO_KPCA_greedy.m],
%   which implements the algorithm proposed in the paper: A Sparse Fixed-Point Online KPCA
%       Extraction Algorithm
% 
% Authors: Prof. J.B.O Souza Filho, P.S.R.Diniz - PEE/POLI/COPPE/EP/FEDERAL
%                   UNIVERSITY OF RIO DE JANEIRO
%   
% E-mails: jbfilho@poli.ufrj.br & diniz@smt.ufrj.br  
%
% Year: 2024 - Version 1.0 (beta)
%
%==================================================================================================

clear all
close all

%-----------------------------------------------------------------------------------------
% >>> Define the input data <<<
%-----------------------------------------------------------------------------------------

%  Read USPS Input Data 
load('usps_all.mat');

% Select the first 100 samples from the digits 1 to 3.
alg1_mat=data(:,1:100,1)';
alg2_mat=data(:,1:100,2)';
alg3_mat=data(:,1:100,3)';
sample_extraction=[alg1_mat; alg2_mat; alg3_mat];

%-----------------------------------------------------------------------------------------
% >>> Define the parameters related to storing algorithm outcomes <<<
%-----------------------------------------------------------------------------------------
sim_ref='example_simulation';

options.NAMES.text='User-defined text for describing the simulations (arbirary)';

folder='figs/';
mkdir(folder)

options.NAMES.name_arq_data=[ 'res_' sim_ref];
options.NAMES.name_arq_fig=['fig_' sim_ref ];

%-----------------------------------------------------------------------------------------
% >>> Define the data pre-processing <<<
%-----------------------------------------------------------------------------------------

% Normalize input data to 0 to 1 range
sample_extraction=double(sample_extraction)/255;

% Define the evaluation set equal to the extraction set.
sample_evaluation=sample_extraction;

%-----------------------------------------------------------------------------------------
% >>> Defining the SFPO-KPCA simulation parameters <<<
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% A) Set of KPCA related hyperparameters.
%
% 1) Type of kernel used: RBF (Gaussian) - the only one implemented.
options.Extraction.kernel_type='RBF';           

% 2) Parameter "w" from Guassian Kernel: RBF(x)=exp(x-xi)/w.
options.Extraction.RBF.kernel_width=8.0;        

% 3) Enable (YES) or disable (NO) the input data mapping centering into the feature space. 
options.Extraction.remove_mean_FS='no';         

% 4) Define the number of kernel components.
options.Extraction.number_of_components=16;     

%-----------------------------------------------------------------------------------------
% B) Set SFPO-KPCA extraction parameters
%
% 1) Define the total of algorithm iterations.
options.Extraction.number_of_iterations=900;    

% 2) Set the standard deviation for the random initialization of A. 
options.Extraction.std_initial_A_values=0.001;  % 

%-----------------------------------------------------------------------------------------
% C) Set SFPO-KPCA dictionary parameters
%
% 1) Define the forgetting factor for R_\beta estimation (\gamma constant).
options.SFPO_KPCA.forgetting_factor=0.999;        

% 2) Define the maximum input data mapping error in FS (\nil constant).
options.SFPO_KPCA.max_subspace_error=0.001;             

% 3) Set the maximum number of dictionary members (if zero, the dictionary
% is unbounded).
options.SFPO_KPCA.dictionary_size_limit=300;    	

% 4) Choose to allow dictionary member replacement or not.
options.SFPO_KPCA.enable_dic_replacement='no';   

%-----------------------------------------------------------------------------------------
% D) Set SFPO-KPCA performance monitoring parameters
%
% 1) Enable the evaluation of ESRE and eigenvalues at 
%    each [options.Evaluate.ESRE_Eigs_periodicity] iteration.
options.Evaluate.ESRE_Eigs_iterations='yes';             

% 2) Enable eigenvalue estimation at each 
%   [options.Evaluate.ESRE_Eigs_periodicity] iteration.
options.Evaluate.eigs='no';                             % 

% 3) Define the periodicty of ESRE and eigenvalues evaluation.
% In the following case, performance curves 
%   with 50 points are produced.
options.Evaluate.ESRE_Eigs_periodicity=...
    floor(options.Extraction.number_of_iterations/50); 

% 4) Produces (if enabled) in [eval] structure a LOG relative 
%   to the inclusion and exclusion of dictionary members.  
options.Log.log_dictionary_management='no';    

%-----------------------------------------------------------------------------------------
% E) Set SFPO-KPCA greedy regression parameters
%
% 1) Define SFPO-KPCA operation mode (=1, sparsification enabled; =0, no
% sparsification - FPO-mode)
options.SPARSE.enable=1;

% 2) Defines the periodicity to apply sparsification steps (=1, every iteration)
% In the following case, it is assumed the same of periodicity of the evaluation curves. 
options.SPARSE.periodicy=options.Evaluate.ESRE_Eigs_periodicity; 

% 3) Define the sparsification algorithm (1=OMP, 2=GOMP, 3=FBP, 4=SDAR, 5=ASDAR)
options.SPARSE.alg=1;

% 4) Define the hyperparameters related to the sparsification algorithm defined
% in [options.SPARSE.alg]

% > Parameters for OMP (options.SPARSE.alg=1)
options.SPARSE.OMP.lim_cols=49;
options.SPARSE.OMP.tolerance=1E-6;

% > Parameters for GOMP (options.SPARSE.alg=2)
% options.SPARSE.GOMP.inc_factor=2;
% options.SPARSE.GOMP.lim_cols=49;
% options.SPARSE.GOMP.tolerance=1E-5;

% > Parameters for FBP (options.SPARSE.alg=3)
% options.SPARSE.FBP.alpha=5;
% options.SPARSE.FBP.beta=2;
% options.SPARSE.FBP.lim_cols=49;
% options.SPARSE.FBP.tolerance=1E-9;

% > Parameters for SDAR (options.SPARSE.alg=4)
% options.SPARSE.SDAR.T=49; % > Number of columns to be selected by SDAR
% options.SPARSE.SDAR.lim_cols=49;
% options.SPARSE.SDAR.J=10; % > Limit of SDAR iterations

% > Parameters for ASDAR (options.SPARSE.alg=5) 
% options.SPARSE.SDAR.warm_restart.active=0;
% options.SPARSE.ASDAR.T=49; % > Number of columns to be selected by SDAR
% options.SPARSE.ASDAR.lim_cols=49;
% options.SPARSE.ASDAR.J=10; % > Limit of SDAR iterations
% options.SPARSE.ASDAR.range_T=[1:49];%[1 5 10 15 20 25 30 35 40 49];
% options.SPARSE.ASDAR.tolerance=0.0001;


%-----------------------------------------------------------------------------------------
% F) Set the parameters for multiple parallel simulations if this resource
% is enabled.
%
% 1) Enable the parallel simulation mode
options.PARALLEL.par_for_enabled=0;

% 2) Enable SFPO-KPCA welcome message
options.PARALLEL.par_welcome_msg_enabled=1;

% 3) Define the number of MATLAB workers for PARFOR
options.PARALLEL.par_for_workers=2;

%-----------------------------------------------------------------------------------------
% H) Set the total number of simulations
options.PARALLEL.num_simul=1;


%-----------------------------------------------------------------------------------------
% >>> Executes SFPO-KPCA <<<
%-----------------------------------------------------------------------------------------

% Welcome message
fprintf('\n>>> SFPO-KPCA (2024) example code \n\n')

% >>> Produces a reference extraction by classical KPCA
fprintf('- Extacting classical KPCA - start')
[ref_par]=f_KPCA_standard_extraction(sample_extraction,options);
fprintf('- finished.\n\n')


% >>> Set some SFPO-KPCA evaluation parameters according to user choices.
if (strcmp(upper(options.Evaluate.ESRE_Eigs_iterations),'NO')==1)
    % Trick - if 'NO', the values of [ESRE] and [Eigs] are only produced 
    %   in the FIRST and LAST iteration  
    options.Evaluate.active='yes';
    options.Evaluate.periodicity=...
        options.Extraction.number_of_iterations+1;
else
    options.Evaluate.active='yes';
    options.Evaluate.periodicity=...
        options.Evaluate.ESRE_Eigs_periodicity;    
end

if (options.PARALLEL.par_for_enabled==1)

    % >>> Run SPFO-KPCA in PARFOR mode
    mypool=parpool('local');
    t_global_inicio=tic;
        parfor (ind_sim=1:options.PARALLEL.num_simul,options.PARALLEL.par_for_workers);
            fprintf('Sim: [%d]:',ind_sim);
                [param(ind_sim),eval(ind_sim)]=f_SFPO_KPCA_greedy...
                 (sample_extraction,sample_evaluation,options);
            
        end
    options.PARALLEL.elapsed_time=toc(t_global_inicio);
    
    delete(mypool);
else

    % >>> Run SPFO-KPCA in STANDARD mode
    t_global_inicio=tic;
        for (ind_sim=1:options.PARALLEL.num_simul);
            fprintf('Sim: [%d]:',ind_sim);
            [param(ind_sim),eval(ind_sim)]=f_SFPO_KPCA_greedy...
                (sample_extraction,sample_evaluation,options);
        end
    options.PARALLEL.elapsed_time=toc(t_global_inicio);
    
end

%=========================================================================
%
%   Pos-extraction evaluation
%
%=========================================================================

%------------------------------------------------------------
% 1) Compute ESRE and EIGs final values for SFPO-KPCA
%------------------------------------------------------------

% >>> Produces matrices [ESRE_matrix] and [EIGS_matrix] summarizing
% values of ESRE and eigenvectors per iteration
[ESRE_matrix,EIGS_matrix,SFPO_iterations]=f_eval_major_inf_struct(param,eval);

% >>> Compute average values for [ESRE_mat] and [EIGS_matrix] 
Average_SFPO_ESRE=mean(ESRE_matrix,1);
Average_SFPO_EIGS=mean(EIGS_matrix,1);

% >>> Store final values for ESRE, EIGs, and the final ITERATION 
final_SFPO_ESRE=Average_SFPO_ESRE(end);
final_SFPO_EIGS=Average_SFPO_EIGS';
final_SFPO_iter=SFPO_iterations(end);

%------------------------------------------------------------
% 2) Compute final values of ESRE and EIGs for Classical KPCA
%------------------------------------------------------------
options.Evaluate.eigs='YES';
[Classical_KPCA_ESRE,Classical_KPCA_EIGS]=...
    f_eval_KPCA(sample_evaluation,ref_par.A,sample_extraction',...
                    ref_par.K_uncentered,options);                
                        
%------------------------------------------------------------
% 3) Produce a summary of SFPO-KPCA extraction
%------------------------------------------------------------

% >>> Compute the median absolute error between SFPO and Classical KPCA for EIGs
median_EIG_abs_error=median(abs(Classical_KPCA_EIGS-final_SFPO_EIGS));

% >>> Compute the median relative error between SFPO and Classical KPCA for EIGs
median_EIG_rel_error=...
    median(abs((Classical_KPCA_EIGS-final_SFPO_EIGS)./Classical_KPCA_EIGS));

% >>> Exibits a brief summary of SFPO-KPCA execution
fprintf('>>> SFPO-KPCA Execution Summary <<<\n\n')
fprintf('1) Final iteration: [%d]\n',final_SFPO_iter)
if (options.SPARSE.enable==1)
    fprintf('2) Sparsification peridicity: [%d]\n',options.SPARSE.periodicy)
else
    fprintf('2) No sparsification applied: FPO-mode\n')
end
fprintf('3) Number of simulations: [%d]\n',options.PARALLEL.num_simul)
fprintf('4) Average SFPO ESRE: [%2.4f]\n',final_SFPO_ESRE);
fprintf('5) Classical KPCA ESRE: [%2.4f]\n',Classical_KPCA_ESRE');
fprintf('6) Median absolute EIGs error: [%3.4f]\n',median_EIG_abs_error);
fprintf('7) Median relative EIGs error: [%3.4f]\n',median_EIG_rel_error);

%------------------------------------------------------------
% 4) Plot SFPO-KPCA and Classical KPCA Eigs
%------------------------------------------------------------

% - Eigenvalues (PLOT)
figure
semilogy(final_SFPO_EIGS,'b--');
hold on
semilogy(Classical_KPCA_EIGS,'r-.')
legend('SFPO-KPCA','Classical KPCA')
grid on
xlabel('Kernel component')
ylabel('Eigenvalue')
title('Eigenvalue curve')

saveas(gcf,[folder options.NAMES.name_arq_fig '_EIGS.jpeg'],'jpeg')
saveas(gcf,[folder options.NAMES.name_arq_fig '_EIGS.fig'],'fig')

%------------------------------------------------------------
% 5) Plot ESRE curve if log ESRE per iteration enabled
%------------------------------------------------------------

if (strcmp(upper(options.Evaluate.ESRE_Eigs_iterations),'YES')==1)
    figure
    loglog(SFPO_iterations,Average_SFPO_ESRE,'b--');
    hold on
    loglog([SFPO_iterations(1) SFPO_iterations(end)],Classical_KPCA_ESRE*[1 1],'r-')
    grid on
    legend('SFPO-KPCA','Classical KPCA')
    xlabel('Iterations')
    ylabel('ESRE value')
    title('ESRE curve')

    saveas(gcf,[folder options.NAMES.name_arq_fig '_ESRE.jpeg'],'jpeg')
    saveas(gcf,[folder options.NAMES.name_arq_fig '_ESRE.fig'],'fig')
end


% >>> Build a summarizing structure to store all the results.
SIM.param=param;
SIM.eval=eval;
SIM.ref_mse=Classical_KPCA_ESRE;
SIM.ref_eigs=Classical_KPCA_EIGS;
SIM.options=options;
SIM.mat_ESRE=ESRE_matrix;
SIM.EIGS_matrix=EIGS_matrix;

save ([folder options.NAMES.name_arq_data '.mat'],'SIM')


