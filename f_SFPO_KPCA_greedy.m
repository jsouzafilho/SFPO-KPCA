
function [PARAM,EVAL]=f_SFPO_KPCA_greedy(mat_data_extraction,mat_data_evaluation,options)
%============================================================================
%
% function [param,eval]=f_SFPO_KPCA_greedy(sample_extraction,sample_evaluation,options)
%
% 	This function implements the SFPO-KPCA algorithm 
%
%	Inputs: mat_data_extraction - input data matrix for KPCA extraction (samples in lines)
%			mat_data_evaluation - input data matrix for KPCA evaluation (samples in lines)
%			options - several configuration parameters of the algorithm (see example script [Sim_f_SFPO_KPCA_grredy.m]).
%
%	Outputs: PARAM(struct) 	- kernel components extracted and related algorithm parameters. 
%            EVAL(struct)	- kernel components related results (ESRE and eigenvalues).		
%
% 	Authors: Prof. J.B.O Souza Filho, P.S.R.Diniz - PEE/POLI/COPPE/EP/FEDERAL
%                   UNIVERSITY OF RIO DE JANEIRO
%   
% E-mails: jbfilho@poli.ufrj.br & diniz@smt.ufrj.br  
%
% Year: 2024 - Version 1.0 (beta)
%============================================================================

    % Start elapsed time monitoring
    t_start=tic;

    % >>> Welcome message
    if (options.PARALLEL.par_welcome_msg_enabled==1)
        fprintf(' ### SFPO-KPCA (2024) - By: JBF & PSR Diniz (PEE/DEL/COPPE/POLI/UFRJ)\n\n - Progress:')
    end
        
	% >>> Set some algorithm flags according to some algorithm options 
	[flag_log_dictionary,flag_eval,flag_remove_mean_FS,...
		flag_log_eigs,flag_dic_replam_en,num_of_alg_steps,max_dict_size,...
            sparse_enabled,num_iter_prog_bar]=f_SKPCA_arg_proc(options);

	% >>> SKPCA initialization procedures
	num_of_ext_samples=size(mat_data_extraction,1); % Data is stored in lines
	dim_of_ext_samples=size(mat_data_extraction,2); 
	    
	% - Define the samples presentation order
	ind_sample_pres_order=f_SKPCA_id_pres_order(num_of_alg_steps,...
        num_of_ext_samples);

    % - Start the dictionary LOG structure (only if enabled)
	if (flag_log_dictionary==1)
		[Dictionary]=f_SKPCA_dictionary_managing_log_start(num_of_alg_steps);
    else
        Dictionary=[];
    end
	
	% - Start the dictionary structure
	[Dictionary]=f_SKPCA_K_m_dictionary_management_start(Dictionary,...
        max_dict_size,dim_of_ext_samples);
	
	% - Get first sample
    iter_cur=1;
    ind_sample_cur=ind_sample_pres_order(iter_cur);
    data_sample_cur=mat_data_extraction(ind_sample_cur,:);
	
    % - Include it into the dictionary LOG (only if enabled)
	if (flag_log_dictionary==1)
		[Dictionary]=f_SKPCA_dictionary_member_inclusion_log(Dictionary,...
			ind_sample_cur,iter_cur);
	end
	
	% - Include first sample into [Dictionary.K_m] and [Dictionary.K_m_inv]
	[Dictionary]=f_SKPCA_K_m_dictionary_member_inclusion([],Dictionary,...
		data_sample_cur,ind_sample_cur,options);
   	
	% - Initialize SKPCA structure
	[SKPCA]=f_SKPCA_param_SKPCA_start(max_dict_size,flag_remove_mean_FS,...
        options);

    % - Update SKPCA structure to consider a new dictionary member  
    [SKPCA]=f_SKPCA_param_SKPCA_update_dictionary_member_inclusion(SKPCA,...
        data_sample_cur,options);
    
    % - Update [KPCA.RB], [KPCA.CB], [KPCA.beta_mean_v] to first sample
    flag_inclusion=1;
    [SKPCA]=f_SKPCA_update_RB_CB_beta_mean(SKPCA,options,...
        flag_inclusion,flag_remove_mean_FS);   

    % - Initialize A matrix	
    SKPCA=f_SKPCA_init_A(SKPCA,Dictionary,max_dict_size,options);

    % - Initialize [EVAL] stucture if enabled (log of ESRE and eigenvalues
    % per iteration PLUS evaluate first sample
    
    if (flag_eval==1)
        % Mount [EVAL] structure
        [EVAL]=f_SKPCA_eval_MSE_eigs_start(flag_log_eigs,options);
        % Evaluate first sample
        [EVAL,~]=f_SKPCA_eval_MSE_eigs(mat_data_evaluation,...
            SKPCA,Dictionary,EVAL,iter_cur,flag_log_eigs,options);
    else
        EVAL=[];
    end
       
    % Sparse flag
    memory_sparse_enabled=sparse_enabled;
    
    % - Iteration loop
    for iter_cur=2:num_of_alg_steps
    
        %-----------------------------------------------------
        % <A> DICTIONARY EVALUATION
        %-----------------------------------------------------
    
        % - Progress bar
        if (rem(iter_cur,num_iter_prog_bar)==0) fprintf('.'); end;
    
        % - Get a new data
        ind_sample_cur=ind_sample_pres_order(iter_cur);
        data_sample_cur=mat_data_extraction(ind_sample_cur,:); % Note data_sample_cur is stored in a line vector

        % - Evaluate ALD criterion about dictionary augmentation
        [SKPCA]=f_KPCA_evaluate_ALD(SKPCA,Dictionary,data_sample_cur,options);
	
        % >>> Set dictionary augmentation FLAGs
        flag_inclusion=0;   
        flag_exclusion=0;
		
        if (SKPCA.ALD.epsilon_current_sq>=options.SFPO_KPCA.max_subspace_error)
           flag_inclusion=1;
           if (options.SFPO_KPCA.dictionary_size_limit==Dictionary.current_size)
                if (flag_dic_replam_en==1)
                    flag_exclusion=1;
                else
                    flag_inclusion=0;   
                end
           end
        end
            
		%-----------------------------------------------------
		% <B> DICTIONARY MEMBER EXCLUSION
		%-----------------------------------------------------
		if (flag_exclusion==1)
 		
			% - Select member to be excluded
			r_m_k=diag(SKPCA.CB);
			[val,ind_exc]=min(r_m_k);

            % Update dictionary log if enabled
			if (flag_log_dictionary==1)
				[Dictionary]=f_SKPCA_dictionary_member_exclusion_log(Dictionary,...
					ind_exc,iter_cur);       
			end

            % Update [K_m] and [K_m_inv] 
			[Dictionary]=f_SKPCA_K_m_dictionary_member_exclusion(Dictionary,...
                ind_exc);
            
			% Update other SFPO parameters 
            [SKPCA]=f_SKPCA_param_SKPCA_update_dictionary_member_exclusion...
                (SKPCA,Dictionary,ind_exc);

        end
   
        %-----------------------------------------------------
        % <C> DICTIONARY MEMBER INCLUSION
        %-----------------------------------------------------
        if (flag_inclusion==1)
        
            % Update dictionary log if enabled
            if (flag_log_dictionary==1)
                [Dictionary]=f_SKPCA_dictionary_member_inclusion_log(Dictionary,...
                    ind_sample_cur,iter_cur);
            end
 
            % Update [K_m] and [K_m_inv] 
        	[Dictionary]=f_SKPCA_K_m_dictionary_member_inclusion(SKPCA,Dictionary,...
            	data_sample_cur,ind_sample_cur,options);
        
            % Update other SFPO parameters
            [SKPCA]=f_SKPCA_param_SKPCA_update_dictionary_member_inclusion(SKPCA,...
                data_sample_cur,options);
				
        end
    
        %-----------------------------------------------------
        % <D> SFPO-KPCA COMPONENT ESTIMATION 
        %-----------------------------------------------------

        % > [CB,RB] estimation
        [SKPCA]=f_SKPCA_update_RB_CB_beta_mean(SKPCA,options,...
            flag_inclusion,flag_remove_mean_FS);   
    
        % > SKPCA component loop
        limit=SKPCA.current_dictionary_size;
        range_col=1:SKPCA.current_dictionary_size;
        
        % - Deflation matrix initialization [Mj]
        Mj=eye(limit);
    
        % - Introduce peridocity of the sparsification procedure
        if (memory_sparse_enabled==1)
           if (mod(iter_cur,options.SPARSE.periodicy)==0)
               sparse_enabled=1;
           else
               sparse_enabled=0;
           end
        end
        
        % - Iteration procedures related to sparse KPCA
        if (sparse_enabled==1)
            
           % - Scaled Cholesky decomposion of CB
            LCT=chol(SKPCA.CB);
            LCT_bar=LCT/max(max(LCT));
        
            % - Producing [U_tilde]
            U=LCT_bar*Dictionary.K_m(range_col,range_col);
        
            D=diag(sqrt(sum(U.*U,2))); 
            D_inv=D^-1;
            U_tilde=U*D_inv;
        end
        
		
        for j=1:options.Extraction.number_of_components
        
            % - FPO kernel 
            a=SKPCA.A_m(j,range_col)';
            r=Dictionary.K_m(range_col,range_col)*a;
            s=SKPCA.CB(range_col,range_col)*r;
        
            alpha_j_tilde=Mj*s;
            alpha_j_bar=alpha_j_tilde/sqrt(alpha_j_tilde'*...
                Dictionary.K_m(range_col,range_col)*alpha_j_tilde);
            
			
			if (sparse_enabled==1)

                % > Sparse part 
                z_j=U*alpha_j_bar;
			
				[b_j_tilde,~,~,~]=...
					f_GREEDY_reg(U_tilde,z_j,options.SPARSE);
				
				b_j=D_inv*b_j_tilde;
			
				b_j_bar=b_j/sqrt(b_j'*...
					Dictionary.K_m(range_col,range_col)*b_j);
            
			else
				b_j_bar=alpha_j_bar;	
			end
			
            % > Deflation
            Mj=Mj-a*r';
               
            SKPCA.A_m(j,1:Dictionary.current_size)=b_j_bar';
            %--------------------------------------------
        end
    
        % Evaluates current component estimate if enabled
        if (flag_eval==1)
            if (mod(iter_cur,options.Evaluate.periodicity)==0)
                [EVAL,~]=f_SKPCA_eval_MSE_eigs(mat_data_evaluation,...
                    SKPCA,Dictionary,EVAL,iter_cur,flag_log_eigs,options);
            end
        end
    end

% - Performs final evaluation
if (flag_eval==1)

    old_choice=options.Evaluate.eigs; % Trick for evaluate eigs 

    options.Evaluate.eigs='YES';
    [EVAL,eigs_vec]=f_SKPCA_eval_MSE_eigs(mat_data_evaluation,...
        SKPCA,Dictionary,EVAL,iter_cur,flag_log_eigs,options);

    EVAL.final_eigenvalues=eigs_vec;
    
    options.Evaluate.eigs=old_choice;
    
end

% Returns final structures
PARAM.SKPCA=SKPCA;
PARAM.Dictionary=Dictionary;

fprintf('\n\n');

% Store time elapsed 
t_start=toc(t_start);
EVAL.sim_time_elapsed=t_start;

end

%======================================================================================
%======================================================================================
%======================================================================================
function [flag_log_dictionary,flag_eval,flag_remove_mean_FS,...
	flag_log_eigs,flag_dic_replam_en,num_of_alg_steps,max_dict_size,...
        sparse_enabled,num_iter_prog_bar]=f_SKPCA_arg_proc(options);

	% Pre-processing of some algorithm options 
	if (strcmp(upper(options.Log.log_dictionary_management),'YES'))
		flag_log_dictionary=1;
	else
		flag_log_dictionary=0;
	end

	if (strcmp(upper(options.Evaluate.active),'YES'))
		flag_eval=1;
	else
		flag_eval=0;
	end

	if (strcmp(upper(options.Extraction.remove_mean_FS),'YES'))
		flag_remove_mean_FS=1;
	else
		flag_remove_mean_FS=0;
	end

	if (strcmp(upper(options.Evaluate.eigs),'YES'))
		flag_log_eigs=1;
	else
		flag_log_eigs=0;
	end

	if (strcmp(upper(options.SFPO_KPCA.enable_dic_replacement),'YES'))
		flag_dic_replam_en=1;
	else
		flag_dic_replam_en=0;
	end
	
	num_of_alg_steps=options.Extraction.number_of_iterations;
    
    if (options.SFPO_KPCA.dictionary_size_limit~=0)
        max_dict_size=options.SFPO_KPCA.dictionary_size_limit;
    else
        max_dict_size=num_of_alg_steps;
    end    
    
    num_iter_prog_bar=ceil(options.Extraction.number_of_iterations/10);
    
    sparse_enabled=options.SPARSE.enable;
    
end
%-----------------------------------------------------------------------------------------
function id_sample_presentation_order=f_SKPCA_id_pres_order(num_of_alg_steps,num_of_samples);
	
	num_col=ceil(num_of_alg_steps/num_of_samples);
	
	dummy_mat=zeros(num_of_samples,num_col);
	
	for ind_c=1:num_col
		dummy_mat(:,ind_c)=randperm(num_of_samples);
	end

	dummy_vec=dummy_mat(:);
	
	id_sample_presentation_order=dummy_vec(1:num_of_alg_steps);

end
%-----------------------------------------------------------------------------------------
function [Dictionary]=f_SKPCA_dictionary_managing_log_start(num_of_alg_steps)

    Dictionary.LOG.inclusion.lentgh=0;
    Dictionary.LOG.inclusion.alg_iteration=zeros(num_of_alg_steps,1); 
	Dictionary.LOG.inclusion.ind_sample_indexes=zeros(num_of_alg_steps,1);
    
    Dictionary.LOG.exclusion.lentgh=0;
    Dictionary.LOG.exclusion.alg_iteration=zeros(num_of_alg_steps,1);
	Dictionary.LOG.exclusion.ind_sample_indexes=zeros(num_of_alg_steps,1);
    
end
%-----------------------------------------------------------------------------------------
function [Dictionary]=f_SKPCA_dictionary_member_inclusion_log(Dictionary,...
        ind_sample,iter_cur);

    Dictionary.LOG.inclusion.lentgh=Dictionary.LOG.inclusion.lentgh+1;

	Dictionary.LOG.inclusion.alg_iteration(Dictionary.LOG.inclusion.lentgh)=...
        iter_cur; 

    Dictionary.LOG.inclusion.ind_sample_indexes(Dictionary.LOG.inclusion.lentgh)=...
        ind_sample;
        
end
%-----------------------------------------------------------------------------------------
function [Dictionary]=f_SKPCA_dictionary_member_exclusion_log(Dictionary,...
        ind_sample,iter_cur);

    Dictionary.LOG.exclusion.lentgh=Dictionary.LOG.exclusion.lentgh+1;

	Dictionary.LOG.exclusion.alg_iteration(Dictionary.LOG.exclusion.lentgh)=...
        iter_cur; 

    Dictionary.LOG.exclusion.ind_sample_indexes(Dictionary.LOG.exclusion.lentgh)=...
        ind_sample;
        
end
%-----------------------------------------------------------------------------------------
function [Dictionary]=f_SKPCA_K_m_dictionary_management_start(Dictionary,...
    max_dict_size,dim_input_data)

	Dictionary.current_size=0;

	Dictionary.K_m=zeros(max_dict_size,max_dict_size);
	
    Dictionary.K_m_inv=zeros(max_dict_size,max_dict_size);
    	
	Dictionary.log_data_sample_index=zeros(max_dict_size,1);
	
	Dictionary.matrix=zeros(dim_input_data,max_dict_size);
	
end
%-----------------------------------------------------------------------------------------
function [SKPCA]=f_SKPCA_param_SKPCA_start(max_dict_size,flag_remove_mean,options)

    % Beta current initialization
    SKPCA.ALD.beta_current=zeros(max_dict_size,1);

    % Kappa current initialization
    SKPCA.ALD.kappa_current=zeros(max_dict_size,1);

    % Mean beta initialization
    SKPCA.beta_mean_v=zeros(max_dict_size,1);

    % Average dictionary A_m
 %   SKPCA.average_dict_A=zeros(max_dict_size,1);
    
    % RB initialization
    SKPCA.RB=zeros(max_dict_size,max_dict_size);
        
    % CB initialization if enabled
    if (flag_remove_mean==1)
		SKPCA.CB=zeros(max_dict_size,max_dict_size);
    end
	
	% Set current dictionary size
	SKPCA.current_dictionary_size=0;
	
end

%-----------------------------------------------------------------------------------------
function [Dictionary]=f_SKPCA_K_m_dictionary_member_inclusion...
    (SKPCA,Dictionary,current_sample,current_sample_index,options);

    %-----------------------------------------------------------------
	% Update dictionary size
	
	n_old_size=Dictionary.current_size;
	Dictionary.current_size=n_old_size+1;
	new_lim=Dictionary.current_size;

    %-----------------------------------------------------------------
	% Update dictionary matrix
    Dictionary.matrix(:,new_lim)=current_sample';

    %-----------------------------------------------------------------
    % [K_m] and [K_m_inv] matrix augmentation

    if (Dictionary.current_size>1)
       
		% K_m augmentation
        range_old=1:n_old_size;
        range_new=1:Dictionary.current_size;
        
        Dictionary.K_m(Dictionary.current_size,Dictionary.current_size)=...
            SKPCA.ALD.kappa_self;
        
        kappa_old=SKPCA.ALD.kappa_current(range_old,1);
        Dictionary.K_m(range_old,Dictionary.current_size)=kappa_old;    % column appending
        Dictionary.K_m(new_lim,range_old)=kappa_old';   % row appending
		Dictionary.K_m(new_lim,new_lim)=SKPCA.ALD.kappa_self;
	
		% K_m_inv matrix augmentation
        beta_cur=SKPCA.ALD.beta_current(range_old,:);
		Dictionary.K_m_inv(range_new,range_new)=...
			[Dictionary.K_m_inv(range_old,range_old) zeros(n_old_size,1);
				zeros(1,n_old_size) 0]+...
					1/SKPCA.ALD.epsilon_current_sq*[-beta_cur; 1]...
                        *[-beta_cur' 1];
        
    else
	
		kappa_self=f_kappa(current_sample,current_sample',options);

        Dictionary.K_m(1,1)=kappa_self;
		
		Dictionary.K_m_inv(1,1)=1/kappa_self;
        
	end

    %-----------------------------------------------------------------
	% Update LOG of indexes of Dictionary samples
	Dictionary.log_data_sample_index(new_lim)=current_sample_index;
	
end
%-----------------------------------------------------------------------------------------
function [Dictionary]=f_SKPCA_K_m_dictionary_member_exclusion...
    (Dictionary,ind_exc);

    %-----------------------------------------------------------------
	% Update dictionary size
	
	old_size=Dictionary.current_size;
	Dictionary.current_size=old_size-1;
	new_size=Dictionary.current_size;
	range_new=1:new_size;

    %-----------------------------------------------------------------
	% Produce indexes of maintained dictionary members
	indices_to_maint=[1:ind_exc-1 (ind_exc+1):old_size];

    %-----------------------------------------------------------------
	% Update dictionary matrix
    Dictionary.matrix(:,range_new)=Dictionary.matrix(:,indices_to_maint);

    %-----------------------------------------------------------------
    % [K_m] and [K_m_inv] matrix shrinkage
	Dictionary.K_m(range_new,range_new)=Dictionary.K_m(indices_to_maint,...
        indices_to_maint);

	% [kappa_old_exchange]
	Dictionary.kappa_old_exchange=Dictionary.K_m(indices_to_maint,ind_exc);
	
	% [K_m_inv]
	K_inv_exchange=Dictionary.K_m_inv([indices_to_maint ind_exc],...
        [indices_to_maint ind_exc]);
    
	norm_u=K_inv_exchange(old_size,old_size);
	u_K_inv=K_inv_exchange(range_new,old_size);
	part_K_inv=K_inv_exchange(range_new,range_new);
	
	Dictionary.K_m_inv(range_new,range_new)=part_K_inv-(u_K_inv*u_K_inv')/norm_u;
	
    %-----------------------------------------------------------------
	% Update LOG of indexes of Dictionary samples
	Dictionary.log_data_sample_index(range_new)=...
		Dictionary.log_data_sample_index(indices_to_maint);
	
end
%-----------------------------------------------------------------------------------------
function [SKPCA]=f_SKPCA_param_SKPCA_update_dictionary_member_inclusion...
    (SKPCA,data_sample_cur,options);
    
    % Update parameters size
    SKPCA.current_dictionary_size=SKPCA.current_dictionary_size+1;
    
    if (SKPCA.current_dictionary_size==1)
        kappa_self=f_kappa(data_sample_cur,data_sample_cur',options);
    else
        kappa_self=SKPCA.ALD.kappa_self;
 	end
    
    %-----------------------------------------------------------------
    % Augmenting [A_m]
    SKPCA.A_m(:,SKPCA.current_dictionary_size)=...
        zeros(options.Extraction.number_of_components,1);
    
	% Update [kappa_k_m_v]
	SKPCA.ALD.kappa_current(SKPCA.current_dictionary_size)=...
        kappa_self;
	
	% Update [beta_k_m_v]
	SKPCA.ALD.beta_current=zeros(SKPCA.current_dictionary_size,1);
	SKPCA.ALD.beta_current(SKPCA.current_dictionary_size,1)=1;
	
	% Update [beta_k_m_v] - Not necessary since pre-started with zero values
	SKPCA.beta_mean_v(SKPCA.current_dictionary_size,1)=0;

	% Update [average_dict_A] - Not necessary since pre-started with zero values
 %   SKPCA.average_dict_A(SKPCA.current_dictionary_size,1)=0;
end
%-----------------------------------------------------------------------------------------
function [SKPCA]=f_SKPCA_param_SKPCA_update_dictionary_member_exclusion...
    (SKPCA,Dictionary,ind_exc);
    
	%-----------------------------------------------------------------
	% This function asssumes that [K_m] and [K_m_inv] are already updated
	
	% - Update dictionary size
    old_size=SKPCA.current_dictionary_size;
    new_size=old_size-1;
    SKPCA.current_dictionary_size=new_size;
	range_col=1:new_size;

	% - Produce indexes of maintained dictionary members
	indices_to_maint=[1:ind_exc-1 (ind_exc+1):old_size];	
                
    % - Readjust A terms due to member exclusion (New K_m_inv*kappa_old_exchange)
    vec_p_old_exchange=Dictionary.K_m_inv(range_col,range_col)*...
        Dictionary.kappa_old_exchange;
            
    A_old_extract=SKPCA.A_m(:,indices_to_maint)';   % |_A_m'_|(1:m-1,:)'
    A_old_vector=SKPCA.A_m(:,ind_exc);              % |_A_m'_|(m,:)'
    SKPCA.A_m(:,range_col)=A_old_extract'+A_old_vector*vec_p_old_exchange';
    
    % - RB cut
    SKPCA.RB(range_col,range_col)=SKPCA.RB(indices_to_maint,indices_to_maint);            
    SKPCA.RB(range_col,old_size)=zeros(new_size,1);        
    SKPCA.RB(old_size,range_col)=zeros(1,new_size);
    SKPCA.RB(old_size,old_size)=0;
    
    % - Mean_vec_med_b
    SKPCA.beta_mean_v(range_col,1)=SKPCA.beta_mean_v(indices_to_maint,1);
    SKPCA.beta_mean_v(old_size,1)=0;
    
    % - Mean average_dict_A
 %   SKPCA.average_dict_A(range_col,1)=SKPCA.average_dict_A(indices_to_maint,1);
 %   SKPCA.average_dict_A(old_size,1)=0;
    
	% - New kappa current
	SKPCA.ALD.kappa_current(range_col,1)=SKPCA.ALD.kappa_current...
		(indices_to_maint,1);
            
    % - New beta_vec
    SKPCA.ALD.beta_current(range_col,1)=Dictionary.K_m_inv(range_col,range_col)*...
		SKPCA.ALD.kappa_current(range_col,1);
	
    % - New epsilon current
    SKPCA.ALD.epsilon_current_sq=SKPCA.ALD.kappa_self-...
        SKPCA.ALD.beta_current(range_col,1)'*SKPCA.ALD.kappa_current(range_col,1);

end
%-----------------------------------------------------------------------------------------
function SKPCA=f_SKPCA_init_A(SKPCA,Dictionary,max_dict_size,options)
    
    % - Initial value for matrix A - assuming only [Dictionary.current_size]
    A_m_part=options.Extraction.std_initial_A_values*...
        randn(options.Extraction.number_of_components,Dictionary.current_size);

	A_m=zeros(options.Extraction.number_of_components,max_dict_size);
	A_m(:,1:Dictionary.current_size)=A_m_part;
	
    % - Impose unitary norm to the previous columns of A in FS 
	K_m=Dictionary.K_m(1:Dictionary.current_size,1:Dictionary.current_size);
	
    for ind_A_m=1:size(A_m,1)
	
		lin_A_m=A_m(ind_A_m,1:Dictionary.current_size);
        
		A_m(ind_A_m,1:Dictionary.current_size)=lin_A_m/sqrt(lin_A_m*K_m*lin_A_m');
    end
    
    SKPCA.A_m=A_m;
end
%-----------------------------------------------------------------------------------------
function [SKPCA]=f_KPCA_evaluate_ALD(SKPCA,Dictionary,data_sample_cur,options);
	
    % - Empirical kernel mapping
    range_col=1:Dictionary.current_size;
    
    SKPCA.ALD.kappa_current(range_col,1)=f_kappa(data_sample_cur,...
        Dictionary.matrix(:,range_col),options);
    
    % - Beta value
    SKPCA.ALD.beta_current(range_col,1)=...
        Dictionary.K_m_inv(range_col,range_col)*SKPCA.ALD.kappa_current(range_col,1);
    
    SKPCA.ALD.beta_current(abs(SKPCA.ALD.beta_current(range_col,1))<10^-9)=0;

    % - Square mapping error
    SKPCA.ALD.kappa_self=f_kappa(data_sample_cur,data_sample_cur',options);
    
    SKPCA.ALD.epsilon_current_sq=SKPCA.ALD.kappa_self-...
        SKPCA.ALD.beta_current(range_col,1)'*SKPCA.ALD.kappa_current(range_col,1);
    	
end
%----------------------------------------------------------------------------------------
function [SKPCA]=f_SKPCA_update_RB_CB_beta_mean(SKPCA,options,...
    flag_inclusion,flag_remove_mean_FS)
    
	var_range=1:SKPCA.current_dictionary_size;
    beta_cur=SKPCA.ALD.beta_current(var_range,1);   
    
    fator=(1-options.SFPO_KPCA.forgetting_factor);
	
    % Infer [average_dict_A_cur] 
%    SKPCA.current_dict_A_rel=sum(abs(SKPCA.A_m(:,var_range)),1);
%    average_dict_A_cur(var_range,1)=SKPCA.current_dict_A_rel;
    
	% Update [beta_mean_v]
    old_beta_mean_v_part=SKPCA.beta_mean_v(var_range,:)*...
		(options.SFPO_KPCA.forgetting_factor);

    SKPCA.beta_mean_v(var_range,:)=old_beta_mean_v_part+...
			fator*beta_cur;
	
	% Update [RB] matrix
	old_RB_part=SKPCA.RB(var_range,var_range)*...
			(options.SFPO_KPCA.forgetting_factor);

    SKPCA.RB(var_range,var_range)=old_RB_part+...
        fator*(beta_cur*beta_cur');
    
    % Update [average_dict_A]
%    old_average_dict_A=SKPCA.average_dict_A(var_range,:)*...
%		(options.SFPO_KPCA.forgetting_factor);

%    SKPCA.average_dict_A(var_range,:)=old_average_dict_A+...
%			fator2*average_dict_A_cur(var_range,:);
    
    % CB matrix if enabled
    if (flag_remove_mean_FS==1)
        SKPCA.CB(var_range,var_range)=SKPCA.RB(var_range,var_range)-...
			SKPCA.beta_mean_v(var_range)*SKPCA.beta_mean_v(var_range)';
    else
        SKPCA.CB(var_range,var_range)=SKPCA.RB(var_range,var_range);
    end

end
%----------------------------------------------------------------------------------------
function [EVAL]=f_SKPCA_eval_MSE_eigs_start(flag_log_eigs,options)

    % Evaluation cycle counter
    EVAL.ind_avaliacoes=0;

    % Set the number of evaluation cycles
    num_eval_cycles=floor(options.Extraction.number_of_iterations/...
        options.Evaluate.periodicity);
    
    % Initialize vector memory of ESRE
    EVAL.ESRE_v=zeros(num_eval_cycles,1);
    
    % Initialize iteration evaluation memory
    EVAL.iteration_log=zeros(num_eval_cycles,1);
    
    % Initialize matrix memory of eigvalues
    if (flag_log_eigs==1)
       EVAL.EIGs_mat=zeros(options.Extraction.number_of_components,...
            num_eval_cycles);
    end
        
end
%----------------------------------------------------------------------------------------
function [EVAL,eigs_vec]=f_SKPCA_eval_MSE_eigs(sample_evaluation,SKPCA,Dictionary,...
    EVAL,iter_cur,flag_log_eigs,options)
    
    % Update evaluation cycle counter
    ind=EVAL.ind_avaliacoes+1;
    EVAL.ind_avaliacoes=ind;
    
    % Define range_col (internal purpose)
    range_col=1:Dictionary.current_size;
    
    % Update iteration log
    EVAL.iteration_log(ind)=iter_cur;
    
    % Infer ESRE and eigenvalues
    [EVAL.ESRE_v(ind),eigs_vec]=...
        f_eval_KPCA(sample_evaluation,SKPCA.A_m(:,range_col)',Dictionary.matrix(:,range_col),...
                     Dictionary.K_m(range_col,range_col),options);
  
    if (flag_log_eigs==1)
        EVAL.EIGs_mat(:,ind)=eigs_vec;
    end
    
end