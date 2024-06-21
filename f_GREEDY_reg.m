%=====================================================================================
function [b,order,error,options]=f_GREDDY_reg(A,y,options)

	if (options.alg==5)
		[b,order,error,options]=f_ASDAR(A,y,options);
	else
		[b,order,error,options]=f_GREDDY_part(A,y,options);
	end
end
%=====================================================================================
function [b_final,selection_index_order,residue_sq,options]=f_GREDDY_part(A,y,options)

	% 
	% ALGORITHM CHOICE: options.alg= 1 (OMP), 2 (GOMP), 3 (FBP), 4 (SDAR); 
	% INCREMENTAL STEP: options.GOMP.inc_factor; options.FBP.alpha (inc_fator)
	% LIMIT OF COLUMNS: options.{OMP,GOMP,FBP}.lim_cols  
	%

	% Starting procedures for incremental/decremental cholesky decomposition
	[L,ATy,A_col_cur_indexes]=f_LLT_start(A);

	% Apply unitary normalization to A columns
	cols_A_norm_factor=diag(1./sqrt(sum(A.*A,1)));
	
	% Define initial residue
	residue=y;

	% Num_lin, Num_col of A
    [num_cols_A]=size(A,2);
	num_lins_A=1;

	% Define candidate index set
	candidate_col_index=ones(num_cols_A,1);
	
	% Define inclusion factor 
	switch options.alg
		case 1 % OMP
		   %------------------------------------------------
			inclusion_factor=1;
			lim_cols=options.OMP.lim_cols;
			limit_of_operarions=options.OMP.lim_cols;
			tolerance=options.OMP.tolerance;
		   %------------------------------------------------
		case 2 % GOMP
		   %------------------------------------------------
			inclusion_factor=options.GOMP.inc_factor;
			lim_cols=options.GOMP.lim_cols;
			limit_of_operarions=options.GOMP.lim_cols;;
			tolerance=options.GOMP.tolerance;
		   %------------------------------------------------
		case 3 % FBP
		   %------------------------------------------------
			inclusion_factor=options.FBP.alpha;
			exclusion_factor=options.FBP.beta;
			lim_cols=options.FBP.lim_cols;
			limit_of_operarions=options.FBP.lim_cols;
			tolerance=options.FBP.tolerance;
		   %------------------------------------------------
		case 4 % SDAR
		   %------------------------------------------------
			inclusion_factor=options.SDAR.T;
			lim_cols=options.SDAR.lim_cols;
			limit_of_operarions=options.SDAR.J; % [limit of iterations]

			tolerance=0;						% Cancel the norm residue error criterion for quiting

			% Options for warm_restart
			if (options.SDAR.warm_restart.active==1)
			
				L=options.SDAR.warm_restart.L;
				ATy=options.SDAR.warm_restart.ATy;
				A_col_cur_indexes=options.SDAR.warm_restart.A_col_cur_indexes;
								
			else
				% > First SDAR iteration

				% Initial bi_di
				bi_di=A'*y/num_lins_A;
				[val,indexes]=sort(abs(bi_di),'descend');
						
				% Initial active set
				inclusion_set=indexes(1:inclusion_factor,1);
			
				% Update L, ATy and A_col_cur_indexes
				[L,ATy,A_col_cur_indexes]=f_LLT_manage(A,L,ATy,A_col_cur_indexes,y,...
					inclusion_set,1);
					
            end	

	end
	
	% Main Loop
	exit_flag=0;
	num_c_cor=0;
	cont_operations=0;
	while (exit_flag==0)
		
		% Update [cont_operations]
		cont_operations=cont_operations+1;
		
		switch options.alg
			case {1,2,3}
				%===================================================================================
				% OMP / GOMP / FBP - Define the inclusion SET
				
				ind_cand_col=find(candidate_col_index==1);
				residue_A_subset_correlation=residue'*A(:,ind_cand_col)*cols_A_norm_factor(ind_cand_col,ind_cand_col);
				[sorted_res_cor_val,sorted_res_cor_indexes]=sort(abs(residue_A_subset_correlation),'descend');
		
				inclusion_set=ind_cand_col(sorted_res_cor_indexes(1:inclusion_factor));

				% Update L, ATy and A_col_cur_indexes
				[L,ATy,A_col_cur_indexes]=f_LLT_manage(A,L,ATy,A_col_cur_indexes,y,...
					inclusion_set,1);

				% Update [candidate_index]
				candidate_col_index(inclusion_set)=0;
		
				% Update [num_c_cor]
				num_c_cor=num_c_cor+inclusion_factor;

				%===================================================================================
		end
		
		% Solve LL' b = A' y  => Lu=A'y => L'b=u => b for the current A
		[b_cur,residue,residue_sq]=b_LLT_solver(A,L,ATy,A_col_cur_indexes,y);

		%===================================================================================
		if (options.alg==3)
			
			% Evaluate Exclusion (FBP)
						
			[asorted_res_cor_val,asorted_res_cor_indexes]=sort(abs(b_cur),'ascend');
		
			ind_cand_exc_col=A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,1);
	
            exclusion_set=ind_cand_exc_col(asorted_res_cor_indexes(1:exclusion_factor));
		
			[L,ATy,A_col_cur_indexes]=f_LLT_manage(A,L,ATy,A_col_cur_indexes,y,exclusion_set,0);

			% Solve LL' b = A' y  => Lu=A'y => L'b=u => b for the current A
			[b_cur,residue,residue_sq]=b_LLT_solver(A,L,ATy,A_col_cur_indexes,y);

			% Update [candidate_index]
			candidate_col_index(exclusion_set)=1;

			% Update [num_c_cor]
			num_c_cor=num_c_cor-exclusion_factor;
	
			%===================================================================================
		
		elseif (options.alg==4)

			% Evaluate Dual (SDAR)
			ind_cand_col_b=A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,1);
			
			candidate_col_index=ones(num_cols_A,1);
			candidate_col_index(ind_cand_col_b)=0;
			
			ind_cand_col_d=find(candidate_col_index==1);
			
			% Determine di
			di=A(:,ind_cand_col_d)'*residue/num_lins_A;
			
			% Update bi_di
			bi_di=zeros(num_cols_A,1);
			
			bi_di(ind_cand_col_b,1)=b_cur;
			bi_di(ind_cand_col_d,1)=di;
			
			% Evaluate ||bi+di||
			[val,indexes]=sort(abs(bi_di),'descend');
			
			% Identify new active set
			new_inclusion_set=indexes(1:inclusion_factor,1);
			
			% : Compare OLD inclusion set with the new one
			
			% Identify candidates for exclusion
			[val_OLD_NOT_NEW,pos_OLD_NOT_NEW]=setdiff(A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,1),...
				new_inclusion_set);
			
			% Identify candidates for inclusion
			[val_NEW_NOT_OLD,pos_NEW_NOT_OLD]=setdiff(new_inclusion_set,...
					A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,1));
                
			exit_flag=1;
			if (~isempty(val_OLD_NOT_NEW))
			
				% Exclude those eliminated in the new selection
				[L,ATy,A_col_cur_indexes]=f_LLT_manage(A,L,ATy,A_col_cur_indexes,y,...
					val_OLD_NOT_NEW,0);
						
				exit_flag=0;
            end
            
            if (~isempty(val_NEW_NOT_OLD))

				% Include those indicated by the new selection
				[L,ATy,A_col_cur_indexes]=f_LLT_manage(A,L,ATy,A_col_cur_indexes,y,...
					val_NEW_NOT_OLD,1);
					
				exit_flag=0;

			end
				
		end
		%===================================================================================
		
		% Exit criterion 
		% > 1 - residue
		if (residue_sq<tolerance)
			exit_flag=1;
		end

		% > 2 - num of columns
		if (num_c_cor>=lim_cols)
			exit_flag=1;
		end
		
		% > 3 - lim of operations
		if (cont_operations>=limit_of_operarions)
			exit_flag=1;
		end
				
	end
	
	% Function outcomes
	
	% Index of selected columns from A
	selection_index_order=A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,1);

	% Produce sparse b coeficients (note: the unormalization is not necessary ) 
	b_final=zeros(num_cols_A,1);
	b_final(selection_index_order,1)=b_cur;
	
	% In case of SDAR, store relevant content for warm_restart
	if (options.alg==4)
		options.SDAR.warm_restart.L=L;
		options.SDAR.warm_restart.ATy=ATy;
		options.SDAR.warm_restart.A_col_cur_indexes=A_col_cur_indexes;
	end

end
%=============================================================================
function [L,ATy,A_col_cur_indexes]=f_LLT_start(A)

		n_col_A=size(A,2);

		% >>> L
		L.data=zeros(n_col_A,n_col_A);
		L.n_lin=0;
		L.n_col=0;

		% >>> ATy
		ATy.data=zeros(n_col_A,1);
		ATy.n_lin=0;
		ATy.n_col=1;

		% >>> Acur
		A_col_cur_indexes.data=zeros(n_col_A,1);
		A_col_cur_indexes.n_lin=0;
		A_col_cur_indexes.n_col=1;

end
%=========================================================================================
function [L,ATy,A_col_cur_indexes]=f_LLT_manage(A,L,ATy,A_col_cur_indexes,y,vec_indices,option)

	num_op=length(vec_indices);

	if (option==1) % Inclusion
		
		for ind_c=1:num_op
			[L,ATy,A_col_cur_indexes]=f_LLT_inc(A,L,ATy,A_col_cur_indexes,y,vec_indices(ind_c));
		end
	else

		for ind_c=1:num_op
			[L,ATy,A_col_cur_indexes]=f_LLT_dec(L,ATy,A_col_cur_indexes,vec_indices(ind_c));
		end
	end
end
%=========================================================================================
function [L,ATy,A_col_cur_indexes]=f_LLT_inc(A,L,ATy,A_col_cur_indexes,y,col_A_index)

		% Retrive A column defined by [col_A_index]
		col_A=A(:,col_A_index);

		% Inclusion
		if (L.n_lin==0)
		
				% First iteration
				L.data(1,1)=sqrt(col_A'*col_A);
				L.n_lin=1;
				L.n_col=1;
				
				ATy.data(1,1)=col_A'*y;
				ATy.n_lin=1;
				ATy.n_col=1;
				
				A_col_cur_indexes.data(1,1)=col_A_index;
				A_col_cur_indexes.n_lin=1;
				
			else
				% Other iterations
				
				%-------------------------------------------------
				% >>> L update
				
				% Square norm of col_A
				sq_norm_col_A=col_A'*col_A;
				
				% Proj A_(k-1)'a_k	
				proj_Ak_old_ak=A(:,A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,1))'*col_A;
				
				% Estimate l_k
				l_k=L.data(1:L.n_lin,1:L.n_col)\proj_Ak_old_ak;

				% Estimate l_kk
				l_kk=sqrt(sq_norm_col_A-l_k'*l_k);
				
				% Build new L matrix L_new=[L_old 0; l_k' l_kk]
				L.n_lin=L.n_lin+1;
				L.data(L.n_lin,1:L.n_col)=l_k';
				L.n_col=L.n_col+1;
				L.data(L.n_col,L.n_col)=l_kk;
				
				%-------------------------------------------------
				% >>> ATy update (projection of col_A into y)
				ATy.n_lin=ATy.n_lin+1;
				ATy.data(ATy.n_lin,:)=col_A'*y;
				
				%-------------------------------------------------
				% >>> A_cur update 
				A_col_cur_indexes.n_lin=A_col_cur_indexes.n_lin+1;
				A_col_cur_indexes.data(A_col_cur_indexes.n_lin,1)=col_A_index;
				%-------------------------------------------------
				
				
		end
end
%=========================================================================================
function [L,ATy,A_col_cur_indexes]=f_LLT_dec(L,ATy,A_col_cur_indexes,ind_col_exc_A)

        %---------------------------------------------------
        % >>> Identifies proper ind_col_exc in list
        ind_col_exc=find(A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,1)==ind_col_exc_A);

		%---------------------------------------------------
		% >>> L update
		
		% Mount permutation matrix over A. The column [ind_col_exc] will turn the last one.
		n_col_A_sel=A_col_cur_indexes.n_lin;
		
		P=diag(ones(n_col_A_sel,1));
		col_indexes=[1:(ind_col_exc-1) (ind_col_exc+1):n_col_A_sel ind_col_exc];
		P=P(:,col_indexes);
		
		% Apply P over L 
		L_perm=P'*L.data(1:L.n_lin,1:L.n_col);
		
		% Set L_perm to a diagonal form by an orthogonal transform: L_perm=R'*Q' => QR=L_perm'
		[Q,R]=qr(L_perm');
		
		% Set new L
		L.n_lin=L.n_lin-1;
		L.n_col=L.n_col-1;

		L.data(1:L.n_lin,1:L.n_col)=R(1:L.n_col,1:L.n_lin)';

		%---------------------------------------------------
		% >>> ATy update
		ATy.n_lin=ATy.n_lin-1;
		ATy.data(1:ATy.n_lin,1)=ATy.data(col_indexes(1:ATy.n_lin),1);

		%---------------------------------------------------
		% >>> A_cur update
		A_col_cur_indexes.n_lin=A_col_cur_indexes.n_lin-1;
		A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,:)=...
        A_col_cur_indexes.data(col_indexes(1:A_col_cur_indexes.n_lin),1);
		
		%---------------------------------------------------
end
%=============================================================================
function [b_cur,residue,residue_sq]=b_LLT_solver(A,L,ATy,A_col_cur_indexes,y)

		% Solve LL' b = A' y  => Lu=A'y => L'b=u => b for the current A
		u_cur=L.data(1:L.n_lin,1:L.n_col)\ATy.data(1:ATy.n_lin,1);
		b_cur=L.data(1:L.n_lin,1:L.n_col)'\u_cur;

		% Determine new residue
		residue=y-A(:,A_col_cur_indexes.data(1:A_col_cur_indexes.n_lin,1))*b_cur;
		residue_sq=residue'*residue;

end
%=====================================================================================
function [b,order,error,options]=f_ASDAR(A,y,options)

	% Starting procedures
	options.alg=4;

	options.SDAR.lim_cols=options.ASDAR.lim_cols;
	options.SDAR.J=options.ASDAR.J;

	% ASDAR loop
	for ind_s=1:length(options.ASDAR.range_T)
    
		if (ind_s==1)
			options.SDAR.warm_restart.active=0;
		else
			options.SDAR.warm_restart.active=1;
		end

		options.SDAR.T=options.ASDAR.range_T(ind_s);

		[b,order,error,options]=f_GREDDY_part(A,y,options);

		if (error<options.ASDAR.tolerance)
			break;
		end
	end
end
%=====================================================================================