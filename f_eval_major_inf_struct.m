function [mse_matrix,eigs_matrix,iterations]=f_eval_major_inf_struct...
    (param,eval)

for ind=1:size(param,2)
    mse_matrix(ind,:)=eval(ind).ESRE_v';
    eigs_matrix(ind,:)=eval(ind).final_eigenvalues';
end

iterations(1,:)=eval(1).iteration_log';

