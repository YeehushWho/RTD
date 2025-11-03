function [MSE_list,Re_list,R_sq_list,MSE_Integer,RE_Integer,R_sq_Integer] = quality_sp(result,X_gt,Omega_gt,Ind_gt)
[n1,n2,n3]=size(result);
MSE_list = zeros(n3,1);
Re_list = zeros(n3,1);
R_sq_list = zeros(n3,1);
for i = 1:n3
    re_slice = result(:,:,i);
    gt_slice = X_gt(:,:,i);
    ind_slice = find(Omega_gt(:,:,i)==1);
    MSE_list(i) = norm(re_slice(ind_slice)-gt_slice(ind_slice),'fro');
    R_sq_list(i) = R_square(re_slice(ind_slice),gt_slice(ind_slice));
    Re_list(i) = MSE_list(i) /norm(gt_slice(ind_slice),'fro');
end
MSE_Integer    = norm(result(Ind_gt)-X_gt(Ind_gt),'fro');
RE_Integer     = MSE_Integer/norm(X_gt(Ind_gt),'fro');
R_sq_Integer = R_square(result(Ind_gt),X_gt(Ind_gt));
end

