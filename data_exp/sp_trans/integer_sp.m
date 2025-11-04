function [temp,Re_list,R_sq_list,RE_Integer,R_sq_Integer] = integer_sp(X_gt_int,odd_tensor_en,even_tensor_en,gt_omega_int,gene_num_list_k,gene_list)
         gt_Ind = find(gt_omega_int==1);
         temp = sptrans_integ(odd_tensor_en,even_tensor_en);
         temp_gtind = temp(gt_Ind);
         temp_gtind(temp_gtind<0) = 0;
         temp_gtind(temp_gtind>1) = 1;
         temp(gt_Ind) = temp_gtind;
         [MSE_list,Re_list,R_sq_list,MSE_Integer,RE_Integer,R_sq_Integer] = quality_sp(temp,X_gt_int,gt_omega_int,gt_Ind);
        for i_met = gene_num_list_k
            fprintf("Gene: %s | RE: %5.3f| R-Square: %5.3f \n",gene_list{i_met},Re_list(i_met),R_sq_list(i_met));
        end
end

