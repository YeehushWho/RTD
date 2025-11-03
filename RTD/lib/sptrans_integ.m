function [int_tensor] = sptrans_integ(odd_tensor,even_tensor)
[n1,n2,n3] = size(odd_tensor);
int_tensor = zeros(2*n1,2*n2,n3);
odd_row_list = [1:2:(2*n1)-1];
even_row_list =[2:2:(2*n1)];
odd_col_list = [1:2:(2*n2)-1];
even_col_list = [2:2:(2*n2)];
int_tensor(odd_row_list,odd_col_list,:) =odd_tensor;
int_tensor(even_row_list,even_col_list,:) =even_tensor;
end

