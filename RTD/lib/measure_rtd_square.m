function [up,down,left,right] = measure_rtd_square(rtd_mask,padding)
%MEASURE_RTD_SQUARE:get a minimum square which can contain the full
%spatio-irregular tensors with the corrsponding mask.
    function [head,tail] = exam_boud(sum_list,padding)
        sum_list(sum_list>0)=1;
        ind_l = find(sum_list==1);
        head = ind_l(1)-padding;
        tail = ind_l(end)+padding;
    end
row_sum = sum(rtd_mask,[2]);
col_sum = sum(rtd_mask,[1]);
[up,down] = exam_boud(row_sum,padding);
[left,right] = exam_boud(col_sum,padding);


end
