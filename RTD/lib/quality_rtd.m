function [psnr_val,Rel_error,Mse_val] = quality_rtd(reftensor, result,Ind)
%==========================================================================
% Evaluates the quality assessment indices for two tensors.
%
% Syntax:
%   [psnr_val,Rel_error,Mse_val] = quality_rtd(reftensor, result,Ind)
%
% Input:
%   reftensor - the reference tensor
%   result - the result
%   Ind - the indices of ragged tensor
% Output:
%   psnr_val  - Peak Signal-to-Noise Ratio (higher is better)
%   Rel_error - Relative Error  (lower is better)
%   Mse_val   - Squared Error   (lower is better)
% by Ye-Xun Hu 2024-07-18
%==========================================================================
res_vec = result(Ind);
ref_vec = reftensor(Ind);
psnr_val = calculate_psnr(ref_vec,res_vec,1);
Mse_val    = norm(res_vec-ref_vec,'fro');
Rel_error   = Mse_val/norm(ref_vec,'fro');
function psnr_value = calculate_psnr(original_signal, recovered_signal, max_val)
    mse = mean((original_signal - recovered_signal).^2);
    if mse == 0
        psnr_value = Inf;
    else
        psnr_value = 10 * log10((max_val^2) / mse);
    end
end
end

