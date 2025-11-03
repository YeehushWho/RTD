function [r_sq] = R_square(X,gt)
%   R-SQUARE 
%   R2 = 1-[|X-gt|_F^2/|gt-mean(gt)|_F^2]

SSR = norm(X-gt,"fro")^2;
gt_mean = mean(gt,'all');
SST = norm(gt-gt_mean,"fro")^2;
r_sq = 1 - (SSR/SST);

end

