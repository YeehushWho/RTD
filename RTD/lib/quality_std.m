function [psnr] = quality_std(imagery1, imagery2)

Nway = size(imagery1);
if length(Nway)>3
    imagery1 = reshape(imagery1,Nway(1),Nway(2),[]);
    imagery2 = reshape(imagery2,Nway(1),Nway(2),[]);
end
Nway = size(imagery1);
psnr = zeros(prod(Nway(3:end)),1);
ssim = psnr;
for i = 1:prod(Nway(3:end))
    psnr(i) = psnr_index(imagery1(:, :, i), imagery2(:, :, i));
end
psnr = mean(psnr);



