clear;
rng(5);
warning('off', 'all');
%% Data selection: 1. pavia  |　2. dcmall
data_selection =  1;

if data_selection==1
    data_name = 'pavia';
else
    data_name = 'dcmall'; 
end
nl=20;
slice_list = [65,45,15];
Nclust =100;

load(['./data_exp/dn/hsi_' data_name '.mat']);
[n1,n2,n3] = size(X_gt_int);
X_ob_int = noise_tensor+X_gt_int;
label_list = unique(sup_labels(:));
gt_lab_int = sup_labels;
label_len = length(label_list);
gt_omega_int = ones(size(X_gt_int));
gt_Ind_int = find(gt_omega_int==1);
Method_list = {'Observed','CPD-PAM','RTD'};
method_len  = length(Method_list);
EN_CP =1;
EN_RTD = 1;

Result_final = cell(method_len,1); 
PSNR_final  = zeros(method_len,1);
RelErr_final = zeros(method_len,1);  
Runtime_final = zeros(method_len,1); 
for i_noise = 1
    fprintf("===============================noise level: %4.2f==================================\n",nl);
    %% observed tensor
    kkk=1;
    fprintf("===================%s==================\n",Method_list{kkk});
    Result_final{kkk,i_noise} = X_ob_int;
    [PSNR_final(kkk,i_noise)] = quality_std(X_gt_int*255,X_ob_int*255);
    RMSE_i   = norm(X_ob_int-X_gt_int,'fro');
    RelErr_final(kkk,i_noise)     = RMSE_i/norm(X_gt_int,'fro');
    fprintf("%s | nl: %4.2d |PSNR: %5.3f | RE: %6.3f \n",Method_list{kkk},nl,PSNR_final(kkk,i_noise),RelErr_final(kkk,i_noise));
    %% CPD-PAM

    kkk=2;
    if EN_CP
        fprintf("===================%s==================\n",Method_list{kkk});
        tol   = 1e-3;
        max_iter    = 2000;
        jjj = 0;
        r = 30;
        fprintf("-------------------Rank:%d---------------\n",r);
        for rho = 10.^[3] 
            jjj = jjj+1;
            start_time = tic;
            ob_Omega = ones(size(X_gt_int));
            [X_CP, res2,iter_rtd,~,~,~] = CP_proximal_dn(X_ob_int,rho,r,tol,max_iter,gt_Ind_int,X_gt_int,gt_omega_int);
            Runtime_final(kkk,i_noise) = toc(start_time);
            Result_final{kkk,i_noise}= X_CP;
            [PSNR_final(kkk,i_noise)]= quality_std(X_gt_int*255,X_CP*255);
            MSE_Integer    = norm(X_CP-X_gt_int,'fro');
            RelErr_final(kkk,i_noise)     = MSE_Integer/norm(X_gt_int,'fro');
            fprintf("%s:　PSNR: %5.3f | Re: %5.3f |Time %5.2f \n",Method_list{kkk},PSNR_final(kkk,i_noise),RelErr_final(kkk,i_noise),Runtime_final(kkk,i_noise)); 
        end
    end
    %% RTD
    kkk=3;
    % Creating parallel pool
    if isempty(gcp('nocreate'))
        parpool; 
    end
    fprintf("===================%s==================\n",Method_list{kkk});
    %---parameters---%
    r=11;
    fprintf("-------------------Rank:%d---------------\n",r);
    for rho = 10.^[1] 
        RTD_seperate = zeros(size(X_gt_int));
        start_time = tic;
        parfor i_labnum = 1:label_len
            lable_id = label_list(i_labnum);
            fprintf("=");
            ind_lab = zeros(size(gt_lab_int));
            ind_lab(gt_lab_int==lable_id)=1;  
            [up,down,left,right] = measure_rtd_square(ind_lab,0);
            ind_lab_cube = repmat(ind_lab,[1,1,n3]);
            ind_lab_cube = ind_lab_cube(up:down,left:right,:);
            len_rtd = down-up+1;
            wid_rtd = right-left+1;


            X_gt = X_gt_int(up:down,left:right,:);
            gt_Ind = find(ind_lab_cube == 1);
            omega_lab_outside = find(ind_lab_cube == 0);
            gt_out_ind = omega_lab_outside;
            X_gt(omega_lab_outside)=0;

            X_ob =X_ob_int(up:down,left:right,:);
            X_ob(omega_lab_outside)=0;
            gt_omega = zeros(size(X_gt));
            gt_omega(gt_Ind) = 1;
            tol         = 1e-3;
            max_iter    = 300; 
            jjj = 0;
            temp_result_for_label = zeros(size(X_gt));
            jjj = jjj+1;
            [X_RTD, res2,iter_rtd,~,~,~] = RTD_proximal_dn(X_ob,rho,r,tol,max_iter,gt_Ind,gt_out_ind,X_gt,gt_omega);
            X_RTD(gt_out_ind) = 0; 
            temp_result_for_label=temp_result_for_label+ X_RTD;

            results{i_labnum}.data = temp_result_for_label;
            results{i_labnum}.coords = {up, down, left, right};
            
        end
        for i_labnum = 1:label_len
            res = results{i_labnum};
            up = res.coords{1}; down = res.coords{2};
            left = res.coords{3}; right = res.coords{4};
            RTD_seperate(up:down,left:right,:) = RTD_seperate(up:down,left:right,:) + res.data;
        end
        Runtime_final(kkk,i_noise) = toc(start_time);
        fprintf("\n");
        Result_final{kkk,i_noise} = RTD_seperate;
        [PSNR_final(kkk,i_noise)]= quality_std(X_gt_int*255,RTD_seperate*255);
        MSE_Integer    = norm(RTD_seperate-X_gt_int,'fro');
        RelErr_final(kkk,i_noise)     = MSE_Integer/norm(X_gt_int,'fro');
        fprintf("%s:　Rank: %d |　PSNR: %5.3f | Re: %5.3f |Time %5.2f \n",Method_list{kkk},r,PSNR_final(kkk,i_noise),RelErr_final(kkk,i_noise),Runtime_final(kkk,i_noise)); 
    end

end
for i_noise = 1
    
    fprintf("===================noise level %4.2f==================\n",nl);
    for i_me = 1:method_len
        fprintf("%8s:　PSNR: %5.3f |  Re: %5.3f | Time: %5.2f \n",Method_list{i_me},PSNR_final(i_me,i_noise),RelErr_final(i_me,i_noise),Runtime_final(i_me,i_noise));
    end
    
end
