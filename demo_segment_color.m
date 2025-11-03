clear;
rng(5);
warning('off', 'all');



data_selection = 1; % 1. chair |　2. sofa
if data_selection==1
    data_name = 'chair';
    r = 10;
else
    data_name = 'sofa';
    r=6;
end

load_path = './data_exp/segmantic/';
load([load_path 'color_segments_' data_name '.mat']); 
Method_list = {'Observed','CPD-PAM','RTD','CP-WOPT'};
method_len  = length(Method_list);
EN_CP =1;
EN_CPWOPT =1;
EN_RTD = 1;


sr_list = [0.4,0.5];
sr_length = length(sr_list);
Result_final = cell(method_len,sr_length);
Result_std_list_final = cell(method_len,sr_length);
PSNR_final  = zeros(method_len,sr_length);
RelErr_final = zeros(method_len,sr_length);
Runtime_final =  zeros(method_len,sr_length);
rank_list = [10,8];
for i_sr = [1:1:sr_length]
    sr = sr_list(i_sr);
    fprintf("===================%4.2f==================\n",sr);
    ob_Omega = ob_omega_list{i_sr,1};
    ob_ind = find(ob_Omega);
       
    ob_ind_std = ob_ind_out_list{i_sr};

    ob_Omega_std = zeros(size(X_gt));
    ob_Omega_std(ob_ind_std)=1;
    ob_Omega_std(ob_ind) = 1;
    ob_ind_int = find(ob_Omega_std==1);

    Ind = gt_Ind;
    X_ob = zeros(size(X_gt));
    X_ob(ob_ind) = X_gt(ob_ind);
    X_ob_std =  zeros(size(X_gt_std));
    X_ob_std(ob_ind_int) =  X_gt_std(ob_ind_int);
    r = rank_list(data_selection);
    %% observed tensor
    kkk=1;
    fprintf("===================%s==================\n",Method_list{kkk});
    Result_final{kkk,i_sr} = X_ob;
    [PSNR_final(kkk,i_sr),RelErr_final(kkk,i_sr)] = quality_rtd(X_gt,X_ob,gt_Ind);
    fprintf("%s | SR: %4.2f |PSNR: %5.3f  | RE: %6.3f \n",Method_list{kkk},sr_list(i_sr),PSNR_final(kkk,i_sr),RelErr_final(kkk,i_sr));
    %% CPD-PAM

    kkk=2;
    if EN_CP
        fprintf("===================%s==================\n",Method_list{kkk});
        tol   = 1e-3;
        max_iter    = 2000;
        jjj = 0;
 
        fprintf("-------------------Rank:%d---------------\n",r);
        for rho = 10.^[3] 
            jjj = jjj+1;
            start_time = tic;
            [X_CP, res2,iter_rtd,~,~,~] =CP_proximal(X_ob,rho,r,tol,max_iter,gt_Ind,ob_ind,X_gt);
            Runtime_final(kkk,i_sr) = toc(start_time);
            X_CP(ob_ind_int) =  X_gt_std(ob_ind_int);
            Result_std_list_final{kkk,i_sr} = X_CP;
            X_CP_rtd = X_CP;
            X_CP_rtd(gt_out_ind) = 0; 
            Result_final{kkk,i_sr} = X_CP_rtd;
            [PSNR_final(kkk,i_sr),RelErr_final(kkk,i_sr)]= quality_rtd(X_gt,X_CP_rtd,gt_Ind);
            fprintf("%s:　PSNR: %5.3f |  Re: %5.3f | Time: %5.2f \n",Method_list{kkk},PSNR_final(kkk,i_sr),RelErr_final(kkk,i_sr),Runtime_final(kkk,i_sr));
        end
    end
    %% RTD
    kkk=3;
    fprintf("===================%s==================\n",Method_list{kkk});
    %---parameters---%
    tol         = 1e-3;
    max_iter    = 300; 
    jjj = 0;

    fprintf("-------------------Rank:%d---------------\n",r);
    for rho = 10.^[1] 
        fprintf("-----------------rho:%f----------------------\n",rho);
        jjj = jjj+1;
        %---main iteration---%
        t_Start = tic;
        [X_RTD, res2,iter_rtd,~,~,~] = RTD_proximal(X_ob,rho,r,tol,max_iter,Ind,gt_out_ind,ob_ind,X_gt);
        runtimejjj = toc(t_Start);
        X_RTD(ob_ind) = X_gt(ob_ind);
        X_RTD(gt_out_ind) = 0; 
        Result_final{kkk,i_sr} = X_RTD;
        Runtime_final(kkk,i_sr) = runtimejjj;
        [PSNR_final(kkk,i_sr),RelErr_final(kkk,i_sr)]= quality_rtd(X_gt,Result_final{kkk,i_sr},gt_Ind);
        fprintf("%s:　PSNR: %5.3f |  Re: %5.3f | Time: %5.2f \n",Method_list{kkk},PSNR_final(kkk,i_sr),RelErr_final(kkk,i_sr),Runtime_final(kkk,i_sr));
    end

    %% CP-WOPT
    if EN_CPWOPT
        kkk=4;

        ncg_opts = ncg('defaults');
        ncg_opts.StopTol = 1.0e-6;
        ncg_opts.RelFuncTol = 1.0e-20;
        ncg_opts.MaxIters = 10^4;
        ncg_opts.DisplayIters = 2000;
        t_Start = tic;
        [X_cpwopt,~,output] = cp_wopt(X_ob, ob_Omega, r,'opt', 'ncg', 'opt_options', ncg_opts);  
        X_cpwopt = full(X_cpwopt);
        Runtime_final(kkk,i_sr) = toc(t_Start);
        X_cpwopt(ob_ind) = X_gt(ob_ind);
        X_cpwopt(gt_out_ind) = 0;
        Result_final{kkk,i_sr} = X_cpwopt;
        [PSNR_final(kkk,i_sr),RelErr_final(kkk,i_sr)]= quality_rtd(X_gt,Result_final{kkk,i_sr},gt_Ind);
        fprintf("%s:　PSNR: %5.3f |  Re: %5.3f | Time: %5.2f \n",Method_list{kkk},PSNR_final(kkk,i_sr),RelErr_final(kkk,i_sr),Runtime_final(kkk,i_sr));

    end
end
for i_sr = [1:1:sr_length]
    sr = sr_list(i_sr);
    fprintf("===================%4.2f==================\n",sr);
    for i_me = 1:method_len
        fprintf("%8s:　PSNR: %5.3f |  Re: %5.3f | Time: %5.2f \n",Method_list{i_me},PSNR_final(i_me,i_sr),RelErr_final(i_me,i_sr),Runtime_final(i_me,i_sr));
    end
    
end
%% close all;
