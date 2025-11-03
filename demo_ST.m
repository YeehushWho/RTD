clear;
rng(5);
warning('off', 'all');

%%  loading data
data_selection=2;
if data_selection==1
    data_name = 'brain';
else
    data_name = 'kidney';
end

load_path  = 'E:\hyx\workspace\matlab\RTD\reply\demo_public_release\data_exp\sp_trans\';
load([load_path 'mouse_' data_name '.mat'])

Method_list = {'Observe','CPD-PAM','CPWOPT','RTD'}; 
method_len = length(Method_list);
sr_list = [0.5];
sr_length = length(sr_list);
R = [20];
tensor_name ={"Odd","Even"};
EN_CP =1;
EN_RTD = 1;
EN_CPWOPT=1;


% original tensor 
odd_gt = sptrans_tensor{1};
even_gt = sptrans_tensor{2};
odd_omega = sptrans_omega{1};
even_omega = sptrans_omega{2};
[n1,n2,n3] = size(odd_gt);
X_gt_int = sptrans_integ(odd_gt,even_gt);
gt_omega_int = sptrans_integ(odd_omega,even_omega);
gt_Ind_int = find(gt_omega_int==1);

% Saving results and metrics
Result_final = cell(method_len,sr_length); 
RelErr_final = zeros(method_len,sr_length);
Rsq_final = zeros(method_len,sr_length);
Runtime_final = zeros(method_len,sr_length);
RE_slice_final = zeros(method_len,sr_length,n3);
Rsq_slice_final = zeros(method_len,sr_length,n3);
%% 实验
% for sr = sr_list
for i_sr =[1]
    sr = sr_list(i_sr);
    fprintf("============= SR %3.1f============\n",sr);
    %% observed tensor
    kkk=1;
    fprintf("===================%s==================\n",Method_list{kkk});
    X_ob_odd  = X_ob_list{1,1};
    X_ob_even  = X_ob_list{1,2};
    [X_ob_int,Re_list,R_sq_list,RE_Integer,R_sq_Integer]= integer_sp(X_gt_int,X_ob_odd,X_ob_even,gt_omega_int,gene_num_list_k,gene_list);
    Result_final{kkk,i_sr} = X_ob_int;
    RelErr_final(kkk,i_sr) = RE_Integer;
    Rsq_final(kkk,i_sr)= R_sq_Integer;
    RE_slice_final(kkk,i_sr,:) = Re_list;
    Rsq_slice_final(kkk,i_sr,:) = R_sq_list;
    fprintf("%-8s |  Overall | RE: %5.3f | R-Square: %5.3f | Time: - \n",Method_list{kkk},RE_Integer,R_sq_Integer);
    %% CPD-PAM
    if EN_CP
        kkk=2;
        fprintf("===================%s==================\n",Method_list{kkk});
        cpd_pam_result = cell(2); 
        cpd_time = zeros(2);
        tol   = 1e-3;
        max_iter    = 200;
        r=R;
        rho = 10.^[3];
        for i_te = [1,2] %　1 odd 2 even
            % create observation
            X_gt =sptrans_tensor{i_te};
            gt_omega = sptrans_omega{i_te};
            gt_out_ind = sptrans_outind{i_te};
            gt_Ind = sptrans_ind{i_te}; 
            ob_Omega = sptrans_ob_omega{i_te};
            ob_ind = find(ob_Omega==1);  
            X_ob = X_ob_list{1,i_te};
            %---main iteration---%
            t_Start = tic;
            [X_CP, res2,iter_rtd,~,~,~] = CP_proximal(X_ob,rho,r,tol,max_iter,gt_Ind,ob_ind,X_gt,gt_omega);
            cpd_time(i_te) = toc(t_Start);
            cpd_pam_result{i_te} =X_CP ;
        end
        [CPD_int,Re_list,R_sq_list,RE_Integer,R_sq_Integer]= integer_sp(X_gt_int,cpd_pam_result{1},cpd_pam_result{2},gt_omega_int,gene_num_list_k,gene_list);
        Result_final{kkk,i_sr} = CPD_int;
        RelErr_final(kkk,i_sr) = RE_Integer;
        Rsq_final(kkk,i_sr)= R_sq_Integer;
        RE_slice_final(kkk,i_sr,:) = Re_list;
        Rsq_slice_final(kkk,i_sr,:) = R_sq_list;
        Runtime_final(kkk,i_sr) = cpd_time(1)+cpd_time(2);
        fprintf("%-8s |  Overall | RE: %5.3f | R-Square: %5.3f | Time: %4.1f \n",Method_list{kkk},RE_Integer,R_sq_Integer,Runtime_final(kkk,i_sr));
    end

    %% RTD 
    kkk=4; 
    fprintf("===================%s==================\n",Method_list{kkk});
    rtd_pam_result = cell(2); 
    rtd_pam_time = zeros(2);
    %---parameters---%
    tol         = 1e-3;
    max_iter    = 200;
    rho = 10.^[2];
    r=R;
    for i_te = [1,2] %　1 odd 2 even
        % create observation
        X_gt =sptrans_tensor{i_te};
        gt_omega = sptrans_omega{i_te};
        gt_out_ind = sptrans_outind{i_te};
        gt_Ind = sptrans_ind{i_te}; 
        ob_Omega = sptrans_ob_omega{i_te};
        ob_ind = find(ob_Omega==1);  
        X_ob = X_ob_list{1,i_te};
        %---main iteration---%
        t_Start = tic;
        [X_RTD, res2,iter_rtd,~,~,~]  = RTD_proximal_gene(X_ob,rho,r,tol,max_iter,gt_Ind,gt_out_ind,ob_ind,X_gt,gt_omega,gene_list);
        rtd_pam_time(i_te) = toc(t_Start);
        X_RTD(ob_ind) = X_gt(ob_ind);
        X_RTD(gt_out_ind) = 0; 
        rtd_pam_result{i_te} =X_RTD ;
    end
    [RTD_int,Re_list,R_sq_list,RE_Integer,R_sq_Integer]= integer_sp(X_gt_int,rtd_pam_result{1},rtd_pam_result{2},gt_omega_int,gene_num_list_k,gene_list);
    Runtime_final(kkk,i_sr) = rtd_pam_time(1)+rtd_pam_time(2);
    RelErr_final(kkk,i_sr) = RE_Integer;
    Rsq_final(kkk,i_sr)= R_sq_Integer;
    RE_slice_final(kkk,i_sr,:) = Re_list;
    Rsq_slice_final(kkk,i_sr,:) = R_sq_list;
    fprintf("%-8s |  Overall | RE: %5.4f | R-Square: %5.3f | Time: %4.1f \n",Method_list{kkk},RE_Integer,R_sq_Integer,Runtime_final(kkk,i_sr));
    %% CP-WOPT 
    % Due to the prolonged time cost of CP-WOPT, its results have been saved in this project.

    kkk=3;
    fprintf("===================%s==================\n",Method_list{kkk});
    load(['./result/sp_trans/cpwopt_sptrans_' data_name '.mat']); % Considering the prolong runtime of CP-WOPT on datasets of ST, we also saved and provided its results.
    [WOPT_int,Re_list,R_sq_list,RE_Integer,R_sq_Integer]= integer_sp(X_gt_int,cp_wopt_result{1},cp_wopt_result{2},gt_omega_int,gene_num_list_k,gene_list);
    Runtime_final(kkk,i_sr) = cp_wopt_time(1)+cp_wopt_time(2);
    RelErr_final(kkk,i_sr) = RE_Integer;
    Rsq_final(kkk,i_sr)= R_sq_Integer;
    RE_slice_final(kkk,i_sr,:) = Re_list;
    Rsq_slice_final(kkk,i_sr,:) = R_sq_list;
    fprintf("%-8s |  Overall | RE: %5.3f | R-Square: %5.3f | Time: %4.1f \n",Method_list{kkk},RE_Integer,R_sq_Integer,Runtime_final(kkk,i_sr));
    
    % To retest this method, please add comments to the specified code section and execute the main code block below.

    % cp_wopt_result = cell(2); 
    % cp_wopt_time = zeros(2);
    % %---parameters---%
    % ncg_opts = ncg('defaults');
    % ncg_opts.StopTol = 1e-6;% 
    % ncg_opts.RelFuncTol = 1.0e-20; 
    % ncg_opts.MaxIters = 10^4; %10^4
    % ncg_opts.DisplayIters = 2000;
    % r=R;
    % for i_te = [1,2] %　1 odd 2 even
    %     X_gt =sptrans_tensor{i_te};
    %     gt_omega = sptrans_omega{i_te};
    %     gt_out_ind = sptrans_outind{i_te};
    %     gt_Ind = sptrans_ind{i_te}; 
    %     ob_Omega = sptrans_ob_omega{i_te};
    %     ob_ind = find(ob_Omega==1);  
    %     X_ob = X_ob_list{1,i_te};
    %     t_Start = tic;
    %     [X_cpwopt,~,output] = cp_wopt(X_ob, ob_Omega, r, 'opt', 'ncg', 'opt_options', ncg_opts); 
    %     X_cpwopt = full(X_cpwopt);
    %     cp_wopt_time(i_te) = toc(t_Start);
    %     X_cpwopt(ob_ind) = X_gt(ob_ind);
    %     X_cpwopt(gt_out_ind) = 0;
    %     X_cpwopt = double(X_cpwopt);
    %     cp_wopt_result{i_te} = X_cpwopt;
    % end
    % [WOPT_int,Re_list,R_sq_list,RE_Integer,R_sq_Integer]= integer_sp(X_gt_int,cp_wopt_result{1},cp_wopt_result{2},gt_omega_int,gene_num_list_k,gene_list);
    % Runtime_final(kkk,i_sr) = cp_wopt_time(1)+cp_wopt_time(2);
    % RelErr_final(kkk,i_sr) = RE_Integer;
    % Rsq_final(kkk,i_sr)= R_sq_Integer;
    % RE_slice_final(kkk,i_sr,:) = Re_list;
    % Rsq_slice_final(kkk,i_sr,:) = R_sq_list;
    % fprintf("%-8s |  Overall | RE: %5.3f | R-Square: %5.3f | Time: %4.1f \n",Method_list{kkk},RE_Integer,R_sq_Integer,Runtime_final(kkk,i_sr));
end
