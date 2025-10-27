clear all;
close all;
addpath('dataset'); 
addpath('Commoncodes');

CluAcc = zeros(6, 1);
accClu = zeros(6, 6);

ALPHA    = 1;
BETA     = 2;
MAX      = 10;
LAMBDA_1 = 1;


allc = [10];
for i = 1 : length(allc)
    nCluster = allc(i);

    all_time_costs = []; 

    lambda1all =  [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3];
    for iparas1 = 1  : length(lambda1all)
        lambda1 = lambda1all(iparas1);
        disp([' lambda1 = ' num2str(lambda1)]);


        lambda2all =  [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3];
        for iparas2 =  1    : length(lambda2all)
            lambda2 = lambda2all(iparas2);


            load yaleb  % YaleB dataset
            num = nCluster*64;  % number of data used for subspace segmentation
            start = 0 ;
            fea = fea(:,start+1:start+num) ;
            gnd = gnd(:,start+1:start+num) ;

            %% Projection
            % PCA
            [ eigvector , eigvalue ] = PCA( fea ) ;
            maxDim = length(eigvalue);
            fea = eigvector' * fea ;
            % fea = [fea ; ones(1,size(fea,2)) ] ;
            for i = 1 : num
                fea(:,i) = fea(:,i) / norm(fea(:,i)) ;
            end
            d = nCluster*6;
            data = fea(1:d,:) ;


            para.knn = 4;
            para.gamma =7;
            para.elpson = 0.001;
            para.aff_type = 'J2'; 
            p = 1.0;
            rc = 0.01;
            rho = 1.1;
            Q = orth(data');
            A = data*Q;  nX = sqrt(sum(data.^2));
            tic;
            [W, value] = ABN_FRMD(data, A, p,lambda1,lambda2,rc,rho,para);
            time_cost = toc;

            all_time_costs(end+1) = time_cost;

            W = Q*W; J = W;
            if strcmp(para.aff_type,'J1')
                L =(abs(J)+abs(J'))/2;
            elseif strcmp(para.aff_type,'J2')
                L=abs(J'*J./(nX'*nX)).^para.gamma;
            elseif strcmp(para.aff_type,'J2_nonorm')
                L=abs(J'*J).^para.gamma;
            end
            W = L;


%             [n1, n2] = size(data); 
%             opts.lambdal1 = 8 * sqrt(max(n1, n2)); 
%             opts.rank = 24; 
%             opts.p = 0.5; 
%             opts.numIter = 150;
%             opts.tol = 1e-3; 
%             [U, V] = RLMF(data, ones(size(data,1),1), opts); 
% 
%             data_new = data * V';
            
    
            W_init = fea(1:d,:);
            H_init = W;
            tic;
            [~, H, info] = h_snmf_l1(data, ALPHA, BETA, LAMBDA_1, W_init, H_init);
            time_cost1 = toc;


            W2 = H;

            %%
            for ic = 1 : size(H,2)
                W2(:,ic) = H(:,ic)/(max(abs(H(:,ic)))+eps) ;
            end
  
            groups = clu_ncut(W2,max(gnd));
            [ACC, NMI, PUR] = ClusteringMeasure(gnd,groups);

            disp([' lambda2 = ' num2str(lambda2), ' choosecluster = ' num2str(nCluster), ' acc = ' num2str(ACC*100), ' nmi = ' num2str(NMI*100), ' PUR= ' num2str(PUR*100), ' time= ' num2str(time_cost)]);
            
            
            CluAcc(iparas2) = ACC*100;

        end

        accClu(:,iparas1)  =  CluAcc;
    end
       

    avg_time = mean(all_time_costs);
    fprintf('Average time cost for nCluster = %d: %.4f seconds\n', nCluster, avg_time);

    eval(['accLGRSmLRRYaleB_' num2str(nCluster) '= accClu']);
    eval(['save accLGRSmLRRYaleB_' num2str(nCluster) ' accLGRSmLRRYaleB_' num2str(nCluster)]);

end
