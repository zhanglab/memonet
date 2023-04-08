use_cells = 1:length(animaltag);
predictors = zscore(tbl{:,use_cells+1},[],2);
responsevar = train_address(use_cells);


numCells = size(predictors,2);
kfold = 10;
train_test_ind = crossvalind('Kfold', responsevar, 10);

%%%% Cross-validation for prediction accuracy

Lambda = [0,logspace(-5,1,20)]; % List of lambdas
t_count= 1;

for t = 3000
    opt_predictors = predictors(LDAind(1:t),:);
    Mdl_CV = cell(kfold,1);
    Beta_CV = nan(size(opt_predictors));
    Bias_CV = nan(size(opt_predictors));
    temp_predictions = cell(kfold,1);
    test_data = cell(kfold,1);
    for cv_ind = 1:kfold
        tim = tic;
        [~,temp_minMESind{cv_ind}] = min(kfoldLoss(fitclinear(opt_predictors(:,train_test_ind ~= cv_ind),responsevar(train_test_ind ~= cv_ind),'ObservationsIn','columns','KFold',10,'Lambda',Lambda,'Learner',mdl.Learner,'Regularization','lasso')));
        temp_model{cv_ind} = selectModels(fitclinear(opt_predictors(:,train_test_ind ~= cv_ind),responsevar(train_test_ind ~= cv_ind),'ObservationsIn','columns','Lambda',Lambda,'Learner',mdl.Learner,'Regularization','lasso'),temp_minMESind{cv_ind});
        Mdl_CV{cv_ind} = struct(temp_model{cv_ind});
        Beta_CV(:,cv_ind) = Mdl_CV{cv_ind}.Beta;
        Bias_CV(cv_ind) = Mdl_CV{cv_ind}.Bias;
        temp_predictions{cv_ind} = predict(temp_model{cv_ind},opt_predictors(:,train_test_ind == cv_ind)');
        test_data{cv_ind} = responsevar(train_test_ind == cv_ind);
        toc(tim)
        disp(['Done with ', num2str(cv_ind), '/10 CV at t = ', num2str(t)])
    end
    temp_predictions = cell2mat(temp_predictions);
    test_data = cell2mat(test_data);
    % Decode dQ from test set
    prediction_matrix = cell(kfold,1);
    for cv_ind = 1:kfold
        prediction_matrix{cv_ind} = [prediction_matrix{cv_ind}, predict(temp_model{cv_ind},opt_predictors(:,train_test_ind == cv_ind)')];
        opt_cv_PA{cv_ind} = confusionmat(responsevar(train_test_ind == cv_ind),predict(temp_model{cv_ind},opt_predictors(:,train_test_ind == cv_ind)'));
    end
    prediction_matrix = cell2mat(prediction_matrix);
    % Calculate decoding accuracy as Pearson correlation coeff
    [DecodingAccuracy,corr_p_val] = corr(prediction_matrix,test_data,'Type','Pearson','Rows','complete');
    
    optconmat = confusionmat(test_data, logical(prediction_matrix));
    
    disp(['Prediction accuracy = ', num2str((optconmat(1,1)+optconmat(2,2))/sum(optconmat(:)))])
%     opt_pred_acc(t_count) = (optconmat(1,1)+optconmat(2,2))/sum(optconmat(:));
%     t_count = t_count+1;
    optfinalbetas = nanmean(Beta_CV,2);
    [optLDAval,optLDAind] = sort(abs(optfinalbetas), 'descend');
end