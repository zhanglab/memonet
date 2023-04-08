Lambda = [0,logspace(-5,1,20)]; % List of lambdas



for shuff = 1:100
    shuffresponsevar = shake(responsevar);
    Mdl_CV = cell(kfold,1);
    Beta_CV = nan(size(predictors,1),kfold);
    Bias_CV = nan(size(predictors,1),kfold);
    temp_predictions = cell(kfold,1);
    test_data = cell(kfold,1);
    for cv_ind = 1:kfold
        tim = tic;
        [~,temp_minMESind{cv_ind}] = min(kfoldLoss(fitclinear(predictors(:,train_test_ind ~= cv_ind),shuffresponsevar(train_test_ind ~= cv_ind),'ObservationsIn','columns','KFold',10,'Lambda',Lambda,'Learner',mdl.Learner,'Regularization','lasso')));
        temp_model{cv_ind} = selectModels(fitclinear(predictors(:,train_test_ind ~= cv_ind),shuffresponsevar(train_test_ind ~= cv_ind),'ObservationsIn','columns','Lambda',Lambda,'Learner',mdl.Learner,'Regularization','lasso'),temp_minMESind{cv_ind});
        Mdl_CV{cv_ind} = struct(temp_model{cv_ind});
        Beta_CV(:,cv_ind) = Mdl_CV{cv_ind}.Beta;
        Bias_CV(cv_ind) = Mdl_CV{cv_ind}.Bias;
        temp_predictions{cv_ind} = predict(temp_model{cv_ind},predictors(:,train_test_ind == cv_ind)');
        test_data_CV{cv_ind} = shuffresponsevar(train_test_ind == cv_ind);
        toc(tim)
    end
    temp_predictions = cell2mat(temp_predictions);
    test_data = cell2mat(test_data_CV');
    % Decode dQ from test set
    prediction_cellmatrix = cell(kfold,1);
    for cv_ind = 1:kfold
        prediction_cellmatrix{cv_ind} = predict(temp_model{cv_ind},predictors(:,train_test_ind == cv_ind)');
    end
    prediction_matrix = cell2mat(prediction_cellmatrix);
    
    conmat = confusionmat(test_data, logical(prediction_matrix));

    ShuffConMat(:,:,shuff) = conmat;
    
    ShuffPredictionAccuracy(shuff) = (conmat(1,1)+conmat(2,2))/sum(conmat(:));

    disp(['Prediction accuracy = ', num2str(ShuffPredictionAccuracy(shuff))])
    
    % mdlBetas = cellfun(@(x) x.Beta, mdl2.Trained, 'uni', false);
    shuffLDAbetas = nanmedian(Beta_CV,2);
    % tgonlyLDAbetas = nanmean(Beta_CV,2);
    
    [shuffval,shuffind] = sort(abs(shuffLDAbetas), 'descend');
    
    ShuffBetas(:,:,shuff) = [shuffval,shuffind];
end