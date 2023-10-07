% this code calculates the prediction accuracy of the classifier by animal
% Code 1 is identical to that used to optimize the top gene number, but has additions to store predictions by animal
% Code 2 summarizes the values



% Code 1:

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
    opt_predictors = predictors(LRind(1:t),:);
    optMdl_CV = cell(kfold,1);
    optBeta_CV = nan(size(opt_predictors));
    optBias_CV = nan(size(opt_predictors));
    opttemp_predictions = cell(kfold,1);
    opttemp_animal_tag = cell(kfold,1);
    opttest_data = cell(kfold,1);
    for cv_ind = 1:kfold
        tim = tic;
        [~,temp_minMESind{cv_ind}] = min(kfoldLoss(fitclinear(opt_predictors(:,train_test_ind ~= cv_ind),responsevar(train_test_ind ~= cv_ind),'ObservationsIn','columns','KFold',10,'Lambda',Lambda,'Learner',mdl.Learner,'Regularization','lasso')));
        opttemp_model{cv_ind} = selectModels(fitclinear(opt_predictors(:,train_test_ind ~= cv_ind),responsevar(train_test_ind ~= cv_ind),'ObservationsIn','columns','Lambda',Lambda,'Learner',mdl.Learner,'Regularization','lasso'),temp_minMESind{cv_ind});
        optMdl_CV{cv_ind} = struct(opttemp_model{cv_ind});
        optBeta_CV(:,cv_ind) = optMdl_CV{cv_ind}.Beta;
        optBias_CV(cv_ind) = optMdl_CV{cv_ind}.Bias;
        opttemp_predictions{cv_ind} = predict(opttemp_model{cv_ind},opt_predictors(:,train_test_ind == cv_ind)');
        opttemp_animal_tag{cv_ind} = animaltag(train_test_ind == cv_ind);
        opttest_data{cv_ind} = responsevar(train_test_ind == cv_ind);
        toc(tim)
        disp(['Done with ', num2str(cv_ind), '/10 CV at t = ', num2str(t)])
    end
    opttemp_predictions_mat = cell2mat(opttemp_predictions);
    opttest_data_mat = cell2mat(opttest_data);
    % Decode dQ from test set
    prediction_matrix = cell(kfold,1);
    for cv_ind = 1:kfold
        prediction_matrix{cv_ind} = [prediction_matrix{cv_ind}, predict(opttemp_model{cv_ind},opt_predictors(:,train_test_ind == cv_ind)')];
        opt_cv_PA{cv_ind} = confusionmat(responsevar(train_test_ind == cv_ind),predict(opttemp_model{cv_ind},opt_predictors(:,train_test_ind == cv_ind)'));
    end
    prediction_matrix = cell2mat(prediction_matrix);
%     % Calculate decoding accuracy as Pearson correlation coeff
%     [DecodingAccuracy,corr_p_val] = corr(prediction_matrix,opttest_data,'Type','Pearson','Rows','complete');
%
    optconmat = confusionmat(opttest_data, logical(prediction_matrix));
    disp(['Prediction accuracy = ', num2str((optconmat(1,1)+optconmat(2,2))/sum(optconmat(:)))])
%     opt_pred_acc(t_count) = (optconmat(1,1)+optconmat(2,2))/sum(optconmat(:));
%     t_count = t_count+1;
    optfinalbetas = nanmean(optBeta_CV,2);
    [optLDAval,optLDAind] = sort(abs(optfinalbetas), 'descend');
    %================================================================
    % Test prediction accuracy per animal
    opttemp_animal_tag_mat = cell2mat(opttemp_animal_tag);
    for a = 1:6
        AnimalPredictions{a} = opttemp_predictions_mat(opttemp_animal_tag_mat==a);
        AnimalLabels{a} = opttest_data_mat(opttemp_animal_tag_mat==a);
        AnimalConfMat{a} = confusionmat(AnimalPredictions{a}, AnimalLabels{a});
    end
    for cv_ind = 1:kfold
        for a = 1:6
            preds = opttemp_predictions{cv_ind}(opttemp_animal_tag{cv_ind}==a);
            labes = opttest_data{cv_ind}(opttemp_animal_tag{cv_ind}==a);
            confmat = confusionmat(preds, labes);
            AnimalPredictionArray{a}(1,cv_ind) = sum(diag(confmat))./sum(confmat(:));
        end
    end
end



% Code 2:

AnimalPredictionArray = cell(1,max(animaltag));
for cv_ind = 1:kfold
for a = 1:6
preds = opttemp_predictions{cv_ind}(opttemp_animal_tag{cv_ind}==a);
labes = opttest_data{cv_ind}(opttemp_animal_tag{cv_ind}==a);
confmat = confusionmat(preds, labes);
AnimalPredictionArray{a}(1,cv_ind) = sum(diag(confmat))./sum(confmat(:));
end
end