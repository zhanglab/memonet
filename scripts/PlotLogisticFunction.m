betas_to_use = optfinalbetas;
bias_to_use = nanmean(Bias_CV(1:10));

genes_to_use = LDAind(1:3000);

data_to_use = predictors(genes_to_use,:);

linear_model_x = bias_to_use + betas_to_use' * data_to_use;

logit_model_x = 1./(1+exp(-linear_model_x));

figure; plot(linear_model_x(~train_address), logit_model_x(~train_address), 'ok', 'markerfacecolor', 'k')
hold on; plot(linear_model_x(train_address), logit_model_x(train_address), 'ok', 'markerfacecolor', 'g')