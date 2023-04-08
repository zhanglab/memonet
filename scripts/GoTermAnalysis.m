

% GOtbl = readtable(uigetfile(), 'Sheet', 'x_vs_0_upreg');

use_cluster_or_tvc = 'cluster'; % cluster-based analysis
% use_cluster_or_tvc = 'tvc'; % train vs control



GOterm = 'inhibitory synapse assembly';

num_top_genes = 40;

switch use_cluster_or_tvc
    case 'cluster'
        term_address = cellfun(@(x) any(strcmpi(x, GOterm)), GOtbl{:,2});
        
        list = GOtbl{term_address,8};
        GO_genes = [];
        for clusters_repped = 1:length(list)
            GO_genes = [GO_genes, regexp(list{clusters_repped}, '/', 'split')];
        end
        GO_genes = unique(GO_genes);
        
        data_to_use = predictors;
        % Remove cluster 4, which was identified as likely contaminants
        data_to_use(:,clustertag==4) = nan;
        clusters_to_consider = [0 1 2 3 5];
        
        gene_address_list = []; GO_exp = []; 
        for g = 1:length(GO_genes)
            ga = logical(cellfun(@(x) strcmpi(x,GO_genes{g}),tbl{:,1}));
            gene_address_list = [gene_address_list; find(ga)];
            clust_count = 1;
            for c = clusters_to_consider
                GO_exp(g,clust_count) = nanmean(data_to_use(ga,find(clustertag==(c))));
                clust_count = clust_count+1;
            end
        end
        %=================================================
        %%% Sort genes by predictor scores
        gene_beta = nan(1,length(GO_genes));
        for g = 1:length(GO_genes)
            gene_beta(g) = LDAval(find(LDAind==gene_address_list(g)));
        end
        [~,abs_sorted_idx] = sort(abs(gene_beta),'descend');
        %
        if length(abs_sorted_idx)>num_top_genes
            abs_sorted_idx = abs_sorted_idx(1:num_top_genes);
        end
        GO_exp_for_im = GO_exp(abs_sorted_idx,:);
        GO_genes = GO_genes(abs_sorted_idx);
        
        %         GO_exp_train = GO_exp_train(sorted_index(1:30),:);
        %         GO_exp_ctrl = GO_exp_ctrl(sorted_index(1:30),:);
        %=================================================
        [goval,goind] = max(GO_exp_for_im,[],2);
        [goval2,goind2] =sort(goind);
        
        greenColorMap = [linspace(0,1,128), ones(1,128)];
        redColorMap = [ones(1, 128), linspace(1, 0, 128)];
        blueColorMap = [ones(1, 128), linspace(1, 0, 128)];
        colorMap = [redColorMap; greenColorMap; blueColorMap]';
        
        figure; t = tiledlayout('flow');
        nexttile([1,3])
        errorbar(0:length(clusters_to_consider)-1,nanmean(GO_exp,1), nanstd(GO_exp,[],1)./sqrt(sum(~isnan(GO_exp),1)),'k', 'linewidth', 2)
        hold on; plot(0:length(clusters_to_consider)-1, zeros(1,length(clusters_to_consider)),'--k')
        set(gca,'XTick', 0:length(clusters_to_consider)-1,'XTickLabel',cellfun(@num2str,mat2cell(clusters_to_consider',ones(1,length(clusters_to_consider)),1), 'uni', false))
        xlim([-0.5 4.5])
        nexttile([2 3])
        imagesc(GO_exp_for_im(goind2,:))
        set(gca, 'YTick', [1:size(GO_exp_for_im,1)], 'YTickLabel', GO_genes(goind2), 'XTick', 1:length(clusters_to_consider),'XTickLabel',cellfun(@num2str,mat2cell(clusters_to_consider',ones(1,length(clusters_to_consider)),1), 'uni', false) )
        title(GOterm)
        set(gca, 'clim', [-0.25 0.25],'colormap', colorMap)
        
        %================================================
        
    case 'tvc'
        term_address = cellfun(@(x) any(strcmpi(x, GOterm)), GOtbl_DGs{:,2});
        
        list = GOtbl_DGs{term_address,8};
        GO_genes = [];
   
        GO_genes = regexp(list{1}, '/', 'split');
        exclude_gene = 'fgf14';
        GO_genes = unique(GO_genes(~cellfun(@(x) strcmpi(x,exclude_gene),GO_genes)));
        
        data_to_use = predictors;
        % Remove cluster 4, which was identified as likely contaminants
        data_to_use(:,clustertag==4) = nan;
        
        gene_address_list = []; GO_exp_train = []; GO_exp_ctrl = [];
        for g = 1:length(GO_genes)
            ga = logical(cellfun(@(x) strcmpi(x,GO_genes{g}),tbl{:,1}));
            gene_address_list = [gene_address_list; find(ga)];
            GO_exp_train(g) = nanmean(data_to_use(ga,find(train_address)));
            GO_exp_ctrl(g) = nanmean(data_to_use(ga,find(~train_address)));
        end
        
        %=================================================
        % Sort genes by sign of top optimized LDA predictors
        gene_LDAidx = nan(1,length(GO_genes));
        for g = 1:length(GO_genes)
            gene_LDAidx(g) = find(LDAind==gene_address_list(g));
        end
        optclass_gene_address_list = gene_address_list(gene_LDAidx<=3000);
        optclass_gene_list = GO_genes(gene_LDAidx<=3000);
        
        [signed_optLDA_val,signed_optLDA_ind] = sort(optfinalbetas, 'descend');
        
        gene_beta = nan(1,length(optclass_gene_address_list));
        for g = 1:length(optclass_gene_address_list)
            gene_beta(g) = signed_optLDA_val(signed_optLDA_ind==gene_LDAidx(g));
        end
        [signed_sorted_val,signed_sorted_idx] = sort(gene_beta,'descend');
        
        pos_indxs = signed_sorted_idx(signed_sorted_val>0);
        pos_genes = optclass_gene_list(pos_indxs);
        if length(pos_indxs)>num_top_genes
            pos_indxs = pos_indxs(1:num_top_genes);
            pos_genes = pos_genes(1:num_top_genes);
        end
        neg_indxs = signed_sorted_idx(signed_sorted_val<0);
        neg_genes = optclass_gene_list(neg_indxs);
        if length(neg_indxs)>num_top_genes
            neg_indxs = neg_indxs(end-num_top_genes:end);
            neg_genes = neg_genes(end-num_top_genes:end);
        end
        
        GO_exp_ctrl_pos = GO_exp_ctrl(pos_indxs);
        GO_exp_ctrl_neg = GO_exp_ctrl(neg_indxs);
        GO_exp_train_pos = GO_exp_train(pos_indxs);
        GO_exp_train_neg = GO_exp_train(neg_indxs);
        
        pos_weighted_data = [GO_exp_ctrl_pos',GO_exp_train_pos'];
        neg_weighted_data = [GO_exp_ctrl_neg',GO_exp_train_neg'];
        
        % Magenta-black-green run
%         greenColorMap = [zeros(1, 128), linspace(0, 1, 128)];
%         redColorMap = [linspace(1, 0, 128), zeros(1, 128)];
%         blueColorMap = [linspace(1, 0, 128), zeros(1, 128)];
%         colorMap = [redColorMap; greenColorMap; blueColorMap]';
%         
        % Magenta-white-green run
        greenColorMap = [linspace(0,1,128), ones(1,128)];
        redColorMap = [ones(1, 128), linspace(1, 0, 128)];
        blueColorMap = [ones(1, 128), linspace(1, 0, 128)];
        colorMap = [redColorMap; greenColorMap; blueColorMap]';
        
        %Red-blue run
%         greenColorMap = [linspace(0, 1, 128), linspace(1,0,128)];
%         redColorMap = [ones(1, 128), linspace(1, 0, 128)];
%         blueColorMap = [linspace(0,1,128), ones(1,128)];
%         colorMap = [redColorMap; greenColorMap; blueColorMap]';


        figure('Name',[GOterm]); t = tiledlayout('flow');
        nexttile([2 1])
        [~,sortorder] = sort(diff(pos_weighted_data,[],2), 'descend'); % sort based on train data
        imagesc(pos_weighted_data(sortorder,:))
        set(gca, 'YTick', [1:length(pos_genes)], 'YTickLabel', pos_genes(sortorder), 'XTick', 1:2,'XTickLabel',{'Ctrl', 'Train'})
        set(gca, 'colormap', colorMap, 'clim',[-0.25 0.25]);
        title('Pos. Weight')
        nexttile([2 1])
        [~,sortorder] = sort(diff(neg_weighted_data,[],2), 'ascend'); % sort based on train data
        imagesc(neg_weighted_data(sortorder,:))
        set(gca, 'colormap', colorMap, 'clim',[-0.25 0.25]); 
        set(gca, 'YTick', [1:length(neg_genes)], 'YTickLabel', neg_genes(sortorder), 'XTick', 1:2,'XTickLabel',{'Ctrl', 'Train'})
        title('Neg. Weight')
end