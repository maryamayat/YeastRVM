
% This function uses SparseBayes software (Version 2.0).
% Predicting Quantitative Traits

function [RVs] = GSr_wheat_RVM(kernel, qlabels, n_folds, kernel_name, trait_no, output_fname)
% results file: output_fname
% headers: trait_no, try_no, kernel_name, nfolds, CC_test  
%
%
% n_folds is the number of folds in CV.

n_indvs = length(qlabels); %number of individuals
rng('default');
for j = 1:10 % 10 times repeat CV
    CVO = cvpartition(n_indvs,'KFold',n_folds);


    % creating RVs vector: saves the number of times that each seq has been
    % chosen as an RV in M rounds (RVMs)
    %
    RVs = zeros(n_indvs,1);

    % creating weigths vector: saves the number of sum of weights for seqs which are RVs 
    % in M rounds (RVMs)
    %
    % weights = zeros(length(original_train_ids),1);

    likelihood_='Gaussian';
    iterations = 1000;

    % CV loop
    CV_CC_train = zeros(n_folds,1);
    CV_CC_test = zeros(n_folds,1);
    y_test_allfolds = zeros(n_indvs,1);
    for i = 1:n_folds
        train_ids = CVO.training(i);
        test_ids = CVO.test(i);

        % creating basis_train and label_train
        %
        basis_train = kernel(train_ids,train_ids);
        outputs_train = qlabels(train_ids);

        % creating basis_test and label_test
        %
        basis_test = kernel(test_ids,train_ids);
        outputs_test = qlabels(test_ids);

        % train an RVM 
        %
        %
        % Set up the options:
        % 
        % - we set the diagnostics level to 2 (reasonable)
        % - we will monitor the progress every 10 iterations
        % 
        OPTIONS		= SB2_UserOptions('iterations',iterations,...
                                  'diagnosticLevel', 2,...
                                  'monitor', 10);
        % Set initial parameter values:
        % SPARSEBAYES will call SB2_PARAMETERSETTINGS itself to obtain an 
        % appropriate default for the noise (and other SETTINGS fields).
        % 
        SETTINGS	= SB2_ParameterSettings('NoiseStd',0.1);

        % Now run the main SPARSEBAYES function
        %
        [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = ...
            SparseBayes(likelihood_, basis_train, outputs_train, OPTIONS, SETTINGS);

        % Manipulate the returned weights
        %
        w_infer	= zeros(CVO.TrainSize(i),1);
        w_infer(PARAMETER.Relevant)	= PARAMETER.Value;


        % test on the training data
        %
        y_train = basis_train*w_infer;

        % test on the test data!
        % 
        y_test = basis_test*w_infer;

        % save this fold y_test for calculating PCC on all data at the end
        %
        y_test_allfolds(test_ids) = y_test;


        % Pearson's linear correlation coefficient on the train data
        %
        CV_CC_train(i) = corr(outputs_train,y_train);

        % Pearson's linear correlation coefficient on the test data
        %
        CV_CC_test(i) = corr(outputs_test,y_test);

        % Updating RVs
        %
        s = find(train_ids==1);
        relevant_seqs = s(PARAMETER.Relevant);
        RVs(relevant_seqs) = RVs(relevant_seqs) + 1;

        % Print
        DIAGNOSTIC.iterations

    end % for loop (i:CV)

    CV_CC_train
    CV_CC_test
    %write the result to the output_fname file
    %
    CC_test = corr(qlabels,y_test_allfolds)
    result_vector = [trait_no j kernel_name n_folds CC_test]; 
    dlmwrite(output_fname, result_vector,'-append');
end % for loop (j:repeat CV)
                       
end



