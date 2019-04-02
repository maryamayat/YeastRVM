
% This function uses SparseBayes software (Version 2.0).
% Ensemble of linear basis RVM
% For identifying most effecting SNPs

function [RVs, CC_ensem] = GSr_yeastFS(kernel, qlabels, M, N)
% M is max number of RVMs in an ensemble.
% N is the size of subsample.

rng('shuffle');

% constructing a vector for summation of all test rounds (1 to M) 
%
num_indivs = length(qlabels);
y_sum_tests = zeros(num_indivs,1);


% creating RVs vector: saves the number of times that each seq has been
% chosen as an RV in M rounds (RVMs)
%
% RVs = zeros(length(original_train_ids),1);
RVs = zeros(N+1,1);

likelihood_='Gaussian';
iterations = 1000;
    
% ensemble loop
for i = 1:M
    i
    % creating a subsample for training (s) with size N 
    %
    % s = datasample(train_ids, N, 'Replace',false);
    s = datasample((1:num_indivs), 500, 'Replace',false);
    
    % creating basis_train and label_train
    %
    basis_train = kernel(s,:);
    outputs_train = qlabels(s);
    
    
    % train an RVM with the subsample s
    %
    %
    % Set up the options:
    % 
    % - we set the diagnostics level to 2 (reasonable)
    % - we will monitor the progress every 10 iterations
    % 
    OPTIONS		= SB2_UserOptions('iterations',iterations,...
							  'diagnosticLevel', 1,...
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
    w_infer	= zeros(N+1,1);
    w_infer(PARAMETER.Relevant)	= PARAMETER.Value;
    

    % test on the training data
    %
	%  y_train = basis_train*w_infer;

    % test on the test data!
    % 
	%  y_test = kernel*w_infer;
    y_test = kernel(setdiff(1:num_indivs,s),:)*w_infer;

    % add this y_test to the summamtion
    %
    y_sum_tests(setdiff(1:num_indivs,s)) = y_sum_tests(setdiff(1:num_indivs,s))+y_test;
    
    % Pearson's linear correlation coefficient on the train data
    %
	%  CC_train = corr(outputs_train,y_train);

    % Pearson's linear correlation coefficient on the test data
    %
    CC_test = corr(qlabels(setdiff(1:num_indivs,s)),y_test)
    
    % Updating RVs
    %
    % relevant_features = s(PARAMETER.Relevant);
    relevant_features = PARAMETER.Relevant;
    RVs(relevant_features) = RVs(relevant_features) + 1;
    
        
end % for loop (ensemble)
ensem_result = y_sum_tests/M;
CC_ensem = corr(qlabels,ensem_result);
                       
end



