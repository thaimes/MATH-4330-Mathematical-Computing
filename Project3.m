% Thomas Haimes
% MATH 4330
% Project 3 - Machine Learning with Support Vector Machines

% Clear console for neatness
clc;
clear;
close all;

% Update this with the PATH for your dataset
filepath = "C:\Users\thoma\Documents\MATLAB\Mathematical Computing\heart_disease\processed.cleveland.data";

% Need to modify this to make any dataset work... 
A = table2array(readtable(filepath, FileType="text"));

[m,n] = size(A);


%% Clean up the data

% Remove all missing data
for i=1:n       
    I = isnan(A(:,i));
    J = find(I == 1);
    A(J,:) = [];
end

% Shift and scale the data
for i = 1:n-1   
    meanvalue = mean(A(:, i));
    A(:, i) = A(:, i) - meanvalue;
    maxvalue = max(abs(A(:, i)));
    A(:, i) = A(:, i)/maxvalue;
end

%% Set up the training and testing sets

k = ceil(m*0.2);
maxiter = 20;
[m,n] = size(A);

% Initialize the accuracy
accuracy_svm =zeros(maxiter, 1);
accuracy_nn = zeros(maxiter, 1);

for iter = 1:maxiter
    I = randi(m, k, 1);
    I = unique(I); % the test set indices
    fulldata = 1:m;

    J = transpose(setdiff(fulldata, I)); % the training set indices
    trainingdata = A(J, :);
    testdata = A(I, :);

    % Set up the SVM, neural net model
    mdl_svm = fitcecoc(trainingdata(:, 1:n-1), trainingdata(:, n)); % fitcsvm -> fitcecoc

    % Compute SVM model parameters
    pred_svm = predict(mdl_svm, testdata(:, 1:n-1)); % Apply the model to the test data

    mdl_nn = fitcnet(trainingdata(:, 1:n-1), trainingdata(:, n)); 

    % Compute neural net model parameters
    pred_nn = predict(mdl_nn, testdata(:, 1:n-1)); % Apply model to the test data

    % Find the accuracy
    errors_svm = length(find(pred_svm ~= testdata(:, n)));
    % Find the error count in the prediction of the SVM model for the
    % 'iter' training/test set. ^^^^

    error_rate_svm = errors_svm/length(I)*100;
    accuracy_svm(iter, 1) = 100 - error_rate_svm;

    errors_nn = length(find(pred_nn ~= testdata(:, n)));
    % Find the error count in the prediction of the neural net model for
    % 'iter' training/test set. ^^^^

    error_rate_nn = errors_nn/length(I)*100;
    accuracy_nn(iter, 1) = 100 - error_rate_nn;
end

% Collect together the performance of SVM for test data sets.
figure(1)
boxplot(accuracy_svm); grid on;
title('SVM Test Data (Multi-Class)')
ylabel('Accuracy (%)');
xlabel('Iteration');

% Collect together the performance of neural net for test data sets.
figure(2)
boxplot(accuracy_nn); grid on;
title('NN Test Data (Multi-Class)')
ylabel('Accuracy (%)');
xlabel('Iteration');