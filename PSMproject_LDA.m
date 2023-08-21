% PSM Decoding LDA project: Zeynep Elif Sab, 19/08/2023

% This code decodes three imagery conditions(ImagPress,ImagFlutt,ImagVibro) from two regions of
% interests: right BA2 and right SII.

%Before you run the code, do some changes here:
% Change this path to the motherfile that contains all the 1st_level_good_bad_Imag_ files of all subjects: 
homefolder = '/Users/sab/MATLAB/decoding/';
% Change these paths to our coregistered ROI masks:
roi_1_mask_path = '/Users/sab/MATLAB/decoding/rois_lda/rRight_BA2.nii'; 
roi_2_mask_path = '/Users/sab/MATLAB/decoding/rois_lda/rRight_SII.nii';
%Also you may need to install MBoxtest and HZmvntest to run this code.
%Links:https://de.mathworks.com/matlabcentral/fileexchange/17931-hzmvntest
%https://www.mathworks.com/matlabcentral/fileexchange/2733-mboxtest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%And that is it!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subj= 10; %number of subjects

X_all = cell(10, 1); %storage for every subject's data matrix, for plotting purposes.
all_accuracy =  []; % Also let's create an array to store all subjects' accuracy scores.

% Storage for all subject's averaged coefficients for every class pairing,later will
% be used for plotting:
all_constants_12 = [];
all_linears_12 = zeros(subj,2); % we have 2 linear coefficients per class pairing.
all_constants_13 = [];
all_linears_13 = zeros(subj,2); 
all_constants_23 = [];
all_linears_23 = zeros(subj,2); 

% Loop over subjects from sub01 to sub10!
for subject_num = 1:subj
    % Convert subject number to a string, taking care of zeros (sub01
    % sub02... sub10)
    subj_str = sprintf('sub%02d', subject_num);
    
    % Selecting the 1st level file of the current subject:
    base_path = [homefolder, '1st_level_good_bad_Imag_', subj_str, '/'];

% There is one beta file per condition per each run! We have six runs
% and we are interested in three conditions: ImagPress,ImagFlutt,ImagVibro.
% The following array contains the beta file numbers of these three conditions 
% respectively as sets of six. 
beta_numbers = [4, 15, 26, 37, 48, 59, 5, 16, 27, 38, 49, 60, 6, 17, 28, 39, 50, 61]; 

X_ROI1 =  []; %empty array for average betas from ROI 1 (rBA2)
X_ROI2 =  []; %empty array for average betas from  ROI 2 (rSII) 

%% Extracting betas from our ROI's!
% Load the mask files and read them in SPM
roi_1_mask = spm_vol(roi_1_mask_path);
roi_1_read = spm_read_vols(roi_1_mask);

roi_2_mask = spm_vol(roi_2_mask_path);
roi_2_read = spm_read_vols(roi_2_mask);

for i = 1:length(beta_numbers) %let's run through the beta files of our conditions!

% Select the current beta file
beta_number = beta_numbers(i);
beta_file_path = [base_path, 'beta_', sprintf('%04d', beta_number), '.nii'];

% Load it and read it!
beta_file = spm_vol(beta_file_path);
beta_read = spm_read_vols(beta_file);

% Get the voxels in our ROI's
roi_1_voxels = find(roi_1_read > 0); %voxels are coded binary, 1 belongs to our roi and 0 are irrelevant voxels. 
roi_2_voxels = find(roi_2_read > 0);

% Extract beta values from all of the voxels in our ROI's
beta_values_roi_1 = beta_read(roi_1_voxels);
beta_values_roi_2 = beta_read(roi_2_voxels);

% Now let's average the beta values for each ROI
average_beta_roi_1 = nanmean(beta_values_roi_1);
average_beta_roi_2 = nanmean(beta_values_roi_2);

% Let's store the average beta value of each iteration! So when the loop is
% over, we'll have 18 beta values in each array.
% These are: in sets of six respectively: ImagPress,ImagFlutt,ImagVibro betas for our
% subject.
X_ROI1 = vertcat(X_ROI1, average_beta_roi_1); 
X_ROI2 = vertcat(X_ROI2, average_beta_roi_2);

end
%% Henze-Zirkler's Test for Multivariate Normality
%Let's test the assumption for multivariate normal distribution.

%Getting the data in the right structure for the test:
HZdata = [ X_ROI1(1:6,:) X_ROI1(7:12,:) X_ROI1(13:18,:); X_ROI2(1:6,:) X_ROI2(7:12,:) X_ROI2(13:18,:)]; 
HZmvntest(HZdata);

%Data analyzed has a normal distribution!

%%  Let's test the similarity of covariance matrices across classes. 
% LDA assumes shared covariance matrices and if that assumption is not met, we
% can consider QDA or GNB with regularizing.

%We'll use MBoxtest, let's get our data in the right structure for it:
Z = [ones(6,1) X_ROI1(1:6,:) X_ROI1(7:12,:) X_ROI1(13:18,:);ones(6,1)*2 X_ROI2(1:6,:) X_ROI2(7:12,:) X_ROI2(13:18,:)]; 
MBoxtest(Z,0.05);

%Covariance matrices are not significantly different.
%LDA it is!
%% Let's make our data prettier for LDA! 
% Y is labels and X is our data with
% two columns for each feature.

Y = [zeros(6,1); ... % ImagPress: "0"
    ones(6,1); ... % ImagFlutt: "1"
    ones(6,1)*2]; % ImagVibro: "2"

X = [X_ROI1  X_ROI2]; 

X_all{subject_num} = X; %remember that we are storing each subject's data matrix for plotting purposes.

%% LDA and cross validation. 
%We have six runs and we will do a leave one run out cross validation. Each
%fold will select one run betas as test data, and the remaining 5 runs betas as train data.
%For each fold, let's store the accuracy and model coefficients! We will
%later average the accuracies to calculate the subject's overall accuracy, and
%model coefficients will be averaged for plotting. 
all_folds_accuracy =  [];

all_folds_constant_12 = [];  
all_folds_linear_12 = zeros(6,2); %we have 6 folds and 2 linear coefficients to store in each fold.
all_folds_constant_13 =  [];
all_folds_linear_13 = zeros(6,2); 
all_folds_constant_23 =  [];
all_folds_linear_23 = zeros(6,2);

for fold = 1:6
   
    test_indices = [fold, fold+6, fold+12]; %these are the indices of the betas in the same run
    train_indices = setdiff(1:18, test_indices); %these are the remaining indices
    X_test = X(test_indices, :);
    Y_test = Y(test_indices);
    X_train = X(train_indices, :);  
    Y_train = Y(train_indices);

discr = fitcdiscr(X_train,Y_train);
predicted = predict(discr, X_test);
fold_accuracy = mean(Y_test == predicted); %accuracy for the current fold

all_folds_accuracy = horzcat(all_folds_accuracy, fold_accuracy);%accuracies of all 6 folds will be stored here

%storing the constant coefficients of all 6 folds, for each class pairing
all_folds_constant_12 = horzcat(all_folds_constant_12, discr.Coeffs(1,2).Const);
all_folds_constant_13 = horzcat(all_folds_constant_13, discr.Coeffs(1,3).Const);
all_folds_constant_23 = horzcat(all_folds_constant_23, discr.Coeffs(2,3).Const);
%storing the linear coefficients of all 6 folds, for each class pairing
all_folds_linear_12(fold,:)=[discr.Coeffs(1,2).Linear(1), discr.Coeffs(1,2).Linear(2)];
all_folds_linear_13(fold,:)=[discr.Coeffs(1,3).Linear(1), discr.Coeffs(1,3).Linear(2)];
all_folds_linear_23(fold,:)=[discr.Coeffs(2,3).Linear(1), discr.Coeffs(2,3).Linear(2)];

end

subject_accuracy= mean(all_folds_accuracy);   
all_accuracy = horzcat(all_accuracy, subject_accuracy); %this is all subjects' accuracy scores.

%averaged constant coefficients of all subjects, for each class pairing
all_constants_12 = horzcat(all_constants_12, mean(all_folds_constant_12));
all_constants_13 = horzcat(all_constants_13, mean(all_folds_constant_13));
all_constants_23 = horzcat(all_constants_23, mean(all_folds_constant_23));
%averaged linear coefficients of all subjects, for each class pairing
all_linears_12(subject_num,:)= [mean(all_folds_linear_12(:,1)), mean(all_folds_linear_12(:,2))];
all_linears_13(subject_num,:)= [mean(all_folds_linear_13(:,1)), mean(all_folds_linear_13(:,2))];
all_linears_23(subject_num,:)= [mean(all_folds_linear_23(:,1)), mean(all_folds_linear_23(:,2))];

end %subject iteration ends!!!

%LDA for all participants ends here! We gathered the average accuracies for
%all participants.

%Also we have gathered the average constant and linear coefficients for
%each participant.

%% One sample t-test
% Now let's see if our classifier performs significantly better than the
% chance level, which is 0.33 for us since we decode 3 conditions.
[h,p,ci,stats] = ttest(all_accuracy, 0.33, 'tail','right');

%Good to know when reporting the t-test results:
mean_acc = mean(all_accuracy);
std_acc = std(all_accuracy);

if h==0 
     fprintf('The classifier does not perform significantly better than chance level');
else
    fprintf('The classifier performs significantly better than chance level');
end

%% Plotting

%Let's plot for the subjects for the best and worst accuracy.
[best_accurate, max_ind] = max(all_accuracy);  
[worst_accurate ,min_ind] = min(all_accuracy); 

n= max_ind; 
%!!!!T!!!!!!!!!!!THIS PLOTS THE PARTICIPANT WITH maximum LDA accuracy!! 
% If you want to see the participant with worst accuracy, change to
%n= min_ind !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

X = X_all{n}; %Selecting the data matrix of the subject of interest.
figure;
hold on
scatter(X(Y==0,1),X(Y==0,2),'c')
scatter(X(Y==1,1),X(Y==1,2),'b')
scatter(X(Y==2,1),X(Y==2,2),'r')
xlabel('Average beta value in rBA2 ') 
ylabel('Average beta value in rSII')
legend({'ImagPress' 'ImagFlutt' 'ImagVibro'},'Location','NW') 
xlim([-2 3])
ylim([-2 3])
if n== max_ind
title(['Highest Accuracy: ', sprintf('%.2f%%', all_accuracy(n) * 100)]);
else
title(['Lowest Accuracy: ', sprintf('%.2f%%', all_accuracy(n) * 100)]);
end

% to plot boundaries, we are now looking at three different lines:
% between category 1 and 2
x1 = X(:,1);
x2 = X(:,2);
f = @(x1,x2) all_constants_12(n) + all_linears_12(n,1)*x1 + all_linears_12(n,2)*x2;
h2 = fimplicit(f,[-2 3 -2 3]);
h2.Color = 'k';
h2.LineWidth = 2;
h2.DisplayName = 'boundary ImagPress/ImagFlutt';

% between category 1 and 3
x1 = X(:,1);
x2 = X(:,2);
f = @(x1,x2) all_constants_13(n) + all_linears_13(n,1)*x1 + all_linears_13(n,2)*x2;
h2 = fimplicit(f,[-2 3 -2 3]);
h2.Color = 'b';
h2.LineWidth = 2;
h2.DisplayName = 'boundary ImagPress/ImagVibro';

% between category 2 and 3
x1 = X(:,1);
x2 = X(:,2);
f = @(x1,x2) all_constants_23(n) + all_linears_23(n,1)*x1 + all_linears_23(n,2)*x2;
h2 = fimplicit(f,[-2 3 -2 3]);
h2.Color = 'r';
h2.LineWidth = 2;
h2.DisplayName = 'boundary ImagFlutt/ImagVibro';

