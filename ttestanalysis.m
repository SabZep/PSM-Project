% ONE SAMPLE t-test
% 
% This script performs 2nd level analysis, i.e, a one sample t-test 

% Initialise SPM 
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');

clear matlabbatch


% results dir
matlabbatch{1}.spm.stats.factorial_design.dir = {'/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_betweenconditions'};

%% Get images for 
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_01/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_02/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_03/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_04/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_05/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_06/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_07/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_08/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_09/beta_0001.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_10/beta_0001.nii,1'
                                                          };
%% Create one-sample t-test
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run',matlabbatch);

% Model estmation 
matlabbatch{1}.spm.stats.fmri_est.spmmat = {'/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_betweenconditions/SPM.mat'};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run',matlabbatch);

clear matlabbatch


%% T-contrast specification

connames = 'Above chance';
convecs  = 1;

% Number of contrasts
numCons  = 1;

% Allocate t-contrast structure
matlabbatch{1}.spm.stats.fmri_est.spmmat = {'/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_betweenconditions/SPM.mat'};

for c = 1:numCons
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name    = connames;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = convecs;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
end 
% run job
spm_jobman('run',matlabbatch);


% clear job variable
clear matlabbatch



%% Run one-sample t-test for cor between fMRI-based RDM & stimulus type RDM 


% results dir
matlabbatch{1}.spm.stats.factorial_design.dir = {'/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_stimulustype'};

%% Get images for 
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_01/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_02/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_03/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_04/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_05/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_06/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_07/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_08/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_09/beta_0002.nii,1'
                                                          '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_10/beta_0002.nii,1'
                                                          };
%% Create one-sample t-test
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run',matlabbatch);

% Model estmation 
matlabbatch{1}.spm.stats.fmri_est.spmmat = {'/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_stimulustype/SPM.mat'};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run',matlabbatch);


matlabbatch{5}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';


spm_jobman('run',matlabbatch);

% clear job variable
clear matlabbatch


% SPM.mat where further manually analysed in the spm12 using the anatomy
% toolbox


