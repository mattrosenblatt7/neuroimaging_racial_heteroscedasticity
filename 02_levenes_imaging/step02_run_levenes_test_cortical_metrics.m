function step02_run_levenes_test_cortical_metrics(modality, covar_regr, site_correction, balanced)
%{
Function to run Levene's test by race for functional connectivity data
    modality, options: mri_y_dti_fa_fs_wm_dst, mri_y_smr_t1_gray_dst
    covar_regr: true, false
    site_correction: true, false
    balanced: true, false - whether to balance # participants by race
%}
addpath(genpath('../../tools'))

fname = [modality, '.csv']; 
% load all data
load_path = '../../data';
save_path = '../output';

% load FA data and restrict to baseline
brain_tb = readtable(fullfile(load_path,fname));  % mri_y_dti_fa_fs_wm_dst.csv
brain_tb = brain_tb( strcmp(brain_tb.eventname, 'baseline_year_1_arm_1'), :);
brain_tb.src_subject_id = cellfun(@(x) ['sub-NDARINV', x(9:end)],...
    brain_tb.src_subject_id, 'UniformOutput', 0);

tb_behavior = readtable(fullfile(load_path,'behavior.csv'));

% merge tables
[tb,ileft,iright] = outerjoin(brain_tb,tb_behavior(:, {'src_subject_id', 'Race'}),'Type','left');

% remove rows with missing race data
missing_race_idx = find(cellfun(@(x) isempty(x), tb.Race));
tb(missing_race_idx, :) = [];

% remove rows with missing FA measure
atlas_tb_idx = 3:153;  % MANUALLY SET WHICH COLUMNS OF TABLE ARE IMAGING MEASURE
missing_fa_idx = find(sum(isnan(tb{:, atlas_tb_idx}), 2)>0);
tb(missing_fa_idx, :) = [];

% REMOVE ROWS FAILING QC
tb_qc = readtable(fullfile(load_path, 'mri_y_qc_incl.csv'));  % load QC table
tb_qc = tb_qc( strcmp(tb_qc.eventname, 'baseline_year_1_arm_1'), :);  % restrict to baseline visit
if contains(fname, '_t1_')
    tb_qc = tb_qc(tb_qc.imgincl_t1w_include==1, :);  % restrict to passed QC only
elseif contains(fname, '_fa_')
    tb_qc = tb_qc(tb_qc.imgincl_dmri_include ==1, :);  % restrict to passed QC only
end
qc_passed_ids = cellfun(@(x) ['sub-NDARINV', x(9:end)],tb_qc.src_subject_id, 'UniformOutput', 0);


% now apply IDs to table so it only includes scans passing QC
tb = tb(ismember(tb.src_subject_id_brain_tb, qc_passed_ids), :);

% make dictionary to map yes/no to 1/0
map = containers.Map({'Asian', 'Black','Non-Black Multiracial Or Other','White'},[1, 2, 3, 4]);

% cell function applies map to every element of cell
race = cell2mat(cellfun(@(x) map(x), tb.Race, 'UniformOutput', false));


% load in covariate data if option is true
if covar_regr
    tb_covar = readtable(fullfile(load_path, 'abcd_covariates.csv'));
    tb_covar.src_subject_id = cellfun(@(x) ['sub-NDARINV', x(9:end)], tb_covar.src_subject_id, 'UniformOutput', false);
    
    % merge covariate data into table
    tb.Properties.VariableNames{1} = 'src_subject_id';  % rename to match covar table
    tb = outerjoin(tb, tb_covar,'Type','Left','Keys', {'src_subject_id'}, 'MergeKeys', 1);
    
    % remove nan
    covar_keys = {'sex', 'age', 'demo_comb_income_v2'};
    missing_covars_idx = find(sum(isnan(tb{:, covar_keys}), 2)>0);
    tb(missing_covars_idx, :) = [];
    race(missing_covars_idx) = [];
elseif site_correction
    tb_covar = readtable(fullfile(load_path, 'abcd_covariates.csv'));
    tb_covar.src_subject_id = cellfun(@(x) ['sub-NDARINV', x(9:end)], tb_covar.src_subject_id, 'UniformOutput', false);
    
    % merge covariate data into table
    tb.Properties.VariableNames{1} = 'src_subject_id';  % rename to match covar table
    tb = outerjoin(tb, tb_covar,'Type','Left','Keys', {'src_subject_id'}, 'MergeKeys', 1);
    
    % remove nan
    covar_keys = {'sex', 'age', 'demo_comb_income_v2'};
    missing_covars_idx = find(sum(isnan(tb{:, covar_keys}), 2)>0);
    tb(missing_covars_idx, :) = [];
    race(missing_covars_idx) = [];    
end

% extract values of interest
imaging_measures = tb{:, 3:153}';

% regress out covariates if option is true
if covar_regr
    C_all = [tb{:, covar_keys}];
    Beta = ((C_all'*C_all)\C_all')*imaging_measures';
    imaging_measures = imaging_measures - (C_all*Beta)';
    save_prefix = 'covar_regr_';
else
    save_prefix = '';  % don't need to append anything to save name
end


% whether to correct for site with ComBat
% https://github.com/Jfortin1/ComBatHarmonization/tree/master/Matlab
if site_correction
    batch = cellfun(@(x) str2num( x((end-1):end) ), tb.site)'; %Batch variable for the scanner id
    imaging_measures = combat(imaging_measures, batch, [], 1);
    
    save_prefix = ['site_corr_', save_prefix];
end

%% Levene's test: compare Black and White

if balanced
    % adjust save name
    save_prefix = ['balanced_', save_prefix];
    
    % restrict edges to these participants
    edges_bw = imaging_measures(:, (race==2) | (race==4));
    race_bw = race(race==2 | race==4);
    
    % subsample so equal number of Black and White participants
    n_b = sum(race_bw==2);  % find number of Black participants
    b_idx = find(race_bw==2);  % find indices w Black participants
    n_w = sum(race_bw==4);  % find number of White participants
    w_idx = find(race_bw==4);  % find indices w White participants
    shuffle_w_idx = w_idx(randperm(n_w));  % shuffle w indices
    new_subsample_idx = [b_idx; shuffle_w_idx(1:n_b)];  % get new sample
    edges_bw = edges_bw(:, new_subsample_idx);
    race_bw = race_bw(new_subsample_idx);
    
else  
    % restrict edges to these participants
    edges_bw = imaging_measures(:, (race==2) | (race==4));
    race_bw = race(race==2 | race==4);
end

% loop over all edges and perform test
p_all_bw = zeros(size(edges_bw, 1), 1);
f_all_bw = zeros(size(edges_bw, 1), 1);
for idx = 1:size(edges_bw, 1)
    
    if mod(idx, 5000)==0
        disp(idx)
    end
    
    [p, stats] = vartestn(edges_bw(idx, :)',race_bw,'TestType','LeveneAbsolute', 'Display','off');
    p_all_bw(idx) = p;
    f_all_bw(idx) = stats.fstat;
end

% compute standard deviation of each edge by race
b_sd = std(edges_bw(:, race_bw==2), [], 2);
w_sd = std(edges_bw(:, race_bw==4), [], 2);
b_var = var(edges_bw(:, race_bw==2), [], 2);
w_var = var(edges_bw(:, race_bw==4), [], 2);

% Dustin's FDR script
[pID,pN] = FDR(p_all_bw, .05);
% first is p value needed to be FDR corrected

%{
Save file with several variables:
    p_all_bw: Levene's test p values for rsfMRI edges in Black vs White
    participants
    pID: less strict p value FDR threshold
    pN: more strict p value FDR threshold
    w_sd(var): standard deviations (variance) of edges in White
    participants
    w_sd(var): standard deviations (variance) of edges in Black
    participants
%}
save(fullfile(save_path, [save_prefix, 'levenes_bw_results_', strrep(fname, '.csv', ''), '.mat']),...
    'p_all_bw', 'f_all_bw', 'pID', 'pN', 'w_sd', 'w_var', 'b_sd', 'b_var')

