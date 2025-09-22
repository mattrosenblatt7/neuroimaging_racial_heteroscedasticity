function step01_run_levenes_test_fc(task_name, covar_regr, site_correction, balanced)
%{
Function to run Levene's test by race for functional connectivity data
    task_name, options: rsfMRI, MID, SST, nBack
    covar_regr: true, false
    site_correction: true, false
    balanced: true, false - whether to balance # participants by race
%}

addpath(genpath('../../tools'))

% load all data
load_path = '../../data';
save_path = '../output';
load(fullfile(load_path,['lowest_motion_', task_name, '.mat']))

% load in table and subset to columns of interest
tb = readtable(fullfile(load_path,'behavior.csv'));
tb = tb(:, {'src_subject_id', 'Race'});

% load in covariate data if option is true
if covar_regr
    tb_covar = readtable(fullfile(load_path, 'abcd_covariates.csv'));
    tb_covar.src_subject_id = cellfun(@(x) ['sub-NDARINV', x(9:end)], tb_covar.src_subject_id, 'UniformOutput', false);
    
    % merge covariate data into table
    tb = outerjoin(tb, tb_covar,'Type','Left','Keys', {'src_subject_id'}, 'MergeKeys', 1);
    
    % remove nan
    covar_keys = {'sex', 'age', 'demo_comb_income_v2'};
    missing_covars_idx = find(sum(isnan(tb{:, covar_keys}), 2)>0);
    tb(missing_covars_idx, :) = [];
elseif site_correction  % need to merge in site info
    tb_covar = readtable(fullfile(load_path, 'abcd_covariates.csv'));
    tb_covar.src_subject_id = cellfun(@(x) ['sub-NDARINV', x(9:end)], tb_covar.src_subject_id, 'UniformOutput', false);
    
    % merge covariate data into table
    tb = outerjoin(tb, tb_covar,'Type','Left','Keys', {'src_subject_id'}, 'MergeKeys', 1);
end

% find table locations matching each matrix subject ID
matching_tb_idx = cellfun(@(x) find(strcmp(tb.src_subject_id, x)), mat_sub_ids_final, 'UniformOutput', 0);
mats_with_missing_info_idx = find(cellfun(@isempty, matching_tb_idx));

% remove matrices with missing data in behavioral table
lowest_motion_vals(mats_with_missing_info_idx) = [];
lowest_motion_connectomes(:, mats_with_missing_info_idx) = [];
mat_sub_ids_final(mats_with_missing_info_idx) = [];
matching_tb_idx(mats_with_missing_info_idx) = [];
matching_tb_idx = cell2mat(matching_tb_idx);  % convert from cell to matrix
tb = tb(matching_tb_idx, :);

% make dictionary to map race to numbers
map = containers.Map({'Asian', 'Black','Non-Black Multiracial Or Other','White'},[1, 2, 3, 4]);

% cell function applies map to every element of cell
race = cell2mat(cellfun(@(x) map(x), tb.Race, 'UniformOutput', false));

% convert connectomes to edges
edges = lowest_motion_connectomes;


% regress out covariates if option is true
if covar_regr
    C_all = [tb{:, covar_keys}, lowest_motion_vals];
    Beta = ((C_all'*C_all)\C_all')*edges';
    edges = edges - (C_all*Beta)';
    save_prefix = 'covar_regr_';
else
    save_prefix = '';  % don't need to append anything to save name
end

% whether to correct for site with ComBat
% https://github.com/Jfortin1/ComBatHarmonization/tree/master/Matlab
if site_correction
    batch = cellfun(@(x) str2num( x((end-1):end) ), tb.site)'; %Batch variable for the scanner id
    edges = combat(edges, batch, [], 1);
    
    save_prefix = ['site_corr_', save_prefix];
end

%% Levene's test: compare Black and White



if balanced
    % adjust save name
    save_prefix = ['balanced_', save_prefix];
    
    % restrict edges to these participants
    edges_bw = edges(:, (race==2) | (race==4));
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
    edges_bw = edges(:, (race==2) | (race==4));
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
save(fullfile(save_path, [save_prefix, 'levenes_bw_results_', task_name, '.mat']),...
    'p_all_bw', 'f_all_bw', 'pID', 'pN', 'w_sd', 'w_var', 'b_sd', 'b_var')
