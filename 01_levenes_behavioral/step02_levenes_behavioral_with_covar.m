%{
Run levene's test for behavioral data to compare variance across groups,
including sex/age/income as covariates
%}

% load all data
load_path = "../output";
tb = readtable(fullfile(load_path,"combined_genon_behaviors.csv"));

% check for missing race data
disp(['Number of participants with missing race data: ', num2str(sum(cellfun(@(x) isempty(x), tb.Race)))])

% get variables of interest
variable_names = tb.Properties.VariableNames';
variable_names = variable_names(8:end);    % ignore info/covariate columns

% loop over all variables of interest
p_behavior = NaN + zeros(length(variable_names), 1);
f_behavior = NaN + zeros(length(variable_names), 1); 
b_sd = NaN + zeros(length(variable_names), 1); 
w_sd = NaN + zeros(length(variable_names), 1); 
N = NaN + zeros(length(variable_names), 1); 
for var_idx = 1:length(variable_names) 
    
    % restrict data to variables of interest
    tb_tmp = tb(:, {'Race', 'sex', 'age', 'demo_comb_income_v2', variable_names{var_idx}});

    % restrict to only variables with data
    good_idx = find(sum(isnan(tb_tmp{:, {'sex', 'age', 'demo_comb_income_v2', variable_names{var_idx}}}), 2)==0);
    tb_tmp = tb_tmp(good_idx, :);

    % covariate regression
    y = tb_tmp{:, variable_names{var_idx} };
    C_all = tb_tmp{:, {'sex', 'age', 'demo_comb_income_v2'}};
    Beta = ((C_all'*C_all)\C_all')*y;
    y = y - (C_all*Beta);

    % run levene's test
    [p, stats] = vartestn(tb_tmp{:, variable_names{var_idx} },tb_tmp{:, 'Race'},'TestType','LeveneAbsolute', 'Display','off');
    p_behavior(var_idx) = p;
    f_behavior(var_idx) = stats.fstat;
    N(var_idx) = stats.df(end) + 2;
    
    % calculate standard deviation by race
    b_sd(var_idx) = nanstd(tb_tmp{strcmp(tb_tmp{:, 'Race'}, 'Black'), variable_names{var_idx} });
    w_sd(var_idx) = nanstd(tb_tmp{strcmp(tb_tmp{:, 'Race'}, 'White'), variable_names{var_idx} });
    
end
results = table(variable_names, p_behavior, f_behavior, b_sd, w_sd, N);
results = sortrows(results,{'f_behavior'},{'descend'});
writetable(results, '../output/levenes_behavioral_results_COVAR.csv');