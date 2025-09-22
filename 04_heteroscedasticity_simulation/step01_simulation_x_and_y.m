%{
This code runs simulations to evaluate MSE as a function of predictor variance and predicted variance 
%}

%% 
clear

%% set variances for heteroskedastic and homoskedastic group x and y data
% empirical avg values which noise will be added to later give more spread to the results
all_x_hom_var_expected = 0.5:0.5:5;
all_x_het_var_expected = all_x_hom_var_expected;
all_y_hom_var_expected = 20.5:0.5:25;
all_y_het_var_expected = all_y_hom_var_expected;

%% other params
% empirical slope for the models
beta1 = 2;

% baseline X_hom
n_hom = 1000;
X_hom = 10*(1:n_hom)/n_hom-10/n_hom;
% add homoskedastic noise variance equal to 1
X_hom = X_hom + randn(1,n_hom);

% baseline X_het
n_het = n_hom;
X_het = 10*(1:n_het)/n_het-10/n_het;
% add heteroskedastic noise with average variance equal to 1
X_het = X_het + normrnd(0,sqrt((1/n_het):(2/n_het):(1/n_het+2/n_het*(n_het-1))));

% normalize to have variance 1 and mean 0
X_hom = (X_hom-mean(X_hom))/std(X_hom,1); % make variance 1 and mean 0
X_het = (X_het-mean(X_het))/std(X_het,1); % make variance 1 and mean 0

%% initialize variables to save
% measured variances, mse, and model parameters
y_het_var = zeros(length(all_x_het_var_expected),length(all_x_hom_var_expected),length(all_y_het_var_expected),length(all_y_hom_var_expected));
y_hom_var = y_het_var;
x_het_var = y_het_var;
x_hom_var = y_het_var;
het_mse = y_het_var;
hom_mse = y_het_var;
all_beta0 = y_het_var;
all_beta1 = y_het_var;

%% run simulation
for x_het_ind = 1:length(all_x_het_var_expected)
    % set expected variance for the heteroskedastic X data
    x_het_var_exp = all_x_het_var_expected(x_het_ind);
    for x_hom_ind = 1:length(all_x_hom_var_expected)
        % set expected variance for the homoskedastic X data
        x_hom_var_exp = all_x_hom_var_expected(x_hom_ind);
        for y_het_ind = 1:length(all_y_het_var_expected)
            % set expected variance for the heteroskedastic Y data
            y_het_var_exp = all_y_het_var_expected(y_het_ind);
            for y_hom_ind = 1:length(all_y_hom_var_expected)
                % set expected variance for the homoskedastic Y data
                y_hom_var_exp = all_y_hom_var_expected(y_hom_ind);

                % set X_het data to have expected variance plus some noise
                X_het = X_het/sqrt(var(X_het,1)/(x_het_var_exp+normrnd(0,0.1,1)));

                % set X_hom data to have expected variance plus some noise
                X_hom = X_hom/sqrt(var(X_hom,1)/(x_hom_var_exp+normrnd(0,0.1,1)));

                % set the expected noise variance to get the right y_het variance based on the empirical model
                noise_het_var_exp = y_het_var_exp - beta1^2*x_het_var_exp;
                % generate sequence of increasing standard deviation such
                % that the expected variance across all noise points added
                % is equal to noise_het_var_exp
                a = n_het*sqrt(noise_het_var_exp/((n_het-1)*sum((X_het-min(X_het)).^2)));
                het_sig = a*(X_het-min(X_het));
                % generate heteroskedastic noise to add to y_het
                noise_het = normrnd(0,het_sig,size(X_het));
                % generate y_het
                y_het = beta1*X_het + noise_het;

                % generate y_hom data
                % set the expected noise variance to get the right y_hom variance based on the empirical model
                noise_hom_var_exp = y_hom_var_exp - beta1^2*x_hom_var_exp;
                % generate homoskedastic noise to add to y_hom such that
                % the expected variance across all noise points added is
                % equal to noise_hom_var_exp
                noise_hom = normrnd(0,sqrt(n_hom/(n_hom-1)*noise_hom_var_exp),size(X_hom));
                % generate y_hom 
                y_hom = beta1*X_hom + noise_hom;

                % create model
                design_mat = [ones(n_hom+n_het,1) [X_hom'; X_het']];
                beta = design_mat\[y_hom'; y_het'];
                pred = (design_mat*beta)';

                % calculate MSE
                het_mse(x_het_ind,x_hom_ind,y_het_ind,y_hom_ind) = mean((y_het-pred((end-n_het+1):end)).^2);
                hom_mse(x_het_ind,x_hom_ind,y_het_ind,y_hom_ind) = mean((y_hom-pred(1:n_hom)).^2);

                % add to saved variables
                y_het_var(x_het_ind,x_hom_ind,y_het_ind,y_hom_ind) = var(y_het,1);
                y_hom_var(x_het_ind,x_hom_ind,y_het_ind,y_hom_ind) = var(y_hom,1);
                x_het_var(x_het_ind,x_hom_ind,y_het_ind,y_hom_ind) = var(X_het,1);
                x_hom_var(x_het_ind,x_hom_ind,y_het_ind,y_hom_ind) = var(X_hom,1);
                all_beta0(x_het_ind,x_hom_ind,y_het_ind,y_hom_ind) = beta(1);
                all_beta1(x_het_ind,x_hom_ind,y_het_ind,y_hom_ind) = beta(2);
            end
        end
    end
end

%% save
all_results = struct;
all_results.y_het_var = y_het_var(:);
all_results.y_hom_var = y_hom_var(:);
all_results.x_het_var = x_het_var(:);
all_results.x_hom_var = x_hom_var(:);
all_results.het_mse = het_mse(:);
all_results.hom_mse = hom_mse(:);
all_results.beta0 = all_beta0(:);
all_results.beta1 = all_beta1(:);
save('results_X_Y.mat','all_results')