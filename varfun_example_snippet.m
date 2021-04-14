n_animals = 10;% number of individuals
n_trials = 40;% number of trials per individual
animal = repmat(1:n_animals,n_trials,1); animal = animal(:);% vector of animal identitities
trial = repmat(1:n_trials,1,n_animals); trial = trial(:);% vector of trial numbers
%useful functions
Logit = @(x) log(x./(1-x));% function for converting from probability to log(odds)
Logistic = @(x) 1./(1+exp(-x));% function for converting from log(odds) to probability
%individual biases
%random normal distribution in logit space: normrnd(mean, standard deviation, rows, columns)
pre_bias_animal = Logistic(normrnd(Logit(0.25),1,n_animals,1));% starting bias (proportion correct)
%generate choices using learning curve parameters
%binomial sample: binomrnd(observerations per datapoint, probabilities)
choice_correct = binornd(1,Logistic(Logit(pre_bias_animal(animal))+0.1*(trial-1)));
data_table = table(logical(choice_correct), ...
                    cellstr(num2str(animal)), ...
                    double(trial), ...
    'VariableNames', {'choice','ID','trial'} ...
    )

%% Calculate mean and normal CI for each animal
pop_mean = varfun( @(x)mean(double(x))+[0,-1,1].*1.95.*std(double(x))/sqrt(length(x)), ...
                   data_table(:,[1,3]), ...
                   'GroupingVariables','trial');
plot(pop_mean.trial, pop_mean.Fun_choice)
