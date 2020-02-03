%% fitting options
options.MaxFunEvals=100000;
options.MaxIter=100000;

%% Input
parameters=nucleus_parameters();
data=load('NucData_12_07.mat','NucData');
exp_data=data.NucData;
parameters.verbose=2;

%% Parameters we choose to fit
parnames={'k0','P_ratio','plaw','N_sat'};
for i=1:numel(parnames)
  par0(i)=getfield(parameters,parnames{i});
end

%% Best fit
% definition of the error function
error_function = @(pars) error_nuc_size(parameters,exp_data,pars,parnames);
% minimization
[BEST_PARS,error] = fminsearch(error_function, par0,options);
% writing the parameters
parameters=include_parameters(parameters,BEST_PARS,parnames);

%% Plotting and saving
parameters.save.prefix='best_fit_stage_';
parameters.plot.name='Best fit';
[CELL_STAGE,flag]=nuclear_growth_sequence(exp_data,parameters);
if flag
    save_cell_stage(CELL_STAGE,parameters.save,exp_data);
    plot_nucleus_sim(CELL_STAGE,parameters.plot,exp_data);
end
