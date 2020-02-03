% Input
parameters=nucleus_parameters();
data=load('NucData_12_07.mat','NucData');
exp_data=data.NucData;
parameters.verbose=2;
% Calling the good stuff
[CELL_STAGE,flag]=nuclear_growth_sequence(exp_data,parameters);
if flag
    save_cell_stage(CELL_STAGE,parameters.save,exp_data);
    plot_nucleus_sim(CELL_STAGE,parameters.plot,exp_data);
end
