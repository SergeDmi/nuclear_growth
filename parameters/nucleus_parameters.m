function [ parameters ] = nucleus_parameters(  )
% Options for modeling nucleus growth
%   just a file to store options
%   Better save different options under a different name
%   (i.e. different file and function name)

%% Experimental parameter for scaling pER to experimental units
% This is just the measured value for pER volume initially !
parameters.Measured_P0=15;

%% Physical parameters
% Not fited
% This is the total quantity of pER in units : ( mu^2 ) ^ (1/plaw)
parameters.P_tot_0=10000000;
% Fited
parameters.k0=100;
parameters.P_ratio=0.1;
parameters.plaw=0.6;
parameters.N_sat=240;

% increase of k0 with Stage
parameters.delta_k0=0.0;
%% Mode selection
parameters.no_division=0;
parameters.include_init=0;

%% Fitting behaviour
parameters.stage_min=2;
parameters.stage_max=9;
%% Program behaviour
parameters.verbose=0;
parameters.do_save=1;

%% saving behaviour
parameters.save.folder='results/';
parameters.save.prefix='nucleus_stage';

%% plotting behaviour
parameters.plot.title='Nucleus growth';

end
