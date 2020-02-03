function [ CELL_STAGE , flag] = nuclear_growth_sequence(  exp_data, parameters  )
% Computes the growth of nuclear size for a sequence of divisions
%   We assumed pER to be consumed by Nucleus and created by cells
%   Based on michaelis-Menten kinetics of Nuclear growth from pER
%
% CELL_STAGE is a structure containing nuclear and pER size
% nucdata : Backup of exp. data (nuclear surface)
% errdata : Backup of exp. data (error on nuclear surface)
% nucsize : simulated nuclear surface
% persize : simulated pER volume
% maxsize : maximum nuclear size
% N_tot : available nuclear surface
% times : simulated times
%
% exp_data : a structure that contains
%				times : time points
%				nucdata : experimentally measured nuclear surface
%				errdata : error on experimental measure of nuclear surface

if nargin < 2
	parameters=nucleus_parameters();
    if nargin < 1
        error('Nuclear Growth Sequence requires experimental data')
    end
end
if ~isfield(parameters,'no_division')
    parameters.no_division=0;
end
if ~isfield(parameters,'include_init')
    parameters.include_init=0;
end

%% Reading experimental results
n_div=numel(exp_data);
initial_size=zeros(n_div,1);
for n=1:n_div
	initial_size(n)=exp_data(n).nucdata(1);
end

%% Initialization
% A flag to know whether the simulation was able to finish
flag=1;
% The structure to save the results
CELL_STAGE(n_div).nucsize=[];
for n=1:n_div
    CELL_STAGE(n).nucdata=exp_data(n).nucdata; % Backup of exp. data
    CELL_STAGE(n).errdata=exp_data(n).errdata; % Backup of exp. data
    CELL_STAGE(n).nucsize=0.0.*exp_data(n).nucdata; % Nuclear surface
    CELL_STAGE(n).persize=0.0.*exp_data(n).nucdata; % PER volume
    CELL_STAGE(n).maxsize=0.0;                      % Max N size
    CELL_STAGE(n).N_tot=0.0;                        % Available N size
    CELL_STAGE(n).time=exp_data(n).times;           % times simulated
end

%% Reading physical parameters
% k0 is the import rate in the absence of saturation
k0=parameters.k0;
% plaw is the power law between PER volume and N surface area
plaw=parameters.plaw;
% P_ratio is the ratio of PER available to make N
P_ratio=parameters.P_ratio;
% P_tot_0 is the initial PER quantity
P_tot_0=parameters.P_tot_0;
% P_sat is the PER saturation
N_sat=parameters.N_sat;
% change of k0 with time
delta_k0=parameters.delta_k0;

%% Other parameters
% Verbosity
verbose=parameters.verbose;

%% Computation
for n=1:n_div

    % Dividing pER
    if ~(parameters.no_division)
        P_tot=P_tot_0/(2.0^(n-1));
    else
        P_tot=P_tot_0;
    end

    % Converting pER to available Nucleus material
    P_available=P_ratio*P_tot;
    N_tot=(P_available)^plaw;

	% Counter only if verbose==1
	if verbose>0 && verbose<2
		disp(['Division ' num2str(n) ' out of ' num2str(n_div)])
    end

	% Initial values of N
	N_init=initial_size(n);

    %% If we include the initial Nuc building in the kinetics
    if parameters.include_init==1

        % We need to compute how long it takes
        % to make a nucleus of size N_init
        ti=get_initial_time(k0,N_tot,N_init,N_sat);

        % If that time is not realisting
        if ~isreal(ti) || 0 > ti
            warning('Not enough material to make nuclei of measured size !')
            flag=0;
            break
        end

        times=exp_data(n).times;
        t=times-min(times)+ti;

        nuc=get_nuclear_size( k0 , N_tot,  N_sat , t );


        per=P_tot-(((nuc-N_init)).^(1/plaw));
        t=t-ti;
    %% If we DO NOT include the initial Nuc building in the kinetics (default)
    else

        times=exp_data(n).times;
        t=times-min(times);

         if 0 > N_tot-N_init
            warning('Not enough material to make nuclei of measured size !')
            flag=0;
            break
        end

        nuc=N_init+get_nuclear_size( k0 , N_tot-N_init, N_sat, t );
        per=P_tot-(((nuc-N_init)).^(1/plaw));

    end

	%% Saving results
	CELL_STAGE(n).nucsize=nuc;
  CELL_STAGE(n).persize=per;
	CELL_STAGE(n).time=t+min(times);
  CELL_STAGE(n).N_tot=N_tot;
  CELL_STAGE(n).P_tot=P_tot;
  CELL_STAGE(n).flag=flag;

	% changing k0
	k0=k0+delta_k0;
end

per_ratio=parameters.Measured_P0/CELL_STAGE(1).persize(1);

for n=1:n_div
    CELL_STAGE(n).rescaled_per=abs((CELL_STAGE(n).persize)*per_ratio);
end


end

%% Wrapper functions
% You should replace those by calls to your own functions
% If you want to test other models of nuclear growth
function [ti]=get_initial_time(k0,N_tot,N,N_sat)
%		get_initial_time : wrapper function
% computes the time necessary to reach P=P_init knowing other parameters
% good thing is, inverting LambertW is easy
ti=(N-N_sat*log(1-N/N_tot))/k0;
end

function [ nuc ] = get_nuclear_size( k0 , N_tot , N_sat, t )
% Size of nucleus at time t following saturated import kinetics
%   k0 is the unsaturated rate of input
%   N_sat is the saturating per
%   N_tot is the material available at t=0
%   t is time
%   Output : nuc : nuclear size
nuc=N_tot-N_sat*get_lambert_w( N_tot * exp( (N_tot - k0*t)/N_sat )/N_sat);
end
