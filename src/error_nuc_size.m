function [error,CELL_STAGE]=error_nuc_size(parameters,exp_data,pars,parnames)
% Simulates nuclear growth sequence and compute the relative error to exp.

%% Writing parameters from pars into parameters
parameters=include_parameters(parameters,pars,parnames);

%% We check the sanity of results
if parameters.k0>0 && parameters.P_ratio>0 && parameters.plaw>=0 && parameters.N_sat>=0
	[CELL_STAGE,flag]=nuclear_growth_sequence(exp_data,parameters);
else
	flag=0;
end

%% Computing error
if flag
    % Initialization
	n_div=numel(CELL_STAGE);

	% Range of stages to fit
	n_min=parameters.stage_min;
	n_max=parameters.stage_max;
	if n_max>n_div
		n_max=n_div;
	end

	ERR=zeros(1,n_div);
    % Total of number of experimental points
    n_pts=0;
    for n=n_min:n_max
				stage=CELL_STAGE(n);
				pts=exp_data(n).nucdata;
        n_pts=n_pts+length(pts);
        % We compute the relative error
        ERR(n)=ERR(n)+sum(((stage.nucsize(:)-pts(:)).^2.0)./(pts(:).^2.0));
    end
    % We then divide by number of experimental points for normalization
		error=sum(ERR(2:end))/n_pts;
else
    % If there was an error in the simulation, we yield a BIG error
	error=10^9.0;
end

end
