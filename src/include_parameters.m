function [ parameters ] = include_parameters( parameters,pars,parnames )
%Take some parameters and modifies them according to an array (pars)
% the name of the parameters to be changed is in parnames
% This is very conveniently non-specific

n_pars=numel(pars);

for i=1:n_pars
  parameters=setfield(parameters,parnames{i},pars(i));
end

end
