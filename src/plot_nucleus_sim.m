function [  ] = plot_nucleus_sim( CELL_STAGE , parameters, exp_data )
% Plots results from nuclear_growth_sequence
%   just that

n_div=numel(CELL_STAGE);

max_sim =zeros(1,n_div);
ntot_sim=zeros(1,n_div);
max_exp =zeros(1,n_div);
err_exp =zeros(1,n_div);

for n=2:n_div
    % Gathering values
    max_exp(n) =max(CELL_STAGE(n).nucdata);
    err_exp(n) =max(CELL_STAGE(n).errdata);
    max_sim(n) =max(CELL_STAGE(n).nucsize);
    ntot_sim(n)=CELL_STAGE(n).N_tot;
    % Finding up to which point the simulation went
    if CELL_STAGE(n).flag
        n_max=n;
    end
end

figure
hold all
for n=2:n_max
    stage=CELL_STAGE(n);
    if nargin>1
        pts=exp_data(n).nucdata;
        err=exp_data(n).errdata;
        t=exp_data(n).times;
        errorbar(t,pts,err)
        %t=stage.time(1)+(1:length(pts))-1;
        %scatter(t,pts)
    end
    plot(stage.time,stage.nucsize,'k','LineWidth',1.5)


end
xlabel('Time')
ylabel('Nuclear Size')
if isfield(parameters,'title')
    title(parameters.title)
end

figure
hold all
errorbar(2:n_max,max_exp(2:n_max),err_exp(2:n_max),'r')
scatter(2:n_max,max_sim(2:n_max),'b')
hold all
plot(2:n_max,ntot_sim(2:n_max),'b')
if isfield(parameters,'title')
    title(parameters.title)
end
xlabel('Stage')
ylabel('N_{max}')

figure
hold all
plot(2:n_max,max_sim(2:n_max)./ntot_sim(2:n_max),'b')
xlabel('Stage')
ylabel('N_{max}/N_{tot}')



figure
hold all
fact=1.0;
for n=1:n_max
    stage=CELL_STAGE(n);
    %if n==1
    %    fact=10.0/max(stage.persize);
    %end

    plot(stage.time,fact*stage.rescaled_per,'k','LineWidth',1.5)


end
xlabel('Time')
ylabel('PER Size')
if isfield(parameters,'title')
    title(parameters.title)
end



end
