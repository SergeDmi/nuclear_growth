function [  ] = save_cell_stage(CELL_STAGE, parameters, exp_data)
% Saves results from nuclear groth sequence
%   CELL_STAGE is a structure containing nusize nucdata errdata persize maxize
%   N_tot P_tot time flag rescaled_per
%

if isfield(parameters,'folder')
  folder=parameters.folder;
else
  folder='results/';
end

if isfield(parameters,'prefix')
  prefix=parameters.prefix;
else
  prefix='nucleus_stage';
end

prefix=[folder prefix];

n_div=numel(CELL_STAGE);

max_sim =zeros(1,n_div);
ntot_sim=zeros(1,n_div);
maxsize =zeros(1,n_div);
errsize =zeros(1,n_div);

simul(n_div).res=[];

fname=[prefix '_simul.txt'];
fid = fopen(fname,'wt');
fprintf(fid,'# time \t nuc_size(exp.) \t err_nuc_size(exp.) \t nuc_size(sim) \t per_size(sim) \n');
fclose(fid);

for i=1:n_div
    simul(i).res=[CELL_STAGE(i).time CELL_STAGE(i).nucdata CELL_STAGE(i).errdata CELL_STAGE(i).nucsize CELL_STAGE(i).rescaled_per];


    dlmwrite(fname,simul(i).res,'delimiter',' ','-append')

    max_sim(i)=max(CELL_STAGE(i).nucsize);
    ntot_sim(i)=CELL_STAGE(i).N_tot;
    maxsize(i)=max(exp_data(i).nucdata);
    errsize(i)=max(exp_data(i).errdata);
end




max_sizes=[(1:n_div)' maxsize' errsize' max_sim' ntot_sim'];
fname=[prefix '_max_sizes.txt'];
fid = fopen(fname,'wt');
fprintf(fid,'# stage \t maxsize(exp) \t error(exp) \t maxsize(sim) \t max_possible_size(sim) \n');
fclose(fid);
dlmwrite(fname,max_sizes,'delimiter',' ','-append');




end
