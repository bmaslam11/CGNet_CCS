nEnsemble = 20;

%% sequential
tic
for i=1:nEnsemble
    fname=['ens-',num2str(i),'.txt'];
    CGsim(fname);
end
t=toc;

%% parallel
tic
parfor i=1:nEnsemble
    fname=['ens-',num2str(i),'.txt'];
    CGsim(fname);
end
t_par=toc;