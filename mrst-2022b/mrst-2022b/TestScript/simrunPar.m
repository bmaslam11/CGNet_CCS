function t_par=simrunPar()
%parallel
    nEnsemble=12;
    tic
    parfor i=1:nEnsemble
        fname=['ens-',num2str(i),'.txt'];
        CGsim(fname);
    end
    t_par=toc;
end


