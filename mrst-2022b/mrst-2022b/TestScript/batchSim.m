c = parcluster();
pathWorking=pwd;

job=batch(c,@simrunPar,1,{}, 'Pool',4,'CurrentFolder',pathWorking,'AdditionalPaths',pathWorking);
%tElapsed=fetchOutputs(job)