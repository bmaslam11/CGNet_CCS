function obj = matchObservedGW(model, states, schedule, observed, varargin)
% Compute mismatch-function edited for gas (CO2) water system

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
opt     = struct('WaterRateWeight',     [] , ...
                 'GasRateWeight',       [] , ...
                 'BHPWeight',           [] , ...
                 'ComputePartials',     false, ...
                 'tStep' ,              [], ...
                 'state',               [],...
                 'from_states',         true,...% can be false for generic models
                 'matchOnlyProducers',  false);
             
opt     = merge_options(opt, varargin{:});

dts   = schedule.step.val;
totTime = sum(dts);

tSteps = opt.tStep;
if isempty(tSteps) %do all
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    numSteps = 1;
    dts = dts(opt.tStep);
end


obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
%   sol_obs = observed{tSteps(step)};
    sol_obs = observed{tSteps(step)};
    nw      = numel(sol_obs.wellSol);
    if opt.matchOnlyProducers
        matchCases = (vertcat(sol.sign) < 0);
    else
        matchCases = true(nw,1);
    end
    qWs_obs = vertcatIfPresent(sol_obs.wellSol, 'qWs', nw);
    qGs_obs = vertcatIfPresent(sol_obs.wellSol, 'qGs', nw);
    bhp_obs = vertcatIfPresent(sol_obs.wellSol, 'bhp', nw);
    status_obs = vertcat(sol_obs.wellSol.status);
    
    [ww, wo, wp] = getWeights(qWs_obs, qGs_obs, bhp_obs, opt);
    
    if opt.ComputePartials
        if(opt.from_states)
            init=true;
            state = model.getStateAD( states{tSteps(step)}, init);
        else
            state = opt.state;
        end
        qWs = model.FacilityModel.getProp(state,'qWs');
        qGs = model.FacilityModel.getProp(state,'qGs');
        bhp = model.FacilityModel.getProp(state,'bhp');
        assert(not(isnumeric(qWs))); 
        status = vertcat(state.wellSol.status);
     else
        state = states{tSteps(step)};
        [qWs, qGs, bhp] = deal( vertcatIfPresent(state.wellSol, 'qWs', nw), ...
                                vertcatIfPresent(state.wellSol, 'qGs', nw), ...
                                vertcatIfPresent(state.wellSol, 'bhp', nw) );
       assert(isnumeric(qWs));
       status = vertcat(state.wellSol.status);
    end

    if all(status == status_obs) 
        if any(~status)
            matchCases = matchCases(status);
        end
    else % problematic status ignore bhp, treat qWs, qGs as zero
        [bhp, bhp_obs] = expandToFull(bhp, bhp_obs, status, status_obs, true);
        [qWs, qWs_obs] = expandToFull(qWs, qWs_obs, status, status_obs, false);
        [qGs, qGs_obs] = expandToFull(qGs, qGs_obs, status, status_obs, false);
    end
    dt = dts(step);
    obj{step} = (dt/(totTime*nnz(matchCases)))*sum( ...
                                (ww*matchCases.*(qWs-qWs_obs)).^2 + ...
                                (wo*matchCases.*(qGs-qGs_obs)).^2 + ...
                                (wp*matchCases.*(bhp-bhp_obs)).^2 );
end
end

%--------------------------------------------------------------------------

function v = vertcatIfPresent(sol, fn, nw)
if isfield(sol, fn)
    v = vertcat(sol.(fn));
    assert(numel(v)==nw);
    v = v(vertcat(sol.status));
else
    v = zeros(nnz(sol.status),1);
end
end

%--------------------------------------------------------------------------

function [v, v_obs] = expandToFull(v, v_obs, status, status_obs, setToZero)
tmp = zeros(size(status));
if isa(v, 'ADI')
    tmp = double2ADI(tmp, v);
end
tmp(status) = v;
v = tmp;
%
tmp = zeros(size(status));
tmp(status_obs) = v_obs;
v_obs = tmp;
if setToZero
    ix = status ~= status_obs;
    v(ix)     = 0;
    v_obs(ix) = 0;
end

end
%--------------------------------------------------------------------------

function  [ww, wo, wp] = getWeights(qWs, qGs, bhp, opt)
ww = opt.WaterRateWeight;
wo = opt.GasRateWeight;
wp = opt.BHPWeight;

rw = sum(abs(qWs)+abs(qGs));

if isempty(ww)
    % set to zero if all are zero
    if sum(abs(qWs))==0
        ww = 0;
    else
        ww = 1/rw;
    end
end

if isempty(wo)
    % set to zero if all are zero
    if sum(abs(qGs))==0
        wo = 0;
    else
        wo = 1/rw;
    end
end

if isempty(wp)
    % set to zero all are same
    dp = max(bhp)-min(bhp);
    if dp == 0
        wp = 0;
    else
        wp = 1/dp;
    end
end
end

