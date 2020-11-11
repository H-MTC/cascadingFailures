function cfm = CF_DC(mpc)
define_constants;
branchNum = size(mpc.branch,1);
busNum = size(mpc.bus, 1);
statusOfLinesPre = ones(branchNum, 1);

cfm = zeros(busNum * 2 + 2 * branchNum, branchNum);
count = 0;
balanceBuses = union(find(mpc.bus(:, PD) ~= 0), mpc.gen(:, GEN_BUS));

while any(statusOfLinesPre - mpc.branch(:, BR_STATUS))
    %%  code area for cases not considering the islanding scenarios.
    statusOfLinesPre = mpc.branch(:, BR_STATUS);
    gra=graph(mpc.branch(logical(mpc.branch(:,BR_STATUS)), F_BUS), mpc.branch(logical(mpc.branch(:,BR_STATUS)), T_BUS));
    digra=digraph(mpc.branch(logical(mpc.branch(:,BR_STATUS)), F_BUS), mpc.branch(logical(mpc.branch(:,BR_STATUS)), T_BUS));
    [whichComp, comSize]=conncomp(gra);
    
    nodeIdx = 1 : busNum;
    branchFlows = zeros(branchNum,1);
    cfc = mpc;
    [maxNN, maxNNidx] = max(comSize);
    nodesInThisComp = nodeIdx(whichComp == maxNNidx);

    isPowerBalance = all(ismember(balanceBuses, nodesInThisComp));
    
    
    
    
    if isPowerBalance
        remoteBuses = setdiff(nodeIdx, nodesInThisComp);
        remoteGenBusIdx = ismember(cfc.gen(:, GEN_BUS), remoteBuses);
        cfc.bus(remoteBuses, BUS_TYPE) = 4;
        cfc.gen(remoteGenBusIdx, GEN_STATUS) = 0;
        if ~any(cfc.bus(nodesInThisComp, BUS_TYPE) == 3)
            [m, imax] = max(cfc.bus(nodesInThisComp, PD));
            cfc.bus(nodesInThisComp(imax), BUS_TYPE) = 3;
        end

        subG = subgraph(digra, nodesInThisComp);
        subEdges = table2array(subG.Edges);
        for tr = 1 : size(subEdges, 1) * size(subEdges, 2)
            subEdges(tr) = nodesInThisComp(subEdges(tr));
        end
        branchIdxInThisSubGraph  = ismember(mpc.branch(:,[F_BUS T_BUS]), subEdges,'rows');
        cfc.branch(:, BR_STATUS) = cfc.branch(:, BR_STATUS) & branchIdxInThisSubGraph;
        cfc = runpf(cfc, mpoption('verbose', 0,'out.all',0));
        branchFlows(branchIdxInThisSubGraph) = cfc.branch(branchIdxInThisSubGraph, 14);
        mpc.branch(:, BR_STATUS) = statusOfLinesPre & abs(branchFlows) < mpc.branch(:, RATE_A);
        Subs = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen) * mpc.baseMVA;
%         sample = [real(Subs);imag(Subs);statusOfLinesPre;branchFlows;cfc.bus(:, VM);cfc.bus(:, VA)];
        sample = [real(Subs);imag(Subs);statusOfLinesPre;branchFlows];
        count = count + 1;
        cfm(:, count) = sample;
    else
        break;
    end

end

cfm = cfm(:, 1 : count);
end

