function fras = redundancyExclusion(mpc)
%   Pick out the fragile lines, considering the situations that the grid may have isolated islands.

define_constants;


fras = [];
gra=graph(mpc.branch(logical(mpc.branch(:,BR_STATUS)), F_BUS), mpc.branch(logical(mpc.branch(:,BR_STATUS)), T_BUS));
[whichComp, comSize]=conncomp(gra);

nodeIdx = 1 : size(mpc.bus, 1);
for co = 1 : length(comSize)
    if comSize(co) > 1
        cfc = mpc;
        branchNum = size(mpc.branch, 1);
        cfc.branch(:, end + 1) = (1:branchNum)';
        
        
        fras.fragileLines = [];
        fras.corrRowsInPTDFUL = [];
        fras.corrRowsInPTDFUR = [];
        fras.flowOfLoads = [];
        fras.corrRowsInPTDF = [];
        fras.isLB = zeros(branchNum, 1);
        fras.isUB = zeros(branchNum, 1);
        fras.maxR = [];
        fras.maxL = [];


        nodesInThisComp = nodeIdx(whichComp == co);
        remoteBuses = setdiff(nodeIdx, nodesInThisComp);
        cfc.bus(remoteBuses, BUS_TYPE) = 4;
        if ~any(cfc.bus(:, BUS_TYPE) == 3)
            continue;
        end
        cfc.gen(ismember(cfc.gen(:, GEN_BUS), remoteBuses), GEN_STATUS) = 0;
        subG = subgraph(gra, nodesInThisComp);
        subEdges = table2array(subG.Edges);
        for tr = 1 : size(subEdges, 1) * size(subEdges, 2)
            subEdges(tr) = nodesInThisComp(subEdges(tr));
        end
        branchIdxInThisSubGraph  = ismember(mpc.branch(:,[F_BUS T_BUS]), subEdges,'rows');
        cfc.branch(:, BR_STATUS) = cfc.branch(:, BR_STATUS) & branchIdxInThisSubGraph;
        re = redundantExclusion(cfc);
        
        fras.fragileLines = [fras.fragileLines;re.fragileLines];
        fras.corrRowsInPTDFUL = [fras.corrRowsInPTDFUL; re.corrRowsInPTDFUL];
        fras.corrRowsInPTDFUR = [fras.corrRowsInPTDFUR; re.corrRowsInPTDFUR];
        fras.flowOfLoads = [fras.flowOfLoads; re.flowOfLoads];
        fras.corrRowsInPTDF = [fras.corrRowsInPTDF; re.corrRowsInPTDF];
        fras.isLB = fras.isLB | re.isLB;
        fras.isUB = fras.isUB | re.isUB;
        fras.maxR = [fras.maxR; re.maxR];
        fras.maxL = [fras.maxL; re.maxL];
        fras.capRatio = re.capRatio;
        
    end
end


end

