function re = redundantExclusion(mpc)

%   find out the lines that may have over flows in for this mpc case whose correspongding grid is connected.
define_constants;

cfc = ext2int(mpc);
%   constraint of power balance
totalLoad = sum(cfc.bus(:, PD));
totalMaxGeneration = sum(cfc.gen(:, PMAX));
isPowerBalance = 0 <=totalLoad && totalLoad <= totalMaxGeneration;

re.capRatio = 1.0;
%   compute the maximum flow on each line and check whether the rated
%   capacity is violated.
PTDF = zeros(size(mpc.branch, 1), size(mpc.bus, 1));
if isPowerBalance
    branchNum = size(mpc.branch, 1);
    genNum = size(mpc.gen, 1);

    
    PTDF(cfc.order.branch.status.on, cfc.order.bus.status.on) = makePTDF(cfc);
    
    
    PTDF_T = [PTDF;-PTDF];
    PTDF_U = PTDF_T(:, mpc.gen(:, GEN_BUS));
    
    
    [sPTDF_U, idx] = sort(PTDF_U, 2, 'descend');

    temp = repmat(mpc.gen(:, PMAX), [1, branchNum * 2]);
    
    tdx = (0 : branchNum * 2 - 1)' * genNum;
    
    selection = temp((repmat(tdx, [1, genNum]) + idx)')';
    
    oper = triu(ones(genNum, genNum));
    totalGenRank = selection * oper;
    
    sufficientGenIdxL = sum(totalGenRank <= totalLoad, 2);
    sufficientGenIdxR = sum(~(totalLoad < totalGenRank), 2);
    
    sol = (selection .* (totalGenRank <= totalLoad))';
    totalGenRank = totalGenRank';
    if sufficientGenIdxL > 0
        sol(sufficientGenIdxL + tdx + 1) = totalLoad - totalGenRank(sufficientGenIdxL + tdx);
    end
    maximumFlows = zeros(branchNum * 2, 1);
    loadsPart = PTDF(:, mpc.bus(:, BUS_I)) * mpc.bus(:, PD);
    for br = 1 : branchNum
        maximumFlows(br) = sPTDF_U(br, :) * sol(:, br);
        maximumFlows(br + branchNum) = sPTDF_U(br + branchNum, :) * sol(:, br + branchNum);
    end
    
    
    positiveOverFlows = maximumFlows(1 : branchNum, 1) > mpc.branch(:, RATE_A) * re.capRatio + loadsPart;
    negativeOverFlows = maximumFlows(1 + branchNum : 2 * branchNum) > mpc.branch(:, RATE_A) * re.capRatio - loadsPart;
    fragileLines = find(or(positiveOverFlows, negativeOverFlows) == 1);
    
    % corrRowsInPTDFU = PTDF_U([fragileLines fragileLines + branchNum], :);
    % flowOfLoads = loadsPart(fragileLines);
    re.fragileLines = mpc.branch(fragileLines, end);
    re.corrRowsInPTDFUL = PTDF_U(fragileLines + branchNum, :);
    re.corrRowsInPTDFUR = PTDF_U(fragileLines, :);
    re.flowOfLoads = loadsPart(fragileLines);
    re.corrRowsInPTDF = PTDF_T([fragileLines fragileLines + branchNum], :);
    re.isLB = negativeOverFlows;
    re.isUB = positiveOverFlows;
    re.maxR = maximumFlows(1:branchNum);
    re.maxL = maximumFlows(1+branchNum:end);
else
    error('违反了功率平衡约束。。。');
end
end

