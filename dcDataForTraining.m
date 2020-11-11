function ACCFM = acDataForTraining(caseData, numOfPowerPlans)
%   @Jeivier, 2020/06/03, Functionality : compute the ac power flow data for training the neural network.

%   As the output of this function, the matrix ACCFM has 2 * busNum + 2 * branchNum + 1 rows, each row represent a single stage of cascading
%   failures of lines, which is represented by a function PL = f(PB, QB, LS). Please refer to the appended paper for more information.

%   initialization of variables.
define_constants; 
original = loadcase(caseData);
mpc=original;

busNum= size(mpc.bus);
branchNum = size(mpc.branch);

busNumWithGens = union(mpc.gen(:, GEN_BUS),[]);
busGenCapPG = zeros(length(busNumWithGens), 1);
busGenCapQG = zeros(length(busNumWithGens), 1);

% numOfPowerPlans = 5000;
ACCFM= zeros(busNum(1) * 2 + 2 * branchNum(1) + 1,numOfPowerPlans * 700);
powerInj = zeros(busNum(1) * 2, 1);
busPG = zeros(busNum(1), 1);
busQG = busPG;



count = 0;
for IP = 1 : numOfPowerPlans
    
    

isGood = false;
re_opf = [];
% warnMsg = [];
while ~isGood
    mpc_opf = original;
    mpc_opf.branch(:,RATE_A) = 0;%  make power flows of lines in opf unlimited.
    
%   according to the data of shannxi province, the ratio of maximum real
%   load over minimum real load is 1.3635, the ratio of average real load
%   over minimum real load is 1.0888.
    mpc_opf.bus(:,PD) = original.bus(:,PD) + 0.3635 * original.bus(:,PD) .* rand(busNum(1), 1);
    mpc_opf.bus(:,QD) = original.bus(:,QD) + 0.3635 * original.bus(:,QD) .* rand(busNum(1), 1);
    
%   convergence of ac algorithm may be not guaranteed for any given power
%   plan by Newton method. So, we need to find a power plan by runopf() to
%   make its result converge.
    warning('') % Clear last warning message
    [re_opf, isGood] = rundcopf(mpc_opf, mpoption('out.all',0));
    [warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg)
        isGood = false;
    end
end

%   finish the instance generation of loads and generators.
%     mpc.bus(:, [PD QD]) = re_opf.bus(:, [PD QD]);
%     mpc.gen(:, [PG QG]) = re_opf.gen(:, [PG QG]);
    
%   compute the real and reactive power of buses.
    for IBNWG = 1 : length(busNumWithGens)
        busGenCapPG(IBNWG) = sum(re_opf.gen(ismember(mpc.gen(:, GEN_BUS), busNumWithGens(IBNWG)), PG));
        busGenCapQG(IBNWG) = sum(re_opf.gen(ismember(mpc.gen(:, GEN_BUS), busNumWithGens(IBNWG)), QG));
    end
    busPG(busNumWithGens) = busGenCapPG;
    busQG(busNumWithGens) = busGenCapQG;
    powerInj(1 : busNum, 1) = -re_opf.bus(:, PD) + busPG;
    powerInj(busNum + 1 : end, 1) = -re_opf.bus(:, QD) + busQG;
    

%   for each instance of bus power injection, compute the processes of
%   cascading failures triggered by each contingency of lines.
    for k = 0:2 % select 1-2 outaged lines from all lines randomly.
        marker = zeros(branchNum(1), 1);
        marker(1:k) = 1;
        while any( marker(1 : branchNum(1) - k))
            for i = 1 : branchNum(1) - 1
                if(marker(i) == 1 && marker(i+1) == 0)
                    marker(i) = 0;
                    marker(i+1) = 1;
                    numOfValue1 = sum(marker(1 : i - 1));
                    marker(1 : i - 1) = 0;
                    marker(1 : numOfValue1) = 1;
                    re_opf.branch(:, BR_STATUS) = ~marker;%    set the selected lines to be outaged.


                    %   give the flawed system to the ac power flow model
                    CF_stateSeries = CF_Evolution(re_opf);%    compute the cascading failure under a give power injection of buses and flawed status of lines.
                    sizeOf = size(CF_stateSeries);
                    if ~isempty(CF_stateSeries)
                        ACCFM(1 : busNum(1) * 2, count + 1 : count + sizeOf(2)) = repmat(powerInj, [1, sizeOf(2)]);
                        ACCFM(busNum(1) * 2 + 1 : end, count + 1: count + sizeOf(2)) = CF_stateSeries;
                        count = count + sizeOf(2);
                    end
                    break;
                end
            end
        end
    end
end
% save('DATA_ACCFM.mat','ACCFM');
ACCFM = ACCFM(:, 1:find(ACCFM(1,:)~=0,1,'last'));
index=find(sum(ACCFM(1:48,1:100),2)==0);
ACCFM(index,:)=[];
% save('DATA_DCCFM.mat','ACCFM');
end

function CF_stateSeries = CF_Evolution(mpc)
define_constants; 
%   基于交流潮流的连锁故障演化程序 @MTC 20191018
%   考虑满足：功率平衡约束(包含系统解列后如果孤岛电网仍然满足功率平衡的情况)。
%   输入节点注入功率向量和线路状态0/1向量，输出连锁故障演化序列。
%   it has been verified that the ac power flow algorithm in matpower6.0
%   can have overloaded lines when power plans goes big or contingency
%   happens. Modificcatins to power plans and lines can be carried out on
%   standard casefile or its struct variable.

%   assumptions: 
%   1) contingencies of lines still preserve the connectivity of the grid; 
%   2) the power flows computed by ac algorithm is based on the standard
%   circumstances, that is, 50Hz frequency of voltages and currents. So,
%   when the cnstraint of power balance is violated under some specific
%   contingency, the quantity of the violations of power balance will be
%   recorded for further purposes. For example, if power generation is
%   greater than demand, the frequency of the system may inncrease, leading
%   to greater power flows on lines.

busNum = size(mpc.bus);
branchNum = size(mpc.branch);
%   check the connectivity of the flawed grid.
G = graph(mpc.branch(logical(mpc.branch(:, BR_STATUS)),F_BUS),mpc.branch(logical(mpc.branch(:, BR_STATUS)),T_BUS));
[comIndex, nodeNumInCom] = conncomp(G);
isConnected = busNum(1) == max(nodeNumInCom);


iterCount = 1;
CF_stateSeries = zeros(2 *  branchNum(1) + 1, branchNum(1));
while isConnected
    CF_stateSeries(1 : branchNum, iterCount) = ~mpc.branch(:, BR_STATUS);
    warning('') % Clear last warning message
    [results, isSucess] = rundcpf(mpc, mpoption('out.all',0));
    [warnMsg, warnId] = lastwarn;
    PMEAN = (results.gen(:,PMAX) + results.gen(:,PMIN)) / 2;
    PHALFRANGE = (results.gen(:,PMAX) - results.gen(:,PMIN)) / 2;
    QMEAN = (results.gen(:,QMAX) + results.gen(:,QMIN)) / 2;
    QHALFRANGE = (results.gen(:,QMAX) - results.gen(:,QMIN)) / 2;
    isFeasible = all(results.gen(:,PG) >= PMEAN - PHALFRANGE  * 1.1) && all(results.gen(:,PG) <= PMEAN + PHALFRANGE  * 1.1) && all(results.gen(:,QG) >= QMEAN - QHALFRANGE  * 1.1) && all(results.gen(:,QG) <= QMEAN + QHALFRANGE  * 1.1);
%     label = max(results.branch(:,14)) > 1e5;
    if ~isempty(warnMsg) || ~isFeasible% || ~isSucess
        CF_stateSeries = [];
        break;
    end
    outagedLines = results.branch(abs(results.branch(:, 14)) > results.branch(:, RATE_B), [F_BUS T_BUS]);
    brokenEdgeNums = find(ismember(mpc.branch(:, [F_BUS T_BUS]), outagedLines, 'rows'));
    mpc.branch(brokenEdgeNums, BR_STATUS) = 0;
    G = graph(mpc.branch(logical(BR_STATUS),F_BUS),mpc.branch(logical(BR_STATUS),T_BUS));
    [comIndex, nodeNumInCom] = conncomp(G);
    isConnected = busNum(1) == max(nodeNumInCom);


    CF_stateSeries(branchNum + 1 : 2 * branchNum, iterCount) = results.branch(:, 14);
    deviation = sum(results.gen(:,PG)) - sum(results.bus(:,PD))-sum(abs(results.branch(:,14)+results.branch(:,16)));
    CF_stateSeries(end, iterCount) = deviation; %delta of real power
    if isempty(brokenEdgeNums) || ~isSucess
        break;
    end
    iterCount = iterCount + 1;
    

end
evoStages = sum(sum(CF_stateSeries, 1) ~= 0);
CF_stateSeries = CF_stateSeries(:,1:evoStages);

end
