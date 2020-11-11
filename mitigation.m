function [initialPoint, optimalPoint, initialBranchFlows, branchFlows, afterBranchFlows] = mitigation(casedata, overFlowData, network)

define_constants;
global contingency branchFlows tr mpc busNum branchNum
branchFlows = [];


%%  系统信息
mpc = loadcase(casedata);
genNum = size(mpc.gen, 1);
busNum = size(mpc.bus, 1);
branchNum = size(mpc.branch, 1);

%%  NNA及相关数据
net = load(network);
tr = net.tr;
data = load(overFlowData);
cfm = data.cfm;



%%  留待之后和被优化过的潮流作比较。
initialBranchFlows = cfm.cfm(cfm.busNum + cfm.branchNum+1:end, selection);


%%  指定需要缓解的故障场景及初始点
selection = 1;
contingency = logical(cfm.cfm(cfm.busNum+1:cfm.busNum + cfm.branchNum, selection));
mpc.branch(:, BR_STATUS) = contingency;

%%  初始点为节点发电量。
initialPoint = cfm.cfm(1:cfm.busNum, selection);


eps = 1e1;
%%  添加线性/非线性约束，使得节点有功和无功注入平衡，注意网损ref75,ref76,ref77。
%%  功率平衡可以由matpower-7.0 Page 33页的方程描述，可以添加为本模型的非线性约束。
A = zeros(2, cfm.busNum);
A(1, 1:cfm.busNum/2) = 1;
A(2, 1:cfm.busNum/2) = -1;
b = [eps;eps];


%%  节点功率注入范围约束，可由固定的负荷和发电机出力范围确定。
lb = initialPoint - abs(initialPoint) * 0.1;  %control the range of mitigations
ub = initialPoint + abs(initialPoint) * 0.1;  %control the range of mitigations


%%  优化参数设置
opts = optimoptions('fmincon','FinDiffRelStep',1e-1);
% opts=[];

%%  优化变量为节点有功和无功发电量
optimalPoint = fmincon(@neuralFun,initialPoint,A,b,[],[],lb,ub, [], opts);




%节点的负荷固定，发电机的出力为变量，利用线性规划求出一个发电计划。
% lbl =zeros(2*busNum+2*genNum, 1);
% ubl = [mpc.bus(:,PD);mpc.gen(:,PG);mpc.bus(:, QD);mpc.gen(:,QG)];
lbl =zeros(2*busNum+2*genNum, 1);
ubl = 1e3*ones(2*busNum+2*genNum, 1);
f = zeros(1, length(lbl));
Aeq = zeros(2*busNum, 2 * busNum + 2 * genNum);
Aeq(1:busNum, 1:busNum) = -eye(busNum);
Aeq(busNum + 1:2*busNum, busNum + genNum + 1:2*busNum + genNum) = -eye(busNum);
for i = 1 : busNum
    genidx = find(mpc.gen(:, 1) == i);
    Aeq(i, busNum + genidx) = 1;
    Aeq(i + busNum, 2*busNum + genNum + genidx) = 1;
end
beq = optimalPoint;

cplxpower = linprog(f, [], [], Aeq, beq, lbl, ubl);
mpc.bus(:,PD) = cplxpower(1:busNum);
mpc.gen(:,PG) = cplxpower(busNum+1:busNum+genNum);
mpc.bus(:,QD) = cplxpower(busNum+genNum+1:2*busNum+genNum);
mpc.gen(:,QG) = cplxpower(2*busNum+genNum+1:2*busNum+2*genNum);
sbus=makeSbus(mpc.baseMVA, mpc.bus, mpc.gen)*mpc.baseMVA;
opt = mpoption('out.all',0,'verbose',0);
re = runpf(mpc, opt);
sbus=makeSbus(re.baseMVA, re.bus, re.gen)*re.baseMVA;
afterBranchFlows = re.branch(:,14);

subplot(2,1,1);bar(1:branchNum, [initialBranchFlows(1:branchNum) branchFlows afterBranchFlows]);
subplot(2,1,2);bar(1:length(initialPoint), [initialPoint optimalPoint]);
% runpf不能从标幺值的电压初值迭代到电压远离标幺值的非正常潮流状态，这是牛顿法的固有缺陷。
end

function loss = neuralFun(x)
global contingency branchFlows tr mpc busNum branchNum
sample = [x; contingency];
normalizedInputs = mapminmax.apply(sample,tr.settingsI);
try
    predicted = predict(tr.net, normalizedInputs);
catch
    predicted = sim(tr.net, normalizedInputs);
end
prediction = double(mapminmax.reverse(predicted, tr.settingsO));
branchFlows=prediction(1:branchNum);
% branchFlows = mapminmax.reverse(predicted, tr.settingsO);
% prediction = double(mapminmax.reverse(predicted, tr.settingsO));
% V = prediction(end-2*busNum+1:end-busNum)+1i*prediction(end-busNum+1:end);
% [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
% [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
% [bus, gen, branch] = pfsoln(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, Ybus, Yf, Yt, V, ref, pv, pq);
% branchFlows = branch(:,14);

% losses = get_losses(mpc);

loss = double(sum(exp(abs(branchFlows(contingency)) ./ mpc.branch(contingency, 6) - 1)));

end

