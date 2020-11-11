function [initialPoint, optimalPoint, initialBranchFlows, branchFlows, afterBranchFlows] = gradDescentAC(casedata, overFlowData, network)

define_constants;
global contingency branchFlows tr mpc busNum branchNum genNum
branchFlows = [];
mpc = loadcase(casedata);
net = load(network);
tr = net.tr;
data = load(overFlowData);
cfm = data.ocfm;

selection = 1000;

initialPoint = cfm.cfm(1:cfm.busNum, selection);
initialBranchFlows = cfm.cfm(cfm.busNum + cfm.branchNum+1:end, selection);

temp = mpc; 
temp.bus(:,[PD QD]) = 0;temp.gen(:,[PG QG]) = temp.gen(:,[PMAX QMAX]);
Submax = makeSbus(temp.baseMVA, temp.bus, temp.gen) * temp.baseMVA;
temp = mpc;temp.gen(:,[PG QG]) = temp.gen(:,[PMIN QMIN]);
Submin = makeSbus(temp.baseMVA, temp.bus, temp.gen) * temp.baseMVA;

lbt = [real(Submin);imag(Submin)];
ubt = [real(Submax);imag(Submax)];
lb0 = find(lbt==0);
ub0 = find(ubt==0);
emptyInjBuses = intersect(lb0,ub0);
remainingBuses = setdiff(1:2*busNum,emptyInjBuses);



eps = 5e0;
% A = zeros(4, cfm.busNum);
% A(1, 1:cfm.busNum/2) = 1;
% A(2, 1:cfm.busNum/2) = -1;
% A(3, cfm.busNum/2 + 1:cfm.busNum) = 1;
% A(4, cfm.busNum/2 + 1:cfm.busNum) = -1;
% b = [eps;eps;eps;eps];

A = zeros(2, cfm.busNum);
A(1, 1:cfm.busNum/2) = 1;
A(2, 1:cfm.busNum/2) = -1;
b = [eps;eps];

lb = cfm.lb - eps;  %control the range of mitigations
ub = cfm.ub + eps;  %control the range of mitigations



lb=lbt(remainingBuses);
ub=ubt(remainingBuses);
% lb = lb - abs(initialPoint) * 0.1;  %control the range of mitigations
% ub = ub + abs(initialPoint) * 0.1;  %control the range of mitigations
contingency = logical(cfm.cfm(cfm.busNum+1:cfm.busNum + cfm.branchNum, selection));
mpc.branch(:, BR_STATUS) = contingency;

genNum = size(mpc.gen, 1);
busNum = size(mpc.bus, 1);
branchNum = size(mpc.branch, 1);

% opts = optimoptions('fmincon','FinDiffRelStep',1e-1,'ConstraintTolerance',5e-1);
opts=[];
optimalPoint = fmincon(@neuralFun,initialPoint,A,b,[],[],lb,ub, [], opts);

% lbr = zeros(2 * busNum + 2 * genNum, 1);
% ubr = [mpc.bus(:, PD);mpc.bus(:, QD);mpc.gen(:, PMAX);mpc.gen(:, QMAX)];
% 
% Aine = zeros(busNum * 2, busNum * 2 + genNum * 2);
% for ent = 1 : length(remainingBuses)
%     Aine(remainingBuses(ent), remainingBuses(ent)) = -1;
%     Aine(remainingBuses(ent), 48 + find([mpc.gen(:, GEN_BUS); 24 + mpc.gen(:, GEN_BUS)] == remainingBuses(ent))) = 1;
% end
% 
% b = zeros(busNum * 2, 1);
% b(remainingBuses) = optimalPoint;
% x0 = (lbr + ubr) / 2;
% optsM = optimoptions('fmincon','FinDiffRelStep',1e-1,'ConstraintTolerance',5e-1, 'MaxFunEvals', 3000);
% 
% optimalPlans = linprog(f,Aine,b,[],[],lbr,ubr,[], optsM);


%%  满足节点的负荷和发电加和后等于上述求解的optimalPoint，然后通过fmincon的目标函数runpf使得潮流越限量最小。
lbr = zeros(2 * busNum + 2 * genNum, 1);
ubr = [mpc.bus(:, PD);mpc.bus(:, QD);mpc.gen(:, PMAX);mpc.gen(:, QMAX)];

Aine = zeros(busNum * 2, busNum * 2 + genNum * 2);
for ent = 1 : length(remainingBuses)
    Aine(remainingBuses(ent), remainingBuses(ent)) = -1;
    Aine(remainingBuses(ent), 48 + find([mpc.gen(:, GEN_BUS); 24 + mpc.gen(:, GEN_BUS)] == remainingBuses(ent))) = 1;
end

b = zeros(busNum * 2, 1);
b(remainingBuses) = optimalPoint;
x0 = (lbr + ubr) / 2;
optsM = optimoptions('fmincon','FinDiffRelStep',1e-1,'ConstraintTolerance',5e-1, 'MaxFunEvals', 3000);
optimalPlans = fmincon(@minOverflow,x0,[],[],Aine,b,lbr,ubr,[], optsM);
mpc.bus(:, 3) = optimalPlans(1:busNum);
mpc.bus(:, 4) = optimalPlans(busNum + 1 : 2 * busNum);
mpc.gen(:, 2) = optimalPlans(2 * busNum + 1 : 2 * busNum + genNum);
mpc.gen(:, 3) = optimalPlans(2 * busNum + genNum + 1 : 2 * busNum + 2 * genNum);
 
mpopt = mpoption('verbose',0, 'out.all',0);
re = runpf(mpc, mpopt);
afterBranchFlows = re.branch(:, 14);

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

loss = double(sum(exp(abs(branchFlows(contingency)) ./ mpc.branch(contingency, 6) - 1)));

end

function obj = minOverflow(x)
global mpc busNum genNum contingency
mpc.bus(:, 3) = x(1:busNum);
mpc.bus(:, 4) = x(busNum + 1 : 2 * busNum);
mpc.gen(:, 2) = x(2 * busNum + 1 : 2 * busNum + genNum);
mpc.gen(:, 3) = x(2 * busNum + genNum + 1 : 2 * busNum + 2 * genNum);
mpopt = mpoption('verbose',0, 'out.all',0);
re = runpf(mpc, mpopt);
branchFlows = re.branch(:, 14);
% obj = double(sum(exp(abs(branchFlows(contingency)) ./ mpc.branch(contingency, 6) - 1)));
obj = double(sum(branchFlows(contingency).^2 / 2500));

end