function [perf, minv, ind] = acpfs(casedata, network, bus, points)
%%画出最小误差和平均误差曲线比较
define_constants;
mpce = loadcase(casedata);
branchNum = size(mpce.branch, 1);
busNum = size(mpce.bus, 1);

temp = mpce; 
temp.bus(:,[PD QD]) = 0;temp.gen(:,[PG QG]) = temp.gen(:,[PMAX QMAX]);
Submax = makeSbus(temp.baseMVA, temp.bus, temp.gen) * temp.baseMVA;
temp = mpce;temp.gen(:,[PG QG]) = temp.gen(:,[PMIN QMIN]);
Submin = makeSbus(temp.baseMVA, temp.bus, temp.gen) * temp.baseMVA;

lb = [real(Submin);imag(Submin)];
ub = [real(Submax);imag(Submax)];
emptyInjBuses = intersect(find(lb==0),find(ub==0));
remainingEntries = setdiff(1 : 2*busNum + branchNum, emptyInjBuses);







net = load(network);
tr = net.tr;
loads = mpce.bus(:, PD);
opt =  mpoption('out.all', 0, 'verbose', 0);

bfs = zeros(branchNum, points);
predictedBfs = bfs;
mpce.branch(:, BR_STATUS) = ones(branchNum, 1);
brokenIndices = [1 10];

mpce.bus(:, PD) = mpce.bus(:, PD) .* (0.6 + 0.4*rand(busNum, 1));
mpce.bus(:, QD) = mpce.bus(:, QD) .* (0.6 + 0.4*rand(busNum, 1));
rate_a=mpce.branch(:,RATE_A);
mpce.branch(:,RATE_A)=0;
loadsV = zeros(points,1);
for pt = 1 : points
    mpc = mpce;
    loadsV(pt) = loads(bus) * (0.6+  pt / points * 0.4);
    mpc.bus(bus, PD) = loads(bus) * (0.6+  pt / points * 0.4);
    mpc = runopf(mpc, opt);
    mpc.branch(brokenIndices, BR_STATUS) = 0;
    mpce.branch(:,RATE_A) = rate_a;
%     results = rundcpf(mpc);
%     bfs(:, pt) = results.branch(:, 14);
    results = CF_AC(mpc);
    bfs(:, pt) = results(busNum*2 + branchNum + 1: busNum*2 + branchNum*2, 1);
    subss = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen) * mpc.baseMVA;
    sampleO = [real(subss);imag(subss); mpc.branch(:, BR_STATUS)];
    sample = sampleO(remainingEntries);
    normalizedInputs = mapminmax.apply(sample,tr.settingsI);
    try
        predicted = predict(tr.net, normalizedInputs);
    catch
        predicted = sim(tr.net, normalizedInputs);
    end
    bfss = mapminmax.reverse(predicted, tr.settingsO);
    predictedBfs(:, pt) = bfss(1:branchNum);
end
ratedCap = repmat(mpc.branch(:, RATE_A), [1, points]);
absErr = abs((bfs - predictedBfs) ./ ratedCap);
% perf = sum(sum(absErr .* absErr)) / size(bfs, 1) / size(bfs, 2);
perf = sum(sum(absErr.*absErr)) / size(bfs, 1) / size(bfs, 2);


[minv, ind] = sort(sum(absErr,2));
% loadlevels=0.6+0.4/points*[1:points];
figure;
plot(loadsV, bfs(ind(1),:), 'r', loadsV, predictedBfs(ind(1), :),'b');
hold on;
plot(loadsV, bfs(ind(end),:), 'r', loadsV, predictedBfs(ind(end), :),'b');
% plot(1:points, predictedBfs(3, :));
% figure;
% plot(1:points, bfs)
end

