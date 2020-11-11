function cfm = DATA_ACN(casedata, bt)

define_constants;
K = 1;
% mpc = loadcase(casedata);
mpc=load(casedata);
mpc=mpc.mpc;
busNum = size(mpc.bus, 1);
branchNum = size(mpc.branch, 1);
genNum = size(mpc.gen, 1);
selectedBranches = cell(K, 1);





% x1=mpc.gen(:,10);y1=mpc.gencost(:,5).*x1.^2+mpc.gencost(:,6).*x1+mpc.gencost(:,7);
% x2=mpc.gen(:,9);y2=mpc.gencost(:,5).*x2.^2+mpc.gencost(:,6).*x2+mpc.gencost(:,7);
% obj=(y1-y2)./(x1-x2);
obj = mpc.gencost(:,5);
%   dc opf for intact mpc

temp = mpc; 
temp.bus(:,[PD QD]) = 0;temp.gen(:,[PG QG]) = temp.gen(:,[PMAX QMAX]);
Submax = makeSbus(temp.baseMVA, temp.bus, temp.gen) * temp.baseMVA;
temp = mpc;temp.gen(:,[PG QG]) = temp.gen(:,[PMIN QMIN]);
Submin = makeSbus(temp.baseMVA, temp.bus, temp.gen) * temp.baseMVA;

lb = [real(Submin);imag(Submin)];
ub = [real(Submax);imag(Submax)];
emptyInjBuses = intersect(find(lb==0),find(ub==0));
remainingEntries = setdiff(1 : 2*busNum + 2 * branchNum, emptyInjBuses);
remainingBuses = setdiff(1:2*busNum,emptyInjBuses);

lb = lb(remainingBuses);
ub = ub(remainingBuses);
samplesCount = 0;
% overFlowCount = 0;
numOfBaseSamples = 2800;
cfm = zeros(length(remainingEntries), 200 * numOfBaseSamples);




loads = mpc.bus(:, PD);
for ls = 1 : numOfBaseSamples
    mpc.bus(:, PD) = loads .* rand(busNum, 1);
    mpc.bus(:, QD) = loads .* rand(busNum, 1);
    mpc = rundcopf(mpc, mpoption('verbose',0,'out.all',0));
    %   contingencies analysis for a power instance.
    for k = 1 : K + 1
        selectedBranches{k, 1} = nchoosek(1:branchNum, k - 1)';
        for con = 1 : size(selectedBranches{k, 1}, 2)
            mpc.branch(:, BR_STATUS) = 1;
            mpc.branch(selectedBranches{k, 1}(:, con), BR_STATUS) = 0;

            cfi = CF_AC(mpc);
            cols = size(cfi, 2);
            cfi = cfi(:,sum(cfi(remainingBuses,:)>repmat(lb,[1, size(cfi,2)])-5e1  & cfi(remainingBuses,:)<repmat(ub,[1, size(cfi,2)])+5e1) == length(remainingBuses));
            if ~isempty(cfi)
                cfm(:, samplesCount + 1 : samplesCount + cols) = cfi(remainingEntries, :);
                samplesCount = samplesCount + cols;
            end
            
        end
    end
    
end

cfm = cfm(:, 1 : samplesCount);
save(['DATA_ACN_CASE118N',num2str(bt),'.mat'],'cfm');
end

function pg = miniLinProg(mpc,obj)
% RATE_A = 6;RATE_B = 7; RATE_C = 8; PD = 3; PMAX = 9; PMIN = 10;

define_constants;
% branchNum = size(mpc.branch, 1);
% busNum = size(mpc.bus, 1);
genNum = size(mpc.gen, 1);
%   compute active and reactive bus power injections

Aeq = ones(1, genNum);
beq = sum(mpc.bus(:, PD));
lb = mpc.gen(:, PMIN) * 0;
ub = mpc.gen(:, PMAX);
options = optimoptions('linprog','Display','off');

[pg,fval,flag.flagL, out] = linprog(obj,[],[],Aeq,beq,lb,ub,lb,options);



end


