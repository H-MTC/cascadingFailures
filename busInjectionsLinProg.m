function [PGO, absDev] = busInjectionsLinProg(mpc, fras)
% RATE_A = 6;RATE_B = 7; RATE_C = 8; PD = 3; PMAX = 9; PMIN = 10;

define_constants;
branchNum = size(mpc.branch, 1);
busNum = size(mpc.bus, 1);
genNum = size(mpc.gen, 1);
%   compute active and reactive bus power injections

Aeq = ones(1, genNum);
beq = sum(mpc.bus(:, PD));
lb = mpc.gen(:, PMIN) * 0;
ub = mpc.gen(:, PMAX);
options = optimoptions('linprog','Display','off');
PGO = [];
absDev = [];
% mpc.branch(:, RATE_A) = mpc.branch(:, RATE_A) * fras.capRatio;
if fras.isLBon
    f = - fras.AL;
    A = [-f;f];
    b = [
        (mpc.branch(fras.branchIdx, RATE_A) - fras.flowOfLoad + fras.bounds)';
        (fras.flowOfLoad - mpc.branch(fras.branchIdx, RATE_A) + fras.bounds)'
        ];
    %%  range of power flows on lines may not cover its rated capacity, so the closest power flows need to be computed.
%     [pg,fval,flag.flagL, out] = quadprog(H,f,[],[],Aeq,beq,lb,ub,lb,[],options);
%     [pg,fval,flag.flagL, out] = linprog(f,A,b,Aeq,beq,lb,ub,lb,options);
    [pg,fval,flag.flagL, out] = linprog(f,A,b,Aeq,beq,lb,ub,lb,options);

    if flag.flagL > 0
        PGO(:,1) = pg;
        absDev(:,1) = abs(abs(f * pg) - (fras.flowOfLoad - mpc.branch(fras.branchIdx, RATE_A) + fras.bounds));
    end
end
if fras.isUBon
    f = fras.AR;
    A = [f;-f];
    b = [
        (fras.flowOfLoad + mpc.branch(fras.branchIdx, RATE_A) + fras.bounds)';
        (-fras.flowOfLoad - mpc.branch(fras.branchIdx, RATE_A) + fras.bounds)'
        ];
    %%  range of power flows on lines may not cover its rated capacity, so the closest power flows need to be computed.
%     [pg, fval, flag.flagR, out] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
    [pg,fval,flag.flagR, out] = linprog(f,A,b,Aeq,beq,lb,ub,lb,options);
    if flag.flagR > 0
    PGO = [PGO, pg];
    absDev = [absDev abs(f * pg - (fras.flowOfLoad + mpc.branch(fras.branchIdx, RATE_A) + fras.bounds)')];
    end
end



end


