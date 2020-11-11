function cfm = DATA_AC(casedata)

define_constants;
K = 2;
mpc = loadcase(casedata);
busNum = size(mpc.bus, 1);
branchNum = size(mpc.branch, 1);
genNum = size(mpc.gen, 1);
selectedBranches = cell(K, 1);

busesWithoutGen = setdiff(1:24,union(mpc.gen(:, GEN_BUS),[]));
busesWithoutLoads = find(mpc.bus(:,PD)==0);
remainingEntries = setdiff(1 : busNum + 2 * branchNum, intersect(busesWithoutGen,busesWithoutLoads));

%   dc opf for intact mpc


samplesCount = 0;
% overFlowCount = 0;
numOfBaseSamples = 1;
cfm = zeros(length(remainingEntries), 5000 * numOfBaseSamples);



loads= mpc.bus(:, PD);

for ls = 1 : numOfBaseSamples
    tic
    mpc.bus(:, PD) = loads .* (0.8 + 0.2 * rand(busNum, 1));
    mpc = runopf(mpc, mpoption('verbose',0,'out.all',0));
    %   contingencies analysis for a power instance.
    for k = 1 : K + 1
        selectedBranches{k, 1} = nchoosek(1:branchNum, k - 1)';
        for con = 1 : size(selectedBranches{k, 1}, 2)
            mpc.branch(:, BR_STATUS) = 1;
            mpc.branch(selectedBranches{k, 1}(:, con), BR_STATUS) = 0;
            
            
            fras = redundancyExclusion(mpc);
            mpc.branch(:, RATE_A) = mpc.branch(:, RATE_A) * fras.capRatio;
            mpcb = mpc;
            mpct = mpc;
            for fl = 1 : length(fras.fragileLines)
                
                fragile.maxL = fras.maxL(fras.fragileLines(fl));
                fragile.maxR = fras.maxR(fras.fragileLines(fl));
                fragile.AR = fras.corrRowsInPTDFUR(fl,:);
                fragile.AL = fras.corrRowsInPTDFUL(fl,:);
                fragile.branchIdx = fras.fragileLines(fl);
                fragile.isLBon = fras.isLB(fras.fragileLines(fl));
                fragile.isUBon = fras.isUB(fras.fragileLines(fl));
                fragile.flowOfLoad = fras.flowOfLoads(fl);
                fragile.bounds = 1e-1;    %   initial deviation bounds away the rated flow capacity of line.
                fragile.capRatio = fras.capRatio;
                
                [PGO, absDev] = busInjectionsLinProg(mpc, fragile);
                if ~isempty(PGO)
                    for ins = 1 : size(PGO, 2)
                        %                     tic
                        mpc.gen(:, PG) = PGO(:, ins);
                        mpcb.bus(:, [PD QD]) = mpc.bus(:, [PD QD]) *0.5;
                        mpcb.gen(:, [PG QG]) = mpc.gen(:, [PG QG]) *0.5;
                        mpct.bus(:, [PD QD]) = mpc.bus(:, [PD QD]) *1.0;
                        mpct.gen(:, [PG QG]) = mpc.gen(:, [PG QG]) *1.0; 
                        mpopt = mpoption('out.all', 0, 'verbose', 1);

                        %   remember to modify the function :
                        %   cpf_flim_event_cb.m and cpf_flim_event.m to
                        %   enforce active(real) branch flow limits
                        mpopt = mpoption(mpopt, 'cpf.enforce_flow_lims',1);

                        results = runcpf(mpcb, mpct, mpopt);
                        cfi = CF_AC(results);
                        cols = size(cfi, 2);
                        cfm(:, samplesCount + 1 : samplesCount + cols) = cfi(remainingEntries, :);
                        samplesCount = samplesCount + cols;
%                         overFlowCount = overFlowCount + cols;
                        %                     size(cfm)
                        %                     toc
                        
                    end
                else
                    mpcb.bus(:, [PD QD]) = mpc.bus(:, [PD QD]) *0.6;
                    mpcb.gen(:, [PG QG]) = mpc.gen(:, [PG QG]) *0.6;
                    mpct.bus(:, [PD QD]) = mpc.bus(:, [PD QD]) *1.2;
                    mpct.gen(:, [PG QG]) = mpc.gen(:, [PG QG]) *1.2; 
                    mpopt = mpoption('out.all', 0, 'verbose', 0);
%                         mpopt = mpoption(mpopt, 'cpf.step', 0.2);
                    mpopt = mpoption(mpopt, 'cpf.enforce_flow_lims',1);

                    results = runcpf(mpcb, mpct, mpopt);
                    cfi = CF_AC(results);
                    cols = size(cfi, 2);
                    cfm(:, samplesCount + 1 : samplesCount + cols) = cfi(remainingEntries, :);
                    
                    samplesCount = samplesCount + cols;
                    %  samplesCount = samplesCount + size(cfm, 2);
                end
            end
            
        end
    end
    toc
end

cfm = cfm(:, 1 : samplesCount);
save('DATA_AC.mat','cfm');

end

