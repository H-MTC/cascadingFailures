function cfm = DATA_DC(casedata)

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
    mpc.bus(:, PD) = loads .* (0.8 + 0.2 * rand(busNum, 1));
    mpc = rundcopf(mpc, mpoption('verbose',0,'out.all',0));
    %   contingencies analysis for a power instance.
    for k = 1 : K + 1
        selectedBranches{k, 1} = nchoosek(1:branchNum, k - 1)';
        for con = 1 : size(selectedBranches{k, 1}, 2)
            mpc.branch(:, BR_STATUS) = 1;
            mpc.branch(selectedBranches{k, 1}(:, con), BR_STATUS) = 0;
            
            
            fras = redundancyExclusion(mpc);
            mpc.branch(:, RATE_A) = mpc.branch(:, RATE_A) * fras.capRatio;
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
                        cfi = CF_DC(mpc);
                        cols = size(cfi, 2);
                        cfm(:, samplesCount + 1 : samplesCount + cols) = cfi(remainingEntries, :);
                        samplesCount = samplesCount + cols;
%                         overFlowCount = overFlowCount + cols;
                        %                     size(cfm)
                        %                     toc
                        
                    end
                else
                    cfi = CF_DC(mpc);
                    cols = size(cfi, 2);
                    cfm(:, samplesCount + 1 : samplesCount + cols) = cfi(remainingEntries, :);
                    
                    samplesCount = samplesCount + cols;
                    %  samplesCount = samplesCount + size(cfm, 2);
                end
            end
            
        end
    end
    
end

cfm = cfm(:, 1 : samplesCount);
save('DATA_DCS.mat','cfm');
end

