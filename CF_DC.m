function cfm = CF_DC(mpc)
define_constants;
branchNum = size(mpc.branch,1);
busNum = size(mpc.bus, 1);
statusOfLinesPre = ones(branchNum, 1);

cfm = zeros(busNum + 2 * branchNum, branchNum);
count = 0;
while any(statusOfLinesPre - mpc.branch(:, BR_STATUS))
    %     cfc = mpc;
    
    
    statusOfLinesPre = mpc.branch(:, BR_STATUS);
    gra=graph(mpc.branch(logical(mpc.branch(:,BR_STATUS)), F_BUS), mpc.branch(logical(mpc.branch(:,BR_STATUS)), T_BUS));
    [whichComp, comSize]=conncomp(gra);
    
    nodeIdx = 1 : busNum;
    branchFlows = zeros(branchNum,1);
    for co = 1 : length(comSize)
        cfc = mpc;
        if comSize(co) > 1
            nodesInThisComp = nodeIdx(whichComp == co);
            remoteBuses = setdiff(nodeIdx, nodesInThisComp);
            cfc.bus(remoteBuses, BUS_TYPE) = 4;
            if ~any(cfc.bus(nodesInThisComp, BUS_TYPE) == 3)
                [m, imax] = max(cfc.bus(nodesInThisComp, PD));
                cfc.bus(nodesInThisComp(imax), BUS_TYPE) = 3;
%                 continue;
            end
            %     cfc.gen(whichComp == temp(comSize==1), GEN_STATUS) = -1;
            remoteGenBusIdx = ismember(cfc.gen(:, GEN_BUS), remoteBuses);
            cfc.gen(remoteGenBusIdx, GEN_STATUS) = 0;
            
            
            %%  Adjusting load level during cascading failures according the change of system frequency has reference on the book : power system stability and control, Page 392,  P. Kundur
            totalLoadInThisSubG = sum(cfc.bus(nodesInThisComp, PD));
            totalGenInThisSubG = sum(cfc.gen(~remoteGenBusIdx, PG));
            D = 1.5 *  totalLoadInThisSubG / totalGenInThisSubG;
            deltaPl = (totalLoadInThisSubG - totalGenInThisSubG) /totalGenInThisSubG;
            df = - deltaPl / D * 50;
            

                cfc.bus(:, PD) = cfc.bus(:, PD) * sum(cfc.gen(~remoteGenBusIdx, PG)) / sum(cfc.bus(:, PD));
                subG = subgraph(gra, nodesInThisComp);
                subEdges = table2array(subG.Edges);
                for tr = 1 : size(subEdges, 1) * size(subEdges, 2)
                    subEdges(tr) = nodesInThisComp(subEdges(tr));
                end
                branchIdxInThisSubGraph  = ismember(mpc.branch(:,[F_BUS T_BUS]), subEdges,'rows');
            if abs(df) < 2.5    %   tolorance of frequency change.
                cfc.branch(:, BR_STATUS) = cfc.branch(:, BR_STATUS) & branchIdxInThisSubGraph;
                cfc = rundcpf(cfc, mpoption('verbose', 0,'out.all',0));
                branchFlows(branchIdxInThisSubGraph) = cfc.branch(branchIdxInThisSubGraph, 14);
                
            else
                branchFlows(branchIdxInThisSubGraph) = cfc.branch(branchIdxInThisSubGraph, RATE_A) * 1.5;   %   the flows during blackouts are set to this scale to make training the neural net easy.
            end
            
        end
    end
    
    
    mpc.branch(:, BR_STATUS) = statusOfLinesPre & abs(branchFlows) < mpc.branch(:, RATE_A);
    sample = [real(makeSbus(mpc.baseMVA, mpc.bus, mpc.gen) * mpc.baseMVA);statusOfLinesPre;branchFlows];
    count = count + 1;
    cfm(:, count) = sample;
end



%
% while any(statusOfLinesPre - mpc.branch(:, BR_STATUS))
%     statusOfLinesPre = mpc.branch(:, BR_STATUS);
%     mpc = rundcpf(mpc, mpoption('verbose',1));
%     mpc.branch(:, BR_STATUS) = statusOfLinesPre & abs(mpc.branch(:, 14)) < mpc.branch(:, RATE_A);
% end

cfm = cfm(:, 1 : count);
end

