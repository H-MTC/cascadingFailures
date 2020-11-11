function [ACCFM, net,err] = Main()
ACCFM = acDataForTraining('case24_ieee_rts',1700);
[net,err] = ACCF_BP(ACCFM);
end

