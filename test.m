function test()

numOfBaseSamples = 4;
M = 4;
parfor (bt = 1 : numOfBaseSamples, M)
    DATA_DC(case24_ieee_rts, bt);
end
end

