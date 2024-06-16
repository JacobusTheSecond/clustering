function [segT0, segHacaExp, cnames] = runHaca(tag)

% The value of tag could be set to any integer between 1 and 14
% which correspods to the trial number of subject 86.
% You can derive similar results as shown in http://humansensing.cs.cmu.edu/projects/aca_more_results.html

%% source
wsSrc = mocSegSrc(tag);
[para, paraH] = stFld(wsSrc, 'para', 'paraH');

%% feature
wsData = mocSegData(wsSrc);
[X, segT] = stFld(wsData, 'X', 'segT');
K = conKnl(conDist(X, X), 'nei', .02);
para.nIni = 1;

%% haca
seg0s = segIni(K, paraH(1));
segHaca = segAlg('haca', [], K, paraH, seg0s, segT);

segHacaExp = segNewS(segHaca(end),wsData.sR, "type", "expand");
segT0 = wsData.segT0;
cnames = wsData.cnames;

end