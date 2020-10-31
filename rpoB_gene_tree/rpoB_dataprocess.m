a=["accolens","afermentans subsp. Afermentans","afermentans subsp. Lipophilum","ammoniagenes","amycolatum","argentoratense","aurimucosum","auris","auriscanis","bovis","callunae","camporealensis","capitovis","confusum","coyleae","cystitidis","diphtheriae","durum","efficiens","falsenii","felinum","flavescens","freneyi","glucuronolyticum","glutamicum","imitans","jeikeium","kroppenstedtii","kutscheri","lipophiloflavum","macginleyi","mastitidis","matruchotii","minutissimum","mucifaciens","mycetoides","phocae","pilosum","propinquum","pseudodiphtheriticum","pseudotuberculosis","renale","riegelii","seminale","simulans","singulare","spheniscorum","striatum","sundsvallense","terpenotabidum","testudinoris","thomssenii","ulcerans","urealyticum","variabile","vitaeruminis","xerosis"];
b=["AY492242","AY492265","AY492266","AY492243","AY492241","AY492249","AY492282","AY492234","AY492244","AY492236","AY492245","AY492246","AY492247","AY492248","AY492250","AY492251","AY492230","AY492252","AP005215","AY492253","AY492254","AY492255","AY492237","AY492256","AP022856","AY492259","AY492231","AY492274","AY492257","AY492260","AY492276","AY492281","AY492238","AY492235","AY492261","AY492262","AY492277","AY492258","AY492279","AY492232","AY492239","AY492240","AY492278","AY492263","AY492264","AY492280","AY492283","AY492267","AY492268","AY492269","AY492284","AY492270","AY492271","AY492275","AY492272","AY492273","AY492233"];

for k=20:length(b)
    pause(2)
    k
c=getgenbank(b(k));
namesx=sprintf('C.%s',a(k));
fastawrite('bbb.fas',namesx,c.Sequence);
end
