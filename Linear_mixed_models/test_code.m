y=[21 19 20 22 14 15 13 16 14 17 15 17 12 11 12 14 16 20 18 19 14 14 14 12]';
brand=[1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2]';
site=[1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3]';

fitlm([brand site],y,'y~x1*x2')

tbl = table(brand,site,y);

lme = fitlme(tbl,'y~brand+(1|site)')

