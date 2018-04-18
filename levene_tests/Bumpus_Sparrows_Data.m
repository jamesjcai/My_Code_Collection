issurvived=false(49,1);
issurvived(1:22)=true;
isdied=~issurvived;

data=[156	245	31.6	18.5	20.5
154	240	30.4	17.9	19.6
153	240	31	18.4	20.6
153	236	30.9	17.7	20.2
155	243	31.5	18.6	20.3
163	247	32	19	20.9
157	238	30.9	18.4	20.2
155	239	32.8	18.6	21.2
164	248	32.7	19.1	21.1
158	238	31	18.8	22
158	240	31.3	18.6	22
160	244	31.1	18.6	20.5
161	246	32.3	19.3	21.8
157	245	32	19.1	20
157	235	31.5	18.1	19.8
156	237	30.9	18	20.3
158	244	31.4	18.5	21.6
153	238	30.5	18.2	20.9
155	236	30.3	18.5	20.1
163	246	32.5	18.6	21.9
159	236	31.5	18	21.5
155	240	31.4	18	20.7
156	240	31.5	18.2	20.6
160	242	32.6	18.8	21.7
152	232	30.3	17.2	19.8
160	250	31.7	18.8	22.5
155	237	31	18.5	20
157	245	32.2	19.5	21.4
165	245	33.1	19.8	22.7
153	231	30.1	17.3	19.8
162	239	30.3	18	23.1
162	243	31.6	18.8	21.3
159	245	31.8	18.5	21.7
159	247	30.9	18.1	19
155	243	30.9	18.5	21.3
162	252	31.9	19.1	22.2
152	230	30.4	17.3	18.6
159	242	30.8	18.2	20.5
155	238	31.2	17.9	19.3
163	249	33.4	19.5	22.8
163	242	31	18.1	20.7
156	237	31.7	18.2	20.3
159	238	31.5	18.4	20.3
161	245	32.1	19.1	20.8
155	235	30.7	17.7	19.6
162	247	31.9	19.1	20.4
153	237	30.6	18.6	20.4
162	245	32.5	18.5	21.1
164	248	32.3	18.8	20.9];

% SOURCE: https://www.amazon.com/Multivariate-Statistical-Methods-Primer-Third/dp/1584884142

%%

dataz=zscore(data,0,1);
m1=median(dataz(issurvived,:));
m2=median(dataz(isdied,:));
d=[vecnorm(dataz(issurvived,:)-m1,2,2); vecnorm(dataz(isdied,:)-m2,2,2)];

figure;
boxplot(d,isdied)
hold on
plot(isdied+1.25,d,'o')
% Figure 1 of Anderson (2006) Distance-based tests for homogeneity of
% multivariate dispersions
p=anova1(d,issurvived,'display','off')
%%

Dx=interdist(dataz(issurvived,:));
Dy=interdist(dataz(~issurvived,:));

[p,F,Dxv,Dyv]=Gijbels_Omelka_test(Dx,Dy);

% F = 4.319
d=[Dxv;Dyv];
figure;
boxplot(d,isdied)
hold on
plot(isdied+1.25,d,'o')
% Figure 1 of Gijbels and Omelka (2013) Testing for Homogeneity of 
% Multivariate Dispersions Using Dissimilarity Measures
