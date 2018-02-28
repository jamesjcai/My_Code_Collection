
load politeness_data.mat

boxplot(T.frequency,{T.gender,T.attitude})


mdl=fitlme(T,'frequency~attitude + (1|subject) + (1|scenario)')

mdl=fitlme(T,'frequency~attitude + gender + (1|subject) + (1|scenario)')