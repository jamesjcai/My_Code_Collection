for k=1:5000; ra=ceil(31*rand(31,1));
   xr=x3(ra);yr=y3(ra);
   cc=corrcoef(xr,yr);rhr(k)=cc(1,2);
   sl(k)=rhr(k)*std(yr)/std(xr);
   in(k)=mean(yr)-sl(k)*mean(xr);
   end
% a "for" loop is done with the variable k running through every single integer
% from 1 to 5000 inclusive; the results are stored in 5000-strong columns
% rhr (correlation coefficients), sl (slopes) and in (intercepts)

Now we can check the mean and standard deviation of the obtained distributions. Commands

>>mean(sl)
>>std(sl)

produce values -1.33 and 0.24, which fits the previous result well but supplements it with a value of the standard deviation.
Similarly mean and std values for the intercept and rho are produced. They are, respectively, 98.2 / 9.0 and -0.704 / 0.095.

To make the picture even more impressive, we present Figure 5 as the result of these commands:

>> subplot(1,2,1)
>> hist(sl,30)
>> subplot(1,2,2)
>> hist(in,30)


>>slh=hist(sl,30)
>>slf=find(slh>=70)
>>sum(slh(slf))
show that 4736 out of 5000 trials fall into histogram bins from 7 to 23. To determine the borders of this area, one finds

>>slbinsize=(max(sl)-min(sl))/30;
>>slleftbound= min(sl)+6*slbinsize
>>slrightbound=max(sl)-7*slbinsize

which produces -1.80 and -0.86 as the right and left boundaries for the slope that hold for 4376/5000=97.4% of the trials.

Similar computations, with 

>> inh=hist(in,30);
>> inff=find(inh>60);
>>sum(inh(inff))
will find the left and right and right boundaries for the intercept that hold for 95.1% of the trials (leaving out 8 bins on the left and 5 bins on the right): 81.7 to 117.4. 

This all can be put into the same Figure 6: 

by, first, defining the three regression lines with

>> y3reg=slope*x3+intercept;
>> y3regleft=slleftbound*x3+inleftbound;
>> y3regright=slrightbound*x3+inrightbound;

and then plotting the four sets onto the same figure:

>> plot(x3,y3,'*k',x3,y3reg,'k',x3,y3regleft,'r',x3,y3regright,'r')
%  x3,y3,'*k' presents student data as black stars
%  x3,y3reg,'k' presents the real regression line in black 
%  x3,y3regleft,'g' and x3,y3regright,'g' present the boundary regressions
%  with green lines

The red lines on Figure 6 show the limits of the regression line for 95% of trials.





