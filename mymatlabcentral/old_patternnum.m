function [xnuc] = patternnum(a,b,c,d)
xnuc=0;
if (a~=c), return; end   % a c are human
if (a~=9 && b~=9 && d~=9)
	       if (a==b)
			if (a==d)
			      xnuc=1;
			else
			      xnuc=2;
			end
	       else
			if (b==d)         % b is chimp, d is zebrafish
				xnuc=3;
			elseif (a==d)
				xnuc=4;
			else
				xnuc=5;
			end		
	       end
elseif(a~=9 && b~=9 && d==9)
				xnuc=6;
end
