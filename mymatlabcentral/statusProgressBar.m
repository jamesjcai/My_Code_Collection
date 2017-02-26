classdef statusProgressBar
    
    properties (SetAccess = private, GetAccess = public)
        
        currentProgress = 0;
        totalSteps;
        sizeStep;  %width to hide to each step
        
        backgroundTextBoxHandles;
        foregroundTextBoxHandles;
        percentLabelHandles;
        
        backgroundPosition;  %[x,y,width,height]
        foregroundPosition;  %[x,y,width,height]
        stop = false;  %will be true once we reached the max
        
    end
    
    
    methods
        
        function obj = statusProgressBar(totalSteps, backgroundHandles, foregroundHandles)
            %constructor
            obj.totalSteps = totalSteps;
            obj.backgroundTextBoxHandles = backgroundHandles;
            obj.foregroundTextBoxHandles = foregroundHandles;
            
            %make background textBox visible
            obj.backgroundPosition = get(backgroundHandles,'position');
            set(backgroundHandles,'visible','on');
            
            %make progress bar textBox visible and init size
            obj.foregroundPosition = [obj.backgroundPosition([1 2]), ...
                0.001, obj.backgroundPosition(4)];
            set(foregroundHandles,'position',obj.foregroundPosition);
            set(foregroundHandles,'visible','on');
                         
            %calculate size of step
            totalWidth = obj.backgroundPosition(3);
            obj.sizeStep = totalWidth/totalSteps;

        end
        
        function delete(obj)
            %destructor
            
            set(obj.foregroundTextBoxHandles,'visible','off');
            set(obj.backgroundTextBoxHandles,'visible','off');
            
        end
        
        function obj = nextStep(aIn)
            obj = aIn;
            
            if obj.stop
                return
            end
            
            obj.currentProgress = aIn.currentProgress + 1;
            
            %define new size of progress bar
            obj.foregroundPosition(3) = obj.foregroundPosition(3) + ...
                obj.sizeStep;
            
            %plot new size of progress bar
            obj = refreshProgressBar(obj);
            
            %plot current percentage inside the progress bar
            obj = refreshProgressLabel(obj);
            
            %check progress
            obj = checkProgress(obj);
            
        end
        
        function obj = diaplayProgress(aIn)
            
            obj = aIn;
            set(obj.foregroundTextBoxHandles,'position', obj.foregroundPosition);
            
        end
        
        
        function obj = checkProgress(aIn)
            %Make sure we can't pass the 100% border line
            
            obj = aIn;
            
            if aIn.currentProgress == aIn.totalSteps
                obj.stop = true;
            end
            
        end
        
        function obj = refreshProgressLabel(aIn)
            
           obj = aIn; 
           currentPosition = get(obj.foregroundTextBoxHandles,'position'); 
           currentWidth = currentPosition(3);
           maxPosition = get(obj.backgroundTextBoxHandles,'position'); 
           maxWidth = maxPosition(3);
           
           ratio = (currentWidth / maxWidth ) * 100.0;
           percent = int16(ratio);
           
           str = sprintf('%d %%',percent);
           set(obj.foregroundTextBoxHandles,'string',str);
           
        end
        
        function obj = refreshProgressBar(aIn)
            %will refresh the plot using the new size
            obj = aIn;
            set(obj.foregroundTextBoxHandles,'position', obj.foregroundPosition);
        end
        
        
        
    end
    
    
end