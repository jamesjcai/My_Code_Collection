function varargout = testStatusProgressBar(varargin)
% TESTSTATUSPROGRESSBAR MATLAB code for testStatusProgressBar.fig
%      TESTSTATUSPROGRESSBAR, by itself, creates a new TESTSTATUSPROGRESSBAR or raises the existing
%      singleton*.
%
%      H = TESTSTATUSPROGRESSBAR returns the handle to a new TESTSTATUSPROGRESSBAR or the handle to
%      the existing singleton*.
%
%      TESTSTATUSPROGRESSBAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTSTATUSPROGRESSBAR.M with the given input arguments.
%
%      TESTSTATUSPROGRESSBAR('Property','Value',...) creates a new TESTSTATUSPROGRESSBAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testStatusProgressBar_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testStatusProgressBar_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testStatusProgressBar

% Last Modified by GUIDE v2.5 28-Dec-2012 10:07:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testStatusProgressBar_OpeningFcn, ...
                   'gui_OutputFcn',  @testStatusProgressBar_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before testStatusProgressBar is made visible.
function testStatusProgressBar_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testStatusProgressBar (see VARARGIN)

% Choose default command line output for testStatusProgressBar
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes testStatusProgressBar wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = testStatusProgressBar_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonStart.
function pushbuttonStart_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myProgress = statusProgressBar(10, ...
    handles.backgroundProgressBar, ...
    handles.movingProgressBar);
handles.myProgress = myProgress;
guidata(hObject, handles);


% --- Executes on button press in pushbuttonNextStep.
function pushbuttonNextStep_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNextStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.myProgress = handles.myProgress.nextStep();
guidata(hObject, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.myProgress);
