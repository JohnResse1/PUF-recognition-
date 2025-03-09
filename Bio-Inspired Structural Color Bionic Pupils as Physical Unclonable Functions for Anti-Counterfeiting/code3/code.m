function varargout = code(varargin)
% CODE M-file for code.fig
%      CODE, by itself, creates a new CODE or raises the existing
%      singleton*.
%
%      H = CODE returns the handle to a new CODE or the handle to
%      the existing singleton*.
%
%      CODE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CODE.M with the given input arguments.
%
%      CODE('Property','Value',...) creates a new CODE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before code_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to code_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help code

% Last Modified by GUIDE v2.5 07-May-2020 17:46:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @code_OpeningFcn, ...
                   'gui_OutputFcn',  @code_OutputFcn, ...
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


% --- Executes just before code is made visible.
function code_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to code (see VARARGIN)

% Choose default command line output for code
handles.output = hObject;
clc;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes code wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = code_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global str
global filenamestr
global I2;
[filename,pathname]=uigetfile({'*.bmp';'*.jpg';'*.gif'},'选择图片');
if isequal(filename,0)
    disp('Users Selected Canceled');
else
str=[pathname,filename];
filenamestr=filename;
im = imread(str);
I2 = imread(str);
axes(handles.axes1);%axes1是坐标轴的标示
imshow(im);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles) %识别
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global str
global template
global mask
global filenamestr
testimage=str;

hmthresh = 0.3;
write = 0;
nname=filenamestr(1:4);
samep=1; %判断是不是要与当前同一个人对比
if samep
    InputPath=['.\',nname,'\']; %同一个人
else
    InputPath='.\0024\'; %不同人
end
if exist(InputPath)
% [result,time] = final1(str)
templatetest=template;
masktest=mask;
tic
shibie();
axes(handles.axes12);
pic=[InputPath,result];
imshow(pic);title('匹配到虹膜');
else
    result='o~o, No match found!';
end
set(handles.text2,'String',result);
t=toc;
disp(['识别用时：',num2str(t)])

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles) %图片运算
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global str
global I2;
%I2=imread('image004.jpg');
% axes(handles.axes2);
% imshow(I2);

eI=edge(I2,'canny', 0.2);
axes(handles.axes3);
imshow(eI);title('canny边缘提取');
% 利用hough变换找到图像中的一个圆
[y0detect,x0detect,Accumulator] = houghcircle(eI,45,4);

axes(handles.axes4);
imshow(I2);
hold on;
for i=1:length(y0detect)   
    plot(x0detect,y0detect,'.r');hold on;
end
% figure;imshow(I2)
axes(handles.axes13);
imshow(Accumulator,[]);
[r,c]=size(I2);

M = circle( c,r,x0detect,y0detect,45);
axes(handles.axes5);
imshow(M,[]);

outI=M.*double(I2);
axes(handles.axes6);
imshow(outI,[]);

outI2=(1-M).*double(I2);
axes(handles.axes7);
imshow(outI2,[]);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles) %定位
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global str
global DIAGPATH % path for writing diagnostic images
%DIAGPATH = 'C:\Documents and Settings\Administrator\桌面\iris';
% DIAGPATH = 'D:\虹膜识别\算法\iris-张冲\template';
DIAGPATH = '.\0023\template';
eyeimage_filename=str;
write=0;
dingwei();

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles) %归一化
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global str
global polar_array
global noise_array
eyeimage_filename=str;
%参数设置
%normalisation parameters
radial_res = 100;
angular_res = 240;
write=0;
% with these settings a 9600 bit iris template is created
guiyihua();

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)%特征提取
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%feature encoding parameters
global str
global polar_array
global noise_array
global template
global mask
eyeimage_filename=str;
nscales=1;
minWaveLength=18;
mult=1; % not applicable if using nscales = 1
sigmaOnf=0.5;
tezhengtiqu()

% --- Executes during object creation, after setting all properties.
function axes10_CreateFcn(hObject, eventdata, handles) %归一化
% hObject    handle to axes10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes10


% --- Executes on button press in pushbutton8.

% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles) %退出
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
close all;


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4
