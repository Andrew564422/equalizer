function varargout = equalizer(varargin)
% EQUALIZER MATLAB code for equalizer.fig
%      EQUALIZER, by itself, creates a new EQUALIZER or raises the existing
%      singleton*.
%
%      H = EQUALIZER returns the handle to a new EQUALIZER or the handle to
%      the existing singleton*.
%
%      EQUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EQUALIZER.M with the given input arguments.
%
%      EQUALIZER('Property','Value',...) creates a new EQUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before equalizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to equalizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help equalizer

% Last Modified by GUIDE v2.5 03-Jan-2022 00:24:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @equalizer_OpeningFcn, ...
                   'gui_OutputFcn',  @equalizer_OutputFcn, ...
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


% --- Executes just before equalizer is made visible.
function equalizer_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to equalizer (see VARARGIN)

% Choose default command line output for equalizer
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes equalizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = equalizer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
                              %% the edit text using to input gain in db%
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g1= get(handles.slider1,'value');
set(handles.text17, 'String',g1);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g2= get(handles.slider2,'value');
set(handles.text18, 'String',g2);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider15_Callback(hObject, eventdata, handles)
% hObject    handle to slider15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g3= get(handles.slider15,'value');
set(handles.text19, 'String',g3);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider14_Callback(hObject, eventdata, handles)
% hObject    handle to slider14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g4= get(handles.slider14,'value');
set(handles.text20, 'String',g4);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider16_Callback(hObject, eventdata, handles)
% hObject    handle to slider16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g5= get(handles.slider16,'value');
set(handles.text21, 'String',g5);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g6= get(handles.slider7,'value');
set(handles.text23, 'String',g6);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g7= get(handles.slider15,'value');
set(handles.text24, 'String',g7);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g8= get(handles.slider9,'value');
set(handles.text25, 'String',g8);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g9= get(handles.slider10,'value');
set(handles.text26, 'String',g9);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

                    %% the edit text using to input the output sample rate%                          
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

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oFs=(2)*str2num(char(get(handles.edit1,'String')));
set(handles.edit1,'string',num2str(oFs));

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oFs=(0.5)* str2num(char(get(handles.edit1,'String')));
set(handles.edit1,'string',num2str(oFs));

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oFs=str2num(char(get(handles.edit1,'String')));
set(handles.edit1,'string',num2str(oFs));

                      %% %choose the sound signal%
% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile({'*.wav'},'File Selector');
handles.fullpathname = strcat(pathname, filename);

set(handles.text3, 'String',handles.fullpathname) %showing fullpathname
guidata(hObject,handles)

% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global player;
global choice ;
global v ;
choice = menu('choose','FIR','IIR');
if (choice == 1)
    fir(hObject , handles)
elseif (choice == 2)
    iir(hObject , handles)
end  
play(player);
pause(10)
play(v);
guidata(hObject,handles)

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global player;
global v;
stop(player);
stop(v);
guidata(hObject,handles)
                             %%  % finit impulse response (FIR)
function fir(hObject, handles)
global player;
global v;
[handles.y,handles.Fs] = audioread(handles.fullpathname);
handles.Y = abs(ifftshift(fft(handles.y)));   %%original signal in freq domain
player = audioplayer(handles.y, handles.Fs);
Fs= handles.Fs
handles.g1=get(handles.slider1,'value');
handles.g2=get(handles.slider2,'value');
handles.g3=get(handles.slider15,'value');
handles.g4=get(handles.slider14,'value');
handles.g5=get(handles.slider16,'value');
handles.g6=get(handles.slider7,'value');
handles.g7=get(handles.slider8,'value');
handles.g8=get(handles.slider9,'value');
handles.g9=get(handles.slider10,'value');
%low pass
cut_off=170; 
orde=16;
a=fir1(orde,cut_off/(handles.Fs/2),'low');
fvtool(a)
y1=handles.g1*filter(a,1,handles.y); %amplify the output signal
x1=filter(a,1,handles.y);
X1= abs(ifftshift(fft(x1)));
% %bandpass1
f1=171;
f2=310;
b1=fir1(orde,[f1/(handles.Fs/2) f2/(handles.Fs/2)],'bandpass');
fvtool(b1)
y2=handles.g2*filter(b1,1,handles.y); %amplify the output signal
x2=filter(b1,1,handles.y);
X2= abs(ifftshift(fft(x2)));
% %bandpass2
f3=311;
f4=600;
b2=fir1(orde,[f3/(handles.Fs/2) f4/(handles.Fs/2)],'bandpass');
fvtool(b2)
y3=handles.g3*filter(b2,1,handles.y);
x3=filter(b2,1,handles.y);
X3= abs(ifftshift(fft(x3)));
% %bandpass3
f4=601;
f5=1000;
b3=fir1(orde,[f4/(handles.Fs/2) f5/(handles.Fs/2)],'bandpass');
fvtool(b3)
y4=handles.g4*filter(b3,1,handles.y);
x4=filter(b3,1,handles.y);
X4= abs(ifftshift(fft(x4)));
% %bandpass4
f5=1001;
f6=3000;
b4=fir1(orde,[f5/(handles.Fs/2) f6/(handles.Fs/2)],'bandpass');
fvtool(b4)
y5=handles.g5*filter(b4,1,handles.y);
x5=filter(b4,1,handles.y);
X5= abs(ifftshift(fft(x5)));
% %bandpass5
f7=3001;
f8=6000;
b5=fir1(orde,[f7/(handles.Fs/2) f8/(handles.Fs/2)],'bandpass');
fvtool(b5)
y6=handles.g6*filter(b5,1,handles.y);
x6=filter(b5,1,handles.y);
X6= abs(ifftshift(fft(x6)));
% %bandpass6
f9=6001;
f10=12000;
b6=fir1(orde,[f9/(handles.Fs/2) f10/(handles.Fs/2)],'bandpass');
fvtool(b6)
y7=handles.g7*filter(b6,1,handles.y);
x7=filter(b6,1,handles.y);
X7= abs(ifftshift(fft(x7)));
% %bandpass7
f11=12001;
f12=14000;
b7=fir1(orde,[f11/(handles.Fs/2) f12/(handles.Fs/2)],'bandpass');
fvtool(b7)
y8=handles.g8*filter(b7,1,handles.y);
x8=filter(b7,1,handles.y);
X8= abs(ifftshift(fft(x8)));
% %bandpass8
f13=14001;
f14=16000;
b8=fir1(orde,[f13/(handles.Fs/2) f14/(handles.Fs/2)],'bandpass');
fvtool(b8)
y9=handles.g9*filter(b8,1,handles.y);
x9=filter(b8,1,handles.y);
X9= abs(ifftshift(fft(x9)));

handles.yT = y1+y2+y3+y4+y5+y6+y7+y8+y9; %Add the amplified-output signals in time domain to form composite signal
handles.YT = abs(ifftshift(fft(handles.yT)));
%plot output signals in time and freq domain
subplot(3,2,1);plot(x1);title('x1');
subplot(3,2,2);plot(X1);title('X1');
subplot(3,2,3);plot(x2);title('x2');
subplot(3,2,4);plot(X2);title('X2');
subplot(3,2,5);plot(x3);title('x3');
subplot(3,2,6);plot(X3);title('X3');
figure
subplot(3,2,1);plot(x4);title('x4');
subplot(3,2,2);plot(X4);title('X4');
subplot(3,2,3);plot(x5);title('x5');
subplot(3,2,4);plot(X5);title('X5');
subplot(3,2,5);plot(x6);title('x6');
subplot(3,2,6);plot(X6);title('X6');
figure
subplot(3,2,1);plot(x7);title('x7');
subplot(3,2,2);plot(X7);title('X7');
subplot(3,2,3);plot(x8);title('x8');
subplot(3,2,4);plot(X8);title('X8');
subplot(3,2,5);plot(x9);title('x9');
subplot(3,2,6);plot(X9);title('X9');

%%Draw and compare the composite signal with the original signal. (in time and frequency)
figure
subplot(2,2,1);
plot(handles.y);title('original in time domain');
subplot(2,2,2);
plot(handles.yT);title('composite in time domain');
subplot(2,2,3);
plot(handles.Y);title('original in freq domain');
subplot(2,2,4);
plot(handles.YT);title('composite in freq domain');
% Play and save the output wave signal. (you can use wavwrite-sound) 
oFs=str2num(char(get(handles.edit1,'String')));
if oFs == 0
    oFs = inputdlg('please enter the output sample rate');
    set(handles.edit1,'string',oFs);
    oFs= str2num(char(get(handles.edit1,'String')));
else
    [P,Q] = rat(oFs/handles.Fs);
    abs(P/Q*handles.Fs-oFs);
end
audiowrite('composit_fir.wav',handles.yT,oFs);
[vx n]=audioread('composit_fir.wav');
VX= abs(ifftshift(fft(vx)));
figure
subplot(2,1,1);
plot(vx);title('after adding output sample rate in time domain');
subplot(2,1,2);
plot(VX);title('after adding output sample rate in time frequency');
v = audioplayer(vx, n);
guidata(hObject,handles)

                             %%   % infinit impulse response (IIR)
function iir(hObject, handles)
global player;
global v;
[handles.y,handles.Fs] = audioread(handles.fullpathname);
handles.Y = abs(ifftshift(fft(handles.y)));   %%original signal in freq domain
player = audioplayer(handles.y, handles.Fs);
Fs= handles.Fs
handles.g1=get(handles.slider1,'value');
handles.g2=get(handles.slider2,'value');
handles.g3=get(handles.slider15,'value');
handles.g4=get(handles.slider14,'value');
handles.g5=get(handles.slider16,'value');
handles.g6=get(handles.slider7,'value');
handles.g7=get(handles.slider8,'value');
handles.g8=get(handles.slider9,'value');
handles.g9=get(handles.slider10,'value');
%low pass
cut_off=170; 
orde=2;
fn=handles.Fs/2;
wn1=170/fn;
[b1,a1]=butter(orde,wn1,'low');    
y1=handles.g1*filter(b1,a1,handles.y);
x1=filter(b1,a1,handles.y); 
X1= abs(ifftshift(fft(x1)));
fvtool(b1,a1)

% %bandpass1
f1=171;
f2=310;
wn2=[f1/fn f2/fn];
[b2,a2]=butter(orde,wn2,'bandpass');
y2=handles.g2*filter(b2,a2,handles.y);
x2=filter(b2,a2,handles.y);
X2= abs(ifftshift(fft(x2)));
fvtool(b2,a2)

% %bandpass2
f3=311;
f4=600;
wn3=[f3/fn f4/fn];
[b3,a3]=butter(orde,wn3,'bandpass');
y3=handles.g3*filter(b3,a3,handles.y);
x3=filter(b3,a3,handles.y);
X3= abs(ifftshift(fft(x3)));
fvtool(b3,a3)

% %bandpass3
f5=601;
f6=1000;
wn4=[f5/fn f6/fn];
[b4,a4]=butter(orde,wn4,'bandpass');
y4=handles.g4*filter(b4,a4,handles.y);
x4=filter(b4,a4,handles.y);
X4= abs(ifftshift(fft(x4)));
fvtool(b4,a4)

% %bandpass4
f7=1001;
f8=3000;
wn5=[f7/fn f8/fn];
[b5,a5]=butter(orde,wn5,'bandpass');
y5=handles.g5*filter(b5,a5,handles.y);
x5=filter(b5,a5,handles.y);
X5= abs(ifftshift(fft(x5)));
fvtool(b5,a5)

% %bandpass5
f9=3001;
f10=6000;
wn6=[f9/fn f10/fn];
[b6,a6]=butter(orde,wn6,'bandpass');
y6=handles.g6*filter(b6,a6,handles.y);
x6=filter(b6,a6,handles.y);
X6= abs(ifftshift(fft(x6)));
fvtool(b6,a6)

% %bandpass6
f11=6001;
f12=12000;
wn7=[f11/fn f12/fn];
[b7,a7]=butter(orde,wn7,'bandpass');
y7=handles.g7*filter(b7,a7,handles.y);
x7=filter(b7,a7,handles.y);
X7= abs(ifftshift(fft(x7)));
fvtool(b7,a7)

% %bandpass7
f13=12001;
f14=14000;
wn8=[f13/fn f14/fn];
[b8,a8]=butter(orde,wn8,'bandpass');
y8=handles.g8*filter(b8,a8,handles.y);
x8=filter(b8,a8,handles.y);
X8= abs(ifftshift(fft(x8)));
fvtool(b8,a8)

% %bandpass8
f15=14001;
f16=16000;
wn9=[f15/fn f16/fn];
[b9,a9]=butter(orde,wn9,'bandpass');
y9=handles.g9*filter(b9,a9,handles.y);
x9=filter(b9,a9,handles.y);
X9= abs(ifftshift(fft(x9)));
fvtool(b9,a9)

handles.yT = y1+y2+y3+y4+y5+y6+y7+y8+y9;  %Add the amplified-output signals in time domain to form composite signal
handles.YT = abs(ifftshift(fft(handles.yT)));
%plot output signal in time and freq domain
subplot(3,2,1);plot(x1);title('x1');
subplot(3,2,2);plot(X1);title('X1');
subplot(3,2,3);plot(x2);title('x2');
subplot(3,2,4);plot(X2);title('X2');
subplot(3,2,5);plot(x3);title('x3');
subplot(3,2,6);plot(X3);title('X3');
figure
subplot(3,2,1);plot(x4);title('x4');
subplot(3,2,2);plot(X4);title('X4');
subplot(3,2,3);plot(x5);title('x5');
subplot(3,2,4);plot(X5);title('X5');
subplot(3,2,5);plot(x6);title('x6');
subplot(3,2,6);plot(X6);title('X6');
figure
subplot(3,2,1);plot(x7);title('x7');
subplot(3,2,2);plot(X7);title('X7');
subplot(3,2,3);plot(x8);title('x8');
subplot(3,2,4);plot(X8);title('X8');
subplot(3,2,5);plot(x9);title('x9');
subplot(3,2,6);plot(X9);title('X9');

%Draw and compare the composite signal with the original signal. (in time and frequency) 
figure
subplot(2,2,1);
plot(handles.y);title('original in time domain');
subplot(2,2,2);
plot(handles.yT);title('composite in time domain');
subplot(2,2,3);
plot(handles.Y);title('original in freq domain');
subplot(2,2,4);
plot(handles.YT);title('composite in freq domain');
% Play and save the output wave signal. (you can use wavwrite-sound) 
oFs=str2num(char(get(handles.edit1,'String')));
if oFs == 0
    oFs = inputdlg('please enter the output sample rate');
    set(handles.edit1,'string',oFs);
    oFs= str2num(char(get(handles.edit1,'String')));
else
    [P,Q] = rat(oFs/handles.Fs);
    abs(P/Q*handles.Fs-oFs);
end
audiowrite('composit_iir.wav',handles.yT,oFs);
[vx n]=audioread('composit_iir.wav');
VX= abs(ifftshift(fft(vx)));
figure
subplot(2,1,1);
plot(vx);title('after adding output sample rate in time domain');
subplot(2,1,2);
plot(VX);title('after adding output sample rate in time frequency');
v = audioplayer(vx, n);
guidata(hObject,handles)
