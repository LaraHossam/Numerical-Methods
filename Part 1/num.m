function varargout = num(varargin)
% NUM MATLAB code for num.fig
%      NUM, by itself, creates a new NUM or raises the existing
%      singleton*.
%
%      H = NUM returns the handle to a new NUM or the handle to
%      the existing singleton*.
%
%      NUM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUM.M with the given input arguments.
%
%      NUM('Property','Value',...) creates a new NUM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before num_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to num_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help num

% Last Modified by GUIDE v2.5 08-Dec-2022 21:11:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @num_OpeningFcn, ...
                   'gui_OutputFcn',  @num_OutputFcn, ...
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


% --- Executes just before num is made visible.
function num_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to num (see VARARGIN)

% Choose default command line output for num
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes num wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = num_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in closed.
function closed_Callback(hObject, eventdata, handles)
% hObject    handle to closed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns closed contents as cell array
%        contents{get(hObject,'Value')} returns selected item from closed


% --- Executes during object creation, after setting all properties.
function closed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to closed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in open.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns open contents as cell array
%        contents{get(hObject,'Value')} returns selected item from open


% --- Executes during object creation, after setting all properties.
function open_CreateFcn(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iterations_Callback(hObject, eventdata, handles)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterations as text
%        str2double(get(hObject,'String')) returns contents of iterations as a double


% --- Executes during object creation, after setting all properties.
function iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function function_box_Callback(hObject, eventdata, handles)
% hObject    handle to function_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of function_box as text
%        str2double(get(hObject,'String')) returns contents of function_box as a double
fun = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function function_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to function_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tolerance as text
%        str2double(get(hObject,'String')) returns contents of tolerance as a double


% --- Executes during object creation, after setting all properties.
function tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fun=get(handles.function_box,'String');
funct_fixed=get(handles.function_box,'String');
S = vectorize(char(handles.fun));
f = str2func(['@(x) ' S]);
es = str2double(get(handles.tolerance,'String'));
max = str2double(get(handles.iterations,'String'));
xu = str2double(get(handles.xupper,'String'));
xl = str2double(get(handles.xlower,'String'));
xold = str2double(get(handles.xold,'String'));
xold2 = str2double(get(handles.xold2,'String'));
if (get(handles.closed,'Value')==2 && get(handles.open,'Value')==1 )
bisection(handles,max,f,es,xu,xl);
end
if (get(handles.closed,'Value')==3 && get(handles.open,'Value')==1)
regula_falsi(handles,max,f,es,xu,xl);
end
if (get(handles.open,'Value')==2 && get(handles.closed,'Value')==1)
fixed_point(handles,max,funct_fixed,f,es,xold);
end
if (get(handles.open,'Value')==3 && get(handles.closed,'Value')==1)
secant(handles,max,f,es,xold,xold2);
end
if (get(handles.open,'Value')==4 && get(handles.closed,'Value')==1)
newton(handles,max,f,es,xold);
end




%set(handles.output_box,'String',handles.fun);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
[file,path] = uigetfile('*.txt');
fullpath = fullfile(path,file);
fid = fopen(fullpath); 
tline = fgetl(fid);
while ischar(tline)
    set(handles.function_box,'String',tline)
    disp(tline)
    tline = fgetl(fid);
end
fclose(fid);
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function xold2_Callback(hObject, eventdata, handles)
% hObject    handle to xold2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xold2 as text
%        str2double(get(hObject,'String')) returns contents of xold2 as a double


% --- Executes during object creation, after setting all properties.
function xold2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xold2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xold_Callback(hObject, eventdata, handles)
% hObject    handle to xold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xold as text
%        str2double(get(hObject,'String')) returns contents of xold as a double


% --- Executes during object creation, after setting all properties.
function xold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xupper_Callback(hObject, eventdata, handles)
% hObject    handle to xupper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xupper as text
%        str2double(get(hObject,'String')) returns contents of xupper as a double


% --- Executes during object creation, after setting all properties.
function xupper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xupper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xlower_Callback(hObject, eventdata, handles)
% hObject    handle to xlower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlower as text
%        str2double(get(hObject,'String')) returns contents of xlower as a double


% --- Executes during object creation, after setting all properties.
function xlower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bisection(handles,max,f,es,xu,xl)
now1 = tic();
% Case bisection:
fxl = f(xl);
fxu = f(xu);



if fxl * fxu < 0
    fprintf('Initial guesses are right. \n');
    xr = (xu + xl) / 2;
ea = 100;
i =0;
str =       sprintf("+----+--------+--------+--------+--------+--------+\n");
str = str + sprintf("| i |   Xl   |   Xu   |  Xr    |  f(Xr)    |   Error       |\n");
str = str + sprintf("+----+--------+--------+--------+--------+--------+\n");
str = str + sprintf("| %d   |%.5f    |%.5f    |%.5f    |%.2f    |%.3f    |\n",i,xl,xu,xr,f(xr),ea);
xrprev=0;
i=1;
array=[xr];
while (ea > es && i ~= max )
    fxr = f(xr);
    fxl = f(xl);
    if fxl * fxr < 0
        xu = xr;
    elseif fxl * fxr > 0
        xl = xr;
    end
    xrprev = xr;
    xr = (xu + xl) / 2;
    
    ea = abs((xrprev - xr) / xr);
    if i>9
      str = str + sprintf("|%d  |%.5f    |%.5f    |%.5f    |%.2f    |%.5f    |\n",i,xl,xu,xr,fxr,ea);  
    else
    str = str + sprintf("| %d   |%.5f    |%.5f    |%.5f    |%.2f    |%.5f    |\n",i,xl,xu,xr,fxr,ea);
    end
    i = i+1;
    array(end+1)=xr;
end
time = toc(now1);
tv = fzero(f,0);
disp(tv);
prec = abs(tv - xr);
str = str + sprintf("+---+--------+--------+--------+--------+--------+\n");
set(handles.output_box,'String',str);
set(handles.root,'String',string(xr));
set(handles.precision,'String',string(prec));
set(handles.execution_time,'String',string(time));

x_axis = [-2:0.1:2];
y_axis= [];
i = 1;
while(i<=41)
    ans = f(x_axis(i));
    y_axis = [y_axis ans];
    i = i+1;
end
size= numel(array);
newy=zeros(1,size);


plot(handles.axes2,x_axis,y_axis,'red','linewidth',2);
 hold on;
stem(array,newy)
hold off;


grid(handles.axes2,'on');
grid(handles.axes2,'minor');
else
    fprintf('Initial guesses are wrong. \n');
    l = msgbox("Initial guesses are wrong. ","Bisection Error","error");
    
end




function regula_falsi(handles,max,f,es,xu,xl)
now1 = tic();
% Case regula falsi:
fxl = f(xl);
fxu = f(xu);

if fxl * fxu < 0
    fprintf('Initial guesses are right. \n');
else
    fprintf('Initial guesses are wrong. \n');
    
end

 xr = xu - (fxu * (xl - xu) / (fxl - fxu));
 array=[xr];
ea = 100;
i =0;
str =       sprintf("+----+--------+--------+--------+--------+--------+\n");
str = str + sprintf("| i |   Xl   |   Xu   |  Xr    |  f(Xr)    |   Error       |\n");
str = str + sprintf("+----+--------+--------+--------+--------+--------+\n");
str = str + sprintf("%d    |%.2f    |%.2f    |%.2f    |%.2f    |%.3f    |\n",i,xl,xu,xr,f(xr),ea);
xrprev=0;
while (ea > es && i ~= max )
    fxr = f(xr);
    fxl = f(xl);
    if fxl * fxr < 0
        xu = xr;
    elseif fxl * fxr > 0
        xl = xr;
    end
    xrprev = xr;
    xr = xu - (fxu * (xl - xu) / (fxl - fxu));
    ea = abs((xrprev - xr) / xr);
    str = str + sprintf("%d    |%.5f    |%.5f    |%.5f    |%.2f    |%.5f    |\n",i,xl,xu,xr,fxr,ea);
    i = i+1;
    array(end+1)=xr
end
time = toc(now1);
tv = fzero(f,0);
disp(tv);
prec = abs(tv - xr);
str = str + sprintf("+---+--------+--------+--------+--------+--------+\n");
set(handles.output_box,'String',str);
set(handles.root,'String',string(xr));
set(handles.precision,'String',string(prec));
set(handles.execution_time,'String',string(time));
x_axis = [-2:0.1:2];
y_axis= [];
i = 1;
while(i<=41)
    ans = f(x_axis(i));
    y_axis = [y_axis ans];
    i = i+1;
end
size= numel(array);
newy=zeros(1,size);

plot(handles.axes2,x_axis,y_axis,'red','linewidth',2);
hold on;
stem(array,newy)
hold off;
grid(handles.axes2,'on');
grid(handles.axes2,'minor');



function fixed_point(handles,max,funct_fixed,f,es,xold)
now1=tic();
syms y;
syms x;
flag=0;
i=0;
ea=100;
 %gx= isolate(funct_fixed, x);


for k=1:length(funct_fixed)
    tf=strcmp( funct_fixed(k),"x");
     if tf==1
        

            newChr = replaceBetween(funct_fixed,k,k,'y');
            m=str2sym(newChr);
            temp=string(isolate(m==0, y))
            s=str2sym(erase(temp,"y == "));
            d = diff(s,x);
            if d==0
                break;
            end
            derv = matlabFunction(d);    
            if abs(derv(xold))<1 
                 flag=1;
                break;
            end
     end
   
                 
end
if flag ==1
     sym x
    gx=matlabFunction(s);
    xnew=gx(xold);
    xprev=0;
    str =       sprintf("+----+--------+--------+--------+--------+\n");
    str = str + sprintf("| i    |   Xi  |  Xi+1    |    Error     |\n");
    str = str + sprintf("+----+--------+--------+--------+--------+\n");
    str = str + sprintf("|%d   |%.5f     |%.5f    |%.3f   |\n",i,xold,xnew,ea);
    i=1;
    while (ea > es && i ~= max )
        xold=xnew;
        xprev=xnew;
        xnew=gx(xold);
        ea = abs((xprev - xnew)/xnew);
        if i> 9
            str = str + sprintf("|%d  |%.5f    |%.5f    |%.5f    |\n",i,xold,xnew,ea);
        else
            str = str + sprintf("|%d   |%.5f    |%.5f    |%.5f    |\n",i,xold,xnew,ea);
        end
      
        i=i+1;
    end
    time = toc(now1);
    tv=fzero(f,0);
    prec = abs(tv - xnew);
   % prec = ea;
    str = str + sprintf("+---+--------+--------+--------+--------+--------+\n");
    set(handles.output_box,'String',str);
    set(handles.root,'String',string(xnew));
    set(handles.precision,'String',string(prec));
    set(handles.execution_time,'String',string(time));
    x_axis = [-10:0.1:10];
    y_axis= [];
    i = 1;
    while(i<=201)
        ans = f(x_axis(i));
        y_axis = [y_axis ans];
        i = i+1;
    end

    plot(handles.axes2,x_axis,y_axis,'red','linewidth',2);
    grid(handles.axes2,'on');
    grid(handles.axes2,'minor');
    
else
  l = msgbox("function does not converge ","Divergion Error","error");
end





        
    
    



function secant(handles,max,f,es,xold,xold2)
now1=tic();
% Case secant:
xnew=xold-((f(xold)*(xold2-xold))/(f(xold2)-f(xold)));
i =0;
ea=100;
str =       sprintf("+----+--------+--------+--------+--------+--------+\n");
str = str + sprintf("| i |   Xi-1   |   Xi  |  Xi+1    |    Error     |\n");
str = str + sprintf("+----+--------+--------+--------+--------+--------+\n");
str = str + sprintf("|%d    |%.5f    |%.5f    |%.5f    |%.3f    |\n",i,xold2,xold,xnew,ea);
xrprev=0;
i=1;
while (ea > es && i ~= max )
   
    xrprev = xnew;
    xold2=xold
    xold=xnew
    xnew=xold-((f(xold)*(xold2-xold))/(f(xold2)-f(xold)));
    ea = abs((xrprev - xnew) / xnew);
    str = str + sprintf("|%d    |%.5f    |%.5f    |%.5f    |%.5f    |\n",i,xold2,xold,xnew,ea);
    i = i+1;
end
time = toc(now1);
tv = fzero(f,0);
disp(tv);
prec = abs(tv - xnew);
str = str + sprintf("+---+--------+--------+--------+--------+--------+\n");
set(handles.output_box,'String',str);
set(handles.root,'String',string(xnew));
set(handles.precision,'String',string(prec));
set(handles.execution_time,'String',string(time));
x_axis = [-10:0.1:10];
y_axis= [];
i = 1;
while(i<=201)
    ans = f(x_axis(i));
    y_axis = [y_axis ans];
    i = i+1;
end

plot(handles.axes2,x_axis,y_axis,'red','linewidth',2);
grid(handles.axes2,'on');
grid(handles.axes2,'minor');


function newton(handles,max,f,es,xold)
now1 = tic();
% Case newton:
syms x;
d = diff(f,x);
derv = matlabFunction(d);
%x(1) = xold;
xnew = xold - (f(xold)/derv(xold));
i =0;
ea = 100;
str =   sprintf("+------+-----------+-----------+-----------+\n");
str = str + sprintf("| i   |     X(i)    |    X(i+1)   |    Error   |\n");
str = str + sprintf("+------+-----------+-----------+-----------+\n");
str = str + sprintf("|%d    |%.5f    |---------------|%.3f    |\n",i,xold,ea);
xprev=0;
i=1;
while(ea > es && i ~= max)
    xold=xnew;
    xprev=xnew;
    xnew =  xold - (f(xold)/derv(xold));
    ea = abs((xprev - xnew)/xnew);
    str = str + sprintf("  |%d    |%.5f    |%.5f    |%.5f    |%.5f  |\n ",i,xold,xnew,ea);
    str = str + sprintf("\n");
    i=i+1;
end
time = toc(now1);
tv = fzero(f,0);
disp(tv);
prec = abs(tv - xnew);
str = str + sprintf("+------+-----------+-----------+-----------+\n");
set(handles.output_box,'String',str);
set(handles.root,'String',string(xnew));
set(handles.precision,'String',string(prec));
set(handles.execution_time,'String',string(time));
x_axis = [-10:0.1:10];
y_axis= [];
j = 1;
while(j<=201)
    ans = f(x_axis(j));
    y_axis = [y_axis ans];
    j = j+1;
end

plot(handles.axes2,x_axis,y_axis,'red','linewidth',2);
grid(handles.axes2,'on');
grid(handles.axes2,'minor');
