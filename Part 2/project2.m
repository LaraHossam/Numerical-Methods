function varargout = project2(varargin)
% PROJECT2 MATLAB code for project2.fig
%      PROJECT2, by itself, creates a new PROJECT2 or raises the existing
%      singleton*.
%
%      H = PROJECT2 returns the handle to a new PROJECT2 or the handle to
%      the existing singleton*.
%
%      PROJECT2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECT2.M with the given input arguments.
%
%      PROJECT2('Property','Value',...) creates a new PROJECT2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before project2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to project2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help project2

% Last Modified by GUIDE v2.5 27-Dec-2022 15:58:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @project2_OpeningFcn, ...
                   'gui_OutputFcn',  @project2_OutputFcn, ...
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


% --- Executes just before project2 is made visible.
function project2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to project2 (see VARARGIN)

% Choose default command line output for project2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes project2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = project2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = get(handles.edit1,'String');
s


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function number_equ_Callback(hObject, eventdata, handles)
% hObject    handle to number_equ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_equ as text
%        str2double(get(hObject,'String')) returns contents of number_equ as a double


% --- Executes during object creation, after setting all properties.
function number_equ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_equ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function system_Callback(hObject, eventdata, handles)
% hObject    handle to system (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of system as text
%        str2double(get(hObject,'String')) returns contents of system as a double


% --- Executes during object creation, after setting all properties.
function system_CreateFcn(hObject, eventdata, handles)
% hObject    handle to system (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.txt');
fullpath = fullfile(path,file);
fileID = fopen(fullpath); 
% Use textscan to read the file, with '%s' as the format specifier to read
% each line as a string
lines = textscan(fileID, '%s', 'Delimiter', '\n');

% Close the file
fclose(fileID);

% Convert the cell array of strings to a regular array
lines = lines{1};


num_iter = lines{1};
set(handles.number_equ,'String',num_iter);

method = lines{2};

if(isequal(method,'Gaussian-elimination '))
    disp("Gaussian elimination method");
    set(handles.menu,'Value',1);
    str ="";
    for i=3:(length(lines))
    str = str + sprintf("%s\n",lines{i}) ;
    end
    set(handles.system,'String',str);
end
if(isequal(method,'Gauss-Seidel '))
    disp("Gauss Seidel method");
    set(handles.menu,'Value',2);
    str ="";
    for i=3:(length(lines)-1)
    str = str + sprintf("%s\n",lines{i}) ;
    end
    set(handles.system,'String',str);
    set(handles.initial,'String',lines{length(lines)});
end
if(isequal(method,'ALL '))
    disp("All of the methods");
    set(handles.menu,'Value',3);
    str ="";
    for i=3:(length(lines)-1)
    str = str + sprintf("%s\n",lines{i}) ;
    end
    set(handles.system,'String',str);
    set(handles.initial,'String',lines{length(lines)});
end





% --- Executes on selection change in menu.
function menu_Callback(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu


% --- Executes during object creation, after setting all properties.
function menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initial_Callback(hObject, eventdata, handles)
% hObject    handle to initial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initial as text
%        str2double(get(hObject,'String')) returns contents of initial as a double


% --- Executes during object creation, after setting all properties.
function initial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iter_Callback(hObject, eventdata, handles)
% hObject    handle to iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iter as text
%        str2double(get(hObject,'String')) returns contents of iter as a double


% --- Executes during object creation, after setting all properties.
function iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tol_Callback(hObject, eventdata, handles)
% hObject    handle to tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tol as text
%        str2double(get(hObject,'String')) returns contents of tol as a double


% --- Executes during object creation, after setting all properties.
function tol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in solve.
function solve_Callback(hObject, eventdata, handles)
% hObject    handle to solve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.num = get(handles.number_equ,'String');
handles.system = get(handles.system,'String');
handles.method = get(handles.menu,'Value'); % 1 for Gauss elimination, 2 for gauss seidel, 3 for all
maxiter = get(handles.iter,'String');
tol = get(handles.tol,'String');
init = get(handles.initial,'String');
handles.system = string(handles.system);
arr = [];
str = "";
%%
% for i=1:str2num(handles.num)
%     if(i==str2num(handles.num))
%             arr = [arr;strip(handles.system(i))];
%     else 
%     arr = [arr;strip(handles.system(i))+","];
%     end
% end
%%
syms a b c d e f g 
disp("Trial here!");
var = [a, b, c, d, e, f, g];
matrix = zeros(str2double(handles.num),str2double(handles.num));
for j=1:str2num(handles.num)
c = 1;
for i=1:str2num(handles.num)
    [other, coeff] = coeffs(str2sym(handles.system(i)),var(j));
    if(~isnan(str2double(char(other(1)))))
        matrix(c,j)=other(1);
        c = c + 1;
    else 
        matrix(c,j)=0;
        c = c + 1;
    end 
end
end
disp("A Matrix")
disp(matrix)
disp("b Matrix")
bMatrix = zeros(str2double(handles.num),1);
c = 1;
for i=1:str2num(handles.num)
    temp = coeffs(str2sym(handles.system(i)),'all');
        bMatrix(c,1)=temp(end,end);
        c = c + 1;
     
end
bMatrix = -1.*bMatrix;
disp(bMatrix)

if(handles.method == 1)
    gauss_elim(handles,matrix,bMatrix);
end
if(handles.method ==  2)
    gauss_seidel(handles,matrix,bMatrix,maxiter,tol,init);
end
if(handles.method ==  3)
    gauss_seidel(handles,matrix,bMatrix,maxiter,tol,init);
    gauss_elim(handles,matrix,bMatrix);

end

function gauss_seidel(handles,A,b,maxiter,tol,init)
now1 = tic();
% Case gauss-seidel:
disp("in gauss-seidel");
str = "";
n = length(str2num(init));
disp(n)
disp(A)
disp(b)

% Set the convergence tolerance and maximum number of iterations
maxiter = str2num(maxiter);
% Perform the Gauss-Seidel iterations
% Set the initial guess for the solution
tokens = split(init, ' ');
x = cellfun(@str2double, tokens);
error  = 100;
iter = 1
disp("This is tolerance")
disp(tol)
res = zeros(n,maxiter);
temp = [];
while(iter < maxiter & error > str2double(tol))
   str = str+ sprintf("%d -->",iter);
   xold = x;
   disp(xold)
   for i = 1:n
      sigma = 0;
      for j = 1:n
         if j ~= i
            sigma = sigma + A(i,j)*x(j);
         end
      end
      x(i) = (b(i) - sigma)/A(i,i);
      temp(i) = abs((x(i)-xold(i))/x(i));
   end
   disp(temp)
   error = max(temp);
   disp(error)
   for k = 1:n
        str = str + sprintf(" %.5f |",x(k));
        res(k,iter) = x(k);
   end
       str = str + sprintf("\n");
       iter = iter + 1;
end

% Curve between the number of iterations and the obtained root value at this 
% iteration for all the methods in the same graph. 
s = "";
for i=1:n
    row = res(i, :);
    s = sprintf("Var %d",i);
    figure('Name',s,'NumberTitle','off')
    plot(row(1:iter));
    hold on;
    legend(s)
end

% Print the solution and number of iterations
time = toc(now1);
set(handles.exec2,'String',string(time));
set(handles.out2,'String',str);
set(handles.num_iter2,'String',string(iter));

% Open a file for writing
fid = fopen('gauss_seidel_output.txt', 'w');

% Write the text to the file
fprintf(fid, '%s', str);

% Close the file
fclose(fid);






function gauss_elim(handles,A,b)
    now1 = tic();
    % Case gauss-elimination:
    disp("in gauss-elimination");
    % Define the coefficient matrix and the right-hand side of the system
    str = "";
% Number of equations
n = size(A, 1);

% Gauss elimination method
for i = 1:n-1
    % Partial pivoting
    [~, pivot] = max(abs(A(i:n, i)));
    pivot = pivot + i - 1;
    if pivot ~= i
        A([i pivot],:) = A([pivot i],:);
        b([i pivot]) = b([pivot i]);
    end
    
    % Elimination
    for j = i+1:n
        factor = A(j,i) / A(i,i);
        A(j,:) = A(j,:) - factor * A(i,:);
        b(j) = b(j) - factor * b(i);
    end
end

% Back substitution
x = zeros(n, 1);
for i = n:-1:1
    x(i) = (b(i) - A(i,i+1:n)*x(i+1:n)) / A(i,i);
    str= str + sprintf(" var (%d) = %f\n",i,x(i));
end

% Display the solution
disp(x)
time = toc(now1);
set(handles.exec,'String',string(time));
set(handles.out,'String',str);
% Open a file for writing
fid = fopen('gauss_elimination_output.txt', 'w');

% Write the text to the file
fprintf(fid, '%s', str);

% Close the file
fclose(fid);


function output_Callback(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output as text
%        str2double(get(hObject,'String')) returns contents of output as a double


% --- Executes during object creation, after setting all properties.
function output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function out_Callback(hObject, eventdata, handles)
% hObject    handle to out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of out as text
%        str2double(get(hObject,'String')) returns contents of out as a double


% --- Executes during object creation, after setting all properties.
function out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function out2_Callback(hObject, eventdata, handles)
% hObject    handle to out2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of out2 as text
%        str2double(get(hObject,'String')) returns contents of out2 as a double


% --- Executes during object creation, after setting all properties.
function out2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to out2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
