function varargout = Decrypt(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Decrypt_OpeningFcn, ...
                   'gui_OutputFcn',  @Decrypt_OutputFcn, ...
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

function Decrypt_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;


guidata(hObject, handles);


function varargout = Decrypt_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



function pushbutton1_Callback(hObject, eventdata, handles)

global img message
[f p]=uigetfile({'*.png';'*jpg';'*.bmp'},'Select the embed image');
BIT=2;
sb=BIT-1;
if isequal(f,0) || isequal(p,0)
    set(handles.text1,'string','Retrive the Data')
else
    I=imread([p f]);    
if size(I,3)==3
    s=I(:,:,1);
    g=I(:,:,2);
    d=I(:,:,3);
else
    s=I;
end
N=s(size(s,1),size(s,2));
k=1;
k1=BIT;
  for i=1:size(s,1)
      for j=1:size(s,2)
          if k<N+1
          B=dec2bin(s(i,j),8);
          if k1>N
             temp=k1-N;
            b(1,k:(k+temp))=B(end-temp:end);
         else
          b(1,k:k1)=B(end-sb:end);
         end
          k=k+BIT;
          k1=k1+BIT;
          end
      end
  end
 i=1;
j=8;
k=1;
while j<=length(b)
B(i,:)=b(k:j);
i=i+1;
k=1+j;
j=j+8;
end
if size(I,3)==3
OI=cat(3,s,g,d);
else
    OI=s;
end
OS = bin2dec(num2str(B))
OO = OS(end);
O = char(OS)';
O1 = OS(end);
Retrieved_Message = ''
O
O1

axes(handles.axes1),imshow(OI,[]); title('Recieved image');
end



function pushbutton2_Callback(hObject, eventdata, handles)

global Stand I1 I2 I3 I4 I11 I31 I21 I41 II Mf1
I11 = imrotate(I11,180);
I21 = imrotate(I21,180);
I31 = imrotate(I31,180);
I41 = imrotate(I41,180);
axes(handles.axes2);
imshow(Mf1);title('Cover image')
I_R = (II-Stand)+[I11 I31 ; I21 I41];

axes(handles.axes3);
imshow(I_R);axis off;
title('Retrieved Image');
