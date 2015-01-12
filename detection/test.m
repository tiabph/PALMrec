filepath = 'E:\a647.tif';
% prompt={'Enter the number of image frames you want to process:'};
imginfoList = imfinfo(filepath);
handles.File = filepath;
handles.TotalImages = length(imginfoList);
defans={num2str(handles.TotalImages)};
% info = inputdlg(prompt, 'Input for process...!', 1, defans);
info = defans;
if ~isempty(info) %&& str2double(info{2,1})==0
    number=str2double(info{1,1});
    [handles.A,handles.ImageNumber]=tiffread(handles.File,[1 number]);
%     pause(eps)
%     if handles.ImageNumber>1;
%        set(handles.slider_image,'visible','on');
%        step=1.0/(handles.ImageNumber-1);
%        set(handles.slider_image,'max',1.0,'min',0,'SliderStep',[step 5*step],'value',0);
%     else
%        set(handles.slider_image,'visible','off');
%     end  
%     set(handles.text10,'string','Particle detection ...')
    handles.V=[];
%     axes(handles.axes1)
%     if handles.parameter.display==1
%         mywaitbar(0,handles.axes2,'');
%     end
    [row column]=size(handles.A(:,:,1));
    handles.data_w=[];
    % handles.data_wavelet=uint16(zeros(row,column,handles.ImageNumber));
    handles.data_wavelet=uint16(zeros(row,column,1));
    %%
    for i=1:handles.ImageNumber  
        data_w=Detection(handles.A(:,:,i),handles.parameter);
        handles.data_wavelet(:,:,1)=uint16(data_w);
        handles.V{i}=FindParticles_para(handles.parameter,data_w,0,0,1,5,5,handles.A(:,:,1));
        IX=handles.V{i}(:,3)>handles.parameter.detection.trash_dim;
        handles.V{i}=handles.V{i}(IX,:);       
%         if handles.parameter.display==1
%         cla(handles.axes1);
%         imshow(handles.A(:,:,i),[handles.contrastLow handles.contrastHigh]);
%         hold on
%         plot(handles.V{i}(:,1),handles.V{i}(:,2),'r*','MarkerSize',2);
%         hold off
%         mywaitbar(i/handles.ImageNumber,handles.axes2,[num2str(floor(i/handles.ImageNumber*100)),'%']);
%         set(handles.slider_image,'value',(i-1)/(handles.ImageNumber-1));
%         end 
%         pause(eps)
%         fstr=strcat('Image #:',num2str(i),'/',num2str(handles.ImageNumber));
%         set(findobj('Tag','text6'),'String',fstr);
        if(mod(i,100)==0)
            i
        end
    end
%     set(handles.Load_std_detection,'enable','on');
%     set(handles.View_detection,'enable','on');
%     set(handles.Edit_particle,'enable','on');
%     set(handles.Particle_number,'enable','on');
%     set(handles.Particle_intensity,'enable','on');
%     set(handles.text10,'string','Detection finished!');   
%     handles.parameter.detection_flag=1;
%     guidata(hObject,handles)
end