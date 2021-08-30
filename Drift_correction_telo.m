
%%% Change correction_type (default is ini), nuc_ch, max_rot, rot_step_size, files, the first for loop and bckgrnd_multiple as from line 93
%%% written by D.S.JOKHUN on 25/06/2016


clear all
clc
correction_type = 'ini'  %Disable this line if rotational correction is to be made with respect to the previous frame
% correction_type = 'pre'  %Disable this line if rotational correction is to be made with respect to the first frame
nuc_ch=1   %set the channel which will be used to mark the nucleus
max_rot = 15 %max expected rotation of the nucleus in degrees
rot_step_size = 0.005 %the step by which the nucleus will be rotated to find the perfect fit, in degrees


centroid_header = {};
rotation_header = {};
nuc_cen=[];
nuc_rot=[];
cell_no=0;

angles=[];
angles(:,1) = -max_rot:rot_step_size:max_rot;

name = 'cir_ctrl'
files = dir ([name,'*.nd2']);
telo_files = dir ([name,'*_rawtelo.xls']);
for f=1:size(files,1);
    f
    filename = files(f).name

    Reader = bfGetReader (filename);
    OmeMeta = Reader.getMetadataStore();

    MetaData.SeriesCount = Reader.getSeriesCount();
    MetaData.TimePoints = OmeMeta.getPixelsSizeT(0).getValue();
    MetaData.Num_of_Ch = OmeMeta.getPixelsSizeC(0).getValue();
    MetaData.Num_of_Pixels_Z = OmeMeta.getPixelsSizeZ(0).getValue();
    MetaData.Num_of_Pixels_X = OmeMeta.getPixelsSizeX(0).getValue();
    MetaData.Num_of_Pixels_Y = OmeMeta.getPixelsSizeY(0).getValue();
    MetaData.Voxel_Size_X = double(OmeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM)); % in um
    MetaData.Voxel_Size_Y = double(OmeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM)); % in um
%    MetaData.Voxel_Size_Z = double(OmeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROM)); % in um
    MetaData.Plane_Origin_X = double(OmeMeta.getPlanePositionX(0,0).value);
    MetaData.Plane_Origin_Y = double(OmeMeta.getPlanePositionY(0,0).value);
    MetaData.Plane_Origin_Z = double(OmeMeta.getPlanePositionZ(0,0).value);

    MetaData.ChannelID = [];
    for ch_count = 1:MetaData.Num_of_Ch ;
        chID_temp = ['   ' char(num2str(ch_count-1)) '.'] ;
        chNAME_temp= [char(OmeMeta.getChannelName(0,ch_count-1))];
        MetaData.ChannelID = [MetaData.ChannelID  chID_temp chNAME_temp];  % (series 0, channel ch_count)
    end
    MetaData

    
    
    %%% loading 3D frames containing nuclear info
    for iSeries = 1:MetaData.SeriesCount;  %%%choosing a specific XY 3D point from the multipoint image
        iSeries 
        XYZ=[];
        Reader.setSeries(iSeries - 1);
        
            
        raw_outline=[];
        translated_outline=[];
        rotated_outline=[];
        nuc_pos_n = [];
        nuc_pos_n(1:MetaData.TimePoints,1:2)=0;
        
        cell_no=cell_no+1;
        centroid_header{1,(2*cell_no)-1} = filename;
        centroid_header{1,2*cell_no} = filename;
        rotation_header{1,cell_no} = filename;
        
        for iT=1:MetaData.TimePoints;
            iT
            for iCh= nuc_ch; %1:MetaData.Num_of_Ch;
                XYZ_temp =uint16([]);
                for iZ=1:MetaData.Num_of_Pixels_Z;
                    iPlane = Reader.getIndex(iZ-1, iCh-1, iT-1) + 1;     %%% The last '1-1' is for timepoint 0 (the 1st timepoint)
                	XYZ_temp(:,:,iZ)= bfGetPlane(Reader, iPlane);
                end
                XYZ{iCh}=XYZ_temp;   %%% 1st element of XYZ will be a 3D matrix of series 1 (XY1) in Channel 1
                                 %%% 2st element of XYZ will be a 3D matrix of series 1 (XY1) in Channel 2                    
            end 
        
        
       %%% {XYZ} currently has the 3D intensity matrix of the nuclear channel at iT(timepoint iT) from iSeries(multipoint i) found in file f 
       %% PERFORM ANALYSIS BELOW!!!
       
       %% thresholding slice by slice
       nuc_bw=[];
       
       if f==1
           bckgrnd_multiple=5;
       end
       if f==2
           bckgrnd_multiple=6;
       end
       if f==3
           bckgrnd_multiple=3;
       end
       if f==4
           bckgrnd_multiple=3;
       end
       if f==5
           bckgrnd_multiple=1.8;
       end
       if f==6
           bckgrnd_multiple=1.5;
       end
       if f==7
           bckgrnd_multiple=2;
       end
       if f==8
           bckgrnd_multiple=1.2;
       end
       if f==9
           bckgrnd_multiple=1.1;
       end
       if f==10
           bckgrnd_multiple=1.2;
       end

       
       parfor count = 1:MetaData.Num_of_Pixels_Z;
           nuc_medfilt2= medfilt2(XYZ{nuc_ch}(:,:,count),'symmetric');  %performs median filter with a 3x3 neighborhood and pads the border with a symetric reflection
           avg_back_iz= (sum(sum(nuc_medfilt2(4:18,4:18))) + sum(sum(nuc_medfilt2(4:18,end-17:end-3))) + sum(sum(nuc_medfilt2(end-17:end-3,4:18))) + sum(sum(nuc_medfilt2(end-17:end-3,end-17:end-3))))/(15*15*4);   % finding the intensities at the corners (plane by plane). Avoiding the last 3 pixels since median filter is not very good at borders
           thresh_lvl_iz = (bckgrnd_multiple*avg_back_iz)/65535;  %adjust the multiple to adjust the threshold level
           nuc_bw (:,:,count)=im2bw(nuc_medfilt2,thresh_lvl_iz);
       end
       
       nuc_CC = bwconncomp(nuc_bw,26);
       real_nuc_bw =[]; 
       real_nuc_bw(1:MetaData.Num_of_Pixels_Y,1:MetaData.Num_of_Pixels_X,1:MetaData.Num_of_Pixels_Z) = 0;
       for obj_count=1:size(nuc_CC.PixelIdxList,2);
            if size(nuc_CC.PixelIdxList{obj_count},1)>5000
                real_nuc_bw(nuc_CC.PixelIdxList{obj_count})=1;
            end
       end
       
       
       real_nuc_pro = (sum(real_nuc_bw,3));
       real_nuc_pro = imfill(real_nuc_pro,'holes');
       real_nuc_pro = imerode (real_nuc_pro,strel('disk',10)); %erodes ~10 pixels from the border to make sure we are using only strong signals for calculations 

       %getting rid of small isolated dots which might appear after eroding
       nuc_CC = bwconncomp(real_nuc_pro,8); 
       real_nuc_pro(1:size(real_nuc_pro,1),1:size(real_nuc_pro,2)) = 0;
       for obj_count=1:size(nuc_CC.PixelIdxList,2);
            if size(nuc_CC.PixelIdxList{obj_count},1)>500
                real_nuc_pro(nuc_CC.PixelIdxList{obj_count})=1;
            end
       end
       
       %expands the field of view so that further processing doesn't throw anything outside the frame
       real_nuc_pro(51:50+size(real_nuc_pro,1),51:50+size(real_nuc_pro,2))=real_nuc_pro;
       real_nuc_pro(1:50,1:end+50)=0;
       real_nuc_pro(51:end+50,1:50)=0;
       
       real_nuc_pro = imfill(real_nuc_pro,'holes');
       
       stats_cen=regionprops(real_nuc_pro>0,'Centroid');  %finds the centroid of the logical (0 and 1) nuclear projection
       
       %saves the position of nuclear centroid for each timepoint
       nuc_cen (iT,(2*cell_no)-1)= ((stats_cen.Centroid(1,1)-50)* MetaData.Voxel_Size_X)+ MetaData.Plane_Origin_X ; % -50 is because the field of view was expanded by 50 pixels. I want the centroid output file to be according to the original image. +MetaData.Plane_Origin_X is to make the origin same as in Imaris (according to the position on the microscope)
       nuc_cen (iT,2*cell_no)= ((stats_cen.Centroid(1,2)-50)* MetaData.Voxel_Size_Y)+ MetaData.Plane_Origin_Y ; 
       
       nuc_pos_n (iT,1)= ((stats_cen.Centroid(1,1)-50)* MetaData.Voxel_Size_X)+ MetaData.Plane_Origin_X ; 
       nuc_pos_n (iT,2)= ((stats_cen.Centroid(1,2)-50)* MetaData.Voxel_Size_Y)+ MetaData.Plane_Origin_Y ;
       
       % translates the projected nucleus such that its centroid is always at the centre of the image before finding the rotation
       translated_nuc_pro = imtranslate(real_nuc_pro,[-(stats_cen.Centroid(1,1)-(size(real_nuc_pro,2)/2)), -(stats_cen.Centroid(1,2)-(size(real_nuc_pro,1)/2))],'FillValues',0);
       nuc_edge = edge(translated_nuc_pro>0,'sobel');
       nuc_edge = imdilate(nuc_edge,strel('disk',10));  %dilates the edge by ~10 pixels on every side
       nuc_edge = bwmorph(nuc_edge,'thin',inf);
       nuc_edge = imdilate(nuc_edge,strel('disk',10));
       nuc_edge = bwmorph(nuc_edge,'thin',inf);
       nuc_edge = imdilate(nuc_edge,strel('disk',10));
                                                        
       if iT==1
           nuc_t0 = nuc_edge;
       end
       

       % nuc_t0 is the reference image
       % nuc_t1 is the image whose rotation is being measured
       
       nuc_t1 = nuc_edge;

       correlation=[];
       correlation(:,1)=angles(:,1);
       parfor angle_count = 1:size(angles,1)
           nuc_t1_rotated = imrotate (nuc_t1,angles(angle_count,1),'crop');
           correlation(angle_count,2)= corr2(nuc_t1_rotated,nuc_t0);
       end
       

       max_corr=max(correlation(:,2));
       relevant_angles=[];
       for count = 1:size(correlation,1)
           if correlation (count,2)==max_corr
               relevant_angles(size(relevant_angles,1)+1,1)=correlation(count,1);
           end
       end
       
       if correction_type == 'ini'
           nuc_rot(iT,cell_no) = -- median (relevant_angles);  %if you have to rotate nuc_t1 by alpha to obtain nuc_t0, it means nuc_t0 was rotated by -alpha
                                                               % the 2nd '-' is to make the angles compliant with what we see on imaris (since imaris was used to track the telomeres)
       end
       
       if correction_type == 'pre'
           if iT == 1
               nuc_rot(iT,cell_no) = -- median (relevant_angles);
           else
               nuc_rot(iT,cell_no) = -- median (relevant_angles)+ nuc_rot(iT-1,cell_no);
           end
       end
                                                               
                                                               
       if iT==1
           raw_outline (:,:,1)= edge(real_nuc_pro>0,'sobel');
           translated_outline (:,:,1) = edge(translated_nuc_pro>0,'sobel');
           rotated_outline(:,:,1)=edge(imrotate(translated_nuc_pro,median(relevant_angles),'crop')>0,'sobel');           
       else
           raw_outline (:,:,size(raw_outline,3)+1)= edge(real_nuc_pro>0,'sobel');
           translated_outline (:,:,size(translated_outline,3)+1) = edge(translated_nuc_pro>0,'sobel');
           rotated_outline(:,:,size(rotated_outline,3)+1)=edge(imrotate(translated_nuc_pro,median(relevant_angles),'crop')>0,'sobel');
       end
       
       
       if correction_type == 'pre'
           nuc_t0=nuc_t1; %sets the current image as the reference for the next loop
       end
     
       
        
       %% PERFORM ANALYSIS ABOVE!!!            
        end
        
        imtool(flip(sum(rotated_outline,3),1));  %flipped along the horizontal axis to match imaris
        
        rawname = [filename,'raw.tif'];
        translname = [filename,'translated.tif'];
        if correction_type == 'pre'
            rotname = [filename, 'rot_pre.tif'];
        else
            rotname = [filename, 'corrected.tif'];
        end
        
%         imwrite(flip(raw_outline(:,:,1),1),rawname, 'WriteMode','overwrite');
%         imwrite(flip(translated_outline(:,:,1),1),translname, 'WriteMode','overwrite');
        imwrite(flip(rotated_outline(:,:,1),1),rotname, 'WriteMode','overwrite');
        for count=2:size(raw_outline,3)
%         imwrite(flip(raw_outline(:,:,count),1),rawname,'WriteMode','append');
%         imwrite(flip(translated_outline(:,:,count),1),translname,'WriteMode','append');
        imwrite(flip(rotated_outline(:,:,count),1),rotname,'WriteMode','append');
        end
        
        
        
        
        raw_telo_pos = [];
        ini_nuc_pos_n = [];
        
        %% Imports telomere positions
        [~,~,raw_telo_pos_xls]=xlsread(telo_files(f).name);
        num_of_tracks = (size(raw_telo_pos_xls,2)-1)/3
        
        %% saves the x-y coordinates in matrix form
        raw_telo_pos(1:size(raw_telo_pos_xls,1)-2,1)=0:size(raw_telo_pos_xls,1)-3; %writes frame number in 1st column
        
% %% smoothing the nuclear trajectory
% nuc_pos_n(:,1) = smooth(nuc_pos_n(:,1),9,'sgolay',3);
% nuc_pos_n(:,2) = smooth(nuc_pos_n(:,2),9,'sgolay',3);  


        ini_nuc_pos_n(1:size(nuc_pos_n,1),1)=nuc_pos_n(1,1);
        ini_nuc_pos_n(1:size(nuc_pos_n,1),2)=nuc_pos_n(1,2);

        %setting initial nuc point to (0,0)
        %all subsequent points will correspond to nuclear translations from initial point
        nuc_pos_n(:,1)=nuc_pos_n(:,1)-ini_nuc_pos_n(1,1);
        nuc_pos_n(:,2)=nuc_pos_n(:,2)-ini_nuc_pos_n(1,2);

        translated_telo_pos=[];
        translated_telo_pos(1:MetaData.TimePoints,1:num_of_tracks*2)=0;
        corrected_telo_pos=[];
        corrected_telo_pos(1:MetaData.TimePoints,1:num_of_tracks*2)=0;
        xls_format=[];
        xls_format(1,(num_of_tracks*3)+1)=0;
        xls_format(3:2+size(raw_telo_pos,1),1)=0:size(raw_telo_pos,1)-1;
        
        %translational correction of telo
        for count_track=1:num_of_tracks
            for count_frame=1:size(raw_telo_pos_xls,1)-2
%             raw_telo_pos(count_frame,count_track*2)=str2double(raw_telo_pos_xls{count_frame+2,(count_track*3)-1}); %fills the x coordinates 
%             raw_telo_pos(count_frame,(count_track*2)+1)=str2double(raw_telo_pos_xls{count_frame+2,count_track*3}); %fills the y coordinates 
            raw_telo_pos(count_frame,count_track*2)=(raw_telo_pos_xls{count_frame+2,(count_track*3)-1}); %fills the x coordinates
            raw_telo_pos(count_frame,(count_track*2)+1)=(raw_telo_pos_xls{count_frame+2,count_track*3}); %fills the y coordinates 
            end
            translated_telo_pos(:,(count_track*2)-1)=(raw_telo_pos(:,count_track*2)-nuc_pos_n(:,1))-ini_nuc_pos_n(:,1); %drift correcting X according to nuclear translations from initial position and then shifting the whole telo trajectory such that it lies within the nucleus with centre (0,0)
            translated_telo_pos(:,count_track*2)=(raw_telo_pos(:,(count_track*2)+1)-nuc_pos_n(:,2))-ini_nuc_pos_n(:,2); %drift correcting Y
        end

% %% smoothing the nuclear rotation
 nuc_rot(1:MetaData.TimePoints,cell_no)= smooth(nuc_rot(1:MetaData.TimePoints,cell_no),17/MetaData.TimePoints,'rlowess');
        
        %rotational correction of telo
        rotated_x=[];
        rotated_y=[];
        parfor count_frame = 1:size(translated_telo_pos,1)
            theta = (pi/180) * (-nuc_rot(count_frame,cell_no))
            rotation_matrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            frame_to_rot = translated_telo_pos (count_frame,:);
            vec_to_rotate=[];
            %writing all the points as column vectors so that rotation matrix can be applied
            for count_track = 1:num_of_tracks
                vec_to_rotate(1,count_track)= frame_to_rot(1,(count_track*2)-1);
                vec_to_rotate(2,count_track)= frame_to_rot(1,count_track*2);
            end   
            rotated_vec = rotation_matrix * vec_to_rotate;
            rotated_x(count_frame,:)=rotated_vec(1,:);
            rotated_y(count_frame,:)=rotated_vec(2,:);
        end
        
        %writing back the rotated vectors in usual format to save
        for count_frame=1:size(translated_telo_pos,1)
            corrected_telo_pos(count_frame,1:2:num_of_tracks*2)=rotated_x(count_frame,:);
            corrected_telo_pos(count_frame,2:2:num_of_tracks*2)=rotated_y(count_frame,:);
        end
        
        translated_xls=xls_format;
        for count_track=1:num_of_tracks
            xls_format(3:end,(count_track*3)-1)=corrected_telo_pos(:,(count_track*2)-1); 
            xls_format(3:end,count_track*3)=corrected_telo_pos(:,count_track*2);
                        
            translated_xls(3:end,(count_track*3)-1)=translated_telo_pos(:,(count_track*2)-1); 
            translated_xls(3:end,count_track*3)=translated_telo_pos(:,count_track*2);
        end

        
        
%         xlswrite(['Corrected_',filename,'.xls'],xls_format);  %saves the corrected telo pos
%         xlswrite(['Translated_',filename,'.xls'],translated_xls);  %saves the translated telo pos
        save(['Corrected_',filename,'.mat'],'xls_format');
        save(['Translated_',filename,'.mat'],'translated_xls');       
        
        figure('Name',filename)
        count =1:num_of_tracks;
        plot(xls_format(3:end,(count*3)-1),xls_format(3:end,count*3))

        
        incre_rot = nuc_rot(2:end,cell_no) - nuc_rot(1:end-1,cell_no);
        figure ('Name',filename)
        plot(0:size(incre_rot,1)-1,incre_rot)
        
    end

end

figure('Name',name)
for plot_count=1:size(nuc_rot,2)
    plot(1:size(nuc_rot,1),nuc_rot(:,plot_count))
    hold on
end
hold off


nuc_rotation=vertcat(rotation_header,num2cell(nuc_rot));
nuc_centroid=vertcat(centroid_header,num2cell(nuc_cen));
% xlswrite(['nuc_rotation_',name,'.xls'],nuc_rotation);
% xlswrite(['nuc_centroid_',name,'.xls'],nuc_centroid);
save(['nuc_rotation_',name,'.mat'],'nuc_rotation');
save(['nuc_centroid_',name,'.mat'],'nuc_centroid');
