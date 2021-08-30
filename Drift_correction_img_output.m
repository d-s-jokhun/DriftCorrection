%%% XY drifts corrects 2d and 3d files_with rotational correction
%%% Change correction_type (default is ini), nuc_ch, max_rot, rot_step_size, files and the first for loop
%%% written by D.S.JOKHUN on 10/06/2016


clear all
clc
correction_type = 'ini'  %Disable this line if rotational correction is to be made with respect to the previous frame
% correction_type = 'pre'  %Disable this line if rotational correction is to be made with respect to the first frame
nuc_ch=1   %set the channel which will be used to mark the nucleus
max_rot = 45 %max expected rotation of the nucleus in degrees
rot_step_size = 1 %the step by which the nucleus will be rotated to find the perfect fit, in degrees


centroid_header = {};
rotation_header = {};
nuc_cen=[];
nuc_rot=[];
cell_no=0;

angles=[];
angles(:,1) = -max_rot:rot_step_size:max_rot;


files = dir ('*30s.tif');
% telo_files = dir ([name,'_rawtelo.xls']);
for f=1:size(files,1);
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
       parfor count = 1:MetaData.Num_of_Pixels_Z;
           nuc_medfilt2= medfilt2(XYZ{nuc_ch}(:,:,count),'symmetric');  %performs median filter with a 3x3 neighborhood and pads the border with a symetric reflection
           avg_back_iz= (sum(sum(nuc_medfilt2(4:18,4:18))) + sum(sum(nuc_medfilt2(4:18,end-17:end-3))) + sum(sum(nuc_medfilt2(end-17:end-3,4:18))) + sum(sum(nuc_medfilt2(end-17:end-3,end-17:end-3))))/(15*15*4);   % finding the intensities at the corners (plane by plane). Avoiding the last 3 pixels since median filter is not very good at borders
           thresh_lvl_iz = (3*avg_back_iz)/65535;  %adjust the multiple to adjust the threshold level
           nuc_bw (:,:,count)=im2bw(nuc_medfilt2,thresh_lvl_iz);
       end
       
       nuc_CC = bwconncomp(nuc_bw,26);
       real_nuc_bw =[]; 
       real_nuc_bw(1:MetaData.Num_of_Pixels_Y,1:MetaData.Num_of_Pixels_X,1:MetaData.Num_of_Pixels_Z) = 0;
       for obj_count=1:size(nuc_CC.PixelIdxList,2);
            if size(nuc_CC.PixelIdxList{obj_count},1)>50000
                real_nuc_bw(nuc_CC.PixelIdxList{obj_count})=1;
            end
       end
       
       
       real_nuc_pro = (sum(real_nuc_bw,3));
       real_nuc_pro = imerode (real_nuc_pro,ones(7)); %erodes ~3 pixels from the border to make sure we are using only strong signals for calculations 
       %getting rid of small isolated dots which might appear after eroding
       nuc_CC = bwconncomp(real_nuc_pro,8); 
       real_nuc_pro(1:size(real_nuc_pro,1),1:size(real_nuc_pro,2)) = 0;
       for obj_count=1:size(nuc_CC.PixelIdxList,2);
            if size(nuc_CC.PixelIdxList{obj_count},1)>500
                real_nuc_pro(nuc_CC.PixelIdxList{obj_count})=1;
            end
       end
       
       
       stats_cen=regionprops(real_nuc_pro>0,'Centroid');  %finds the centroid of the logical (0 and 1) nuclear projection
       
       %saves the position of nuclear centroid for each timepoint
       nuc_cen (iT,(2*cell_no)-1)= stats_cen.Centroid(1,1)* MetaData.Voxel_Size_X ; 
       nuc_cen (iT,2*cell_no)= stats_cen.Centroid(1,2)* MetaData.Voxel_Size_Y ;
       
       nuc_pos_n (iT,1)= stats_cen.Centroid(1,1)* MetaData.Voxel_Size_X ; 
       nuc_pos_n (iT,2)= stats_cen.Centroid(1,2)* MetaData.Voxel_Size_Y ;
       
       % translates the raw channel such that the centroid of the nucleus is always at the centre of the image
       % translates the raw projected nucleus as well
       translated = imtranslate(XYZ{nuc_ch},[-(stats_cen.Centroid(1,1)-(MetaData.Num_of_Pixels_X/2)), -(stats_cen.Centroid(1,2)-(MetaData.Num_of_Pixels_Y/2))],'FillValues',0);
       
       translated_nuc_pro = imtranslate(real_nuc_pro,[-(stats_cen.Centroid(1,1)-(MetaData.Num_of_Pixels_X/2)), -(stats_cen.Centroid(1,2)-(MetaData.Num_of_Pixels_Y/2))],'FillValues',0);  %will be used to find rotation
       nuc_edge = edge(translated_nuc_pro>0,'sobel');
       nuc_edge = imdilate(nuc_edge,ones(11));
                                                        
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
       nuc_rot(iT,cell_no) = - median (relevant_angles);  %if you have to rotate nuc_t1 by alpha to obtain nuc_t0, it means nuc_t0 was rotated by -alpha
       
       
       corrected=imrotate(translated,median(relevant_angles),'crop');      
       corrected_pro=uint16(sum(corrected,3)/MetaData.Num_of_Pixels_Z);
       
       corrected_name = [filename,'_correct_t',num2str(iT),'.tif'];
       corrected_pro_name = [filename,'_correct2d_t',num2str(iT),'.tif'];
       
%         %saving 3d stacks
%         imwrite(corrected(:,:,1),corrected_name, 'WriteMode','overwrite');
%         for count=2:size(corrected,3)
%         imwrite(corrected(:,:,count),corrected_name,'WriteMode','append');
%         end
        %saving projected image
        imwrite(corrected_pro,corrected_pro_name, 'WriteMode','overwrite');
       
       
       if iT==1
           rotated_outline(:,:,1)=edge(imrotate(translated_nuc_pro,median(relevant_angles),'crop')>0,'sobel');           
       else
           rotated_outline(:,:,size(rotated_outline,3)+1)=edge(imrotate(translated_nuc_pro,median(relevant_angles),'crop')>0,'sobel');
       end

       
       
       if correction_type == 'pre'
           nuc_t0=nuc_t1; %sets the current image as the reference for the next loop
       end
     
       
        
       %% PERFORM ANALYSIS ABOVE!!!            
        end
        
        imtool(sum(rotated_outline,3))
        

    end

end

figure
for plot_count=1:size(nuc_rot,2)
    plot(1:size(nuc_rot,1),nuc_rot(:,plot_count))
    hold on
end
hold off


nuc_rotation=vertcat(rotation_header,num2cell(nuc_rot));
nuc_centroid=vertcat(centroid_header,num2cell(nuc_cen));
xlswrite('nuc_rotation.xls',nuc_rotation);
xlswrite('nuc_centroid.xls',nuc_centroid);
% save('nuc_rotation_rec.mat','nuc_rotation');
% save('nuc_centroid_rec.mat','nuc_centroid');
