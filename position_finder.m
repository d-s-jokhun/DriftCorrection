
%%% Change correction_type (default is ini), nuc_ch, max_rot, rot_step_size, files, the first for loop and bckgrnd_multiple as from line 93
%%% written by D.S.JOKHUN on 25/06/2016


clear all
clc
nuc_ch=1   %set the channel which will be used to mark the nucleus

centroid_header = {};
rotation_header = {};
nuc_cen=[];
cell_no=0;


name = 'rec_4FPS_004'
files = dir ([name,'*.nd2']);

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
%     MetaData.Plane_Origin_X = double(OmeMeta.getPlanePositionX(0,0).value);
%     MetaData.Plane_Origin_Y = double(OmeMeta.getPlanePositionY(0,0).value);
%     MetaData.Plane_Origin_Z = double(OmeMeta.getPlanePositionZ(0,0).value);
MetaData.Plane_Origin_X = 0;
MetaData.Plane_Origin_Y = 0;
MetaData.Plane_Origin_Z = 0;




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
        
        cell_no=cell_no+1;
        centroid_header{1,(2*cell_no)-1} = filename;
        centroid_header{1,2*cell_no} = filename;
        
        
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
           bckgrnd_multiple=1.6;
       end
       if f==9
           bckgrnd_multiple=1.1;
       end
       if f==10
           bckgrnd_multiple=1.2;
       end

       
       for count = 1:MetaData.Num_of_Pixels_Z;
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
                
                translated_nuc_pro = imtranslate(real_nuc_pro,[-(stats_cen.Centroid(1,1)-(size(real_nuc_pro,2)/2)), -(stats_cen.Centroid(1,2)-(size(real_nuc_pro,1)/2))],'FillValues',0);
                if iT==1
                    raw_outline (:,:,1)= edge(real_nuc_pro>0,'sobel');
                    translated_outline (:,:,1) = edge(translated_nuc_pro>0,'sobel');
                else
                    raw_outline (:,:,size(raw_outline,3)+1)= edge(real_nuc_pro>0,'sobel');
                    translated_outline (:,:,size(translated_outline,3)+1) = edge(translated_nuc_pro>0,'sobel');
                end


        
       %% PERFORM ANALYSIS ABOVE!!!            
        end
        
        imtool(flip(sum(translated_outline,3),1));  %flipped along the horizontal axis to match imaris
        
        rawname = [filename,'raw.tif'];
        translname = [filename,'translated.tif'];
        
%         imwrite(flip(raw_outline(:,:,1),1),rawname, 'WriteMode','overwrite');
        imwrite(flip(translated_outline(:,:,1),1),translname, 'WriteMode','overwrite');
        for count=2:size(raw_outline,3)
%         imwrite(flip(raw_outline(:,:,count),1),rawname,'WriteMode','append');
        imwrite(flip(translated_outline(:,:,count),1),translname,'WriteMode','append');
        end
        

        
    end

end


nuc_centroid=vertcat(centroid_header,num2cell(nuc_cen));
xlswrite(['nuc_centroid_',name,'.xls'],nuc_centroid);
% save(['nuc_centroid_',name,'.mat'],'nuc_centroid');
