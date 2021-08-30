
%%% written by D.S.JOKHUN on 23/01/2017


clear all
clc

nuc_ch=1   %set the channel which will be used to mark the nucleus


centroid_header = {};
nuc_cen=[];
cell_no=0;


name = 'LMNA'
files = dir ([name,'*.nd2']);


for f=5%1:size(files,1);
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
        
                 
        nuc_edge_stack=[];
        nuc_pos_n = [];
        nuc_pos_n(1:MetaData.TimePoints,1:2)=0;
        
        cell_no=cell_no+1;
        centroid_header{1,(2*cell_no)-1} = filename;
        centroid_header{1,2*cell_no} = filename;
        

        
        for iT=1:50:MetaData.TimePoints;
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
           bckgrnd_multiple=2;
       end
       if f==2
           bckgrnd_multiple=1.2;
       end
       if f==3
           bckgrnd_multiple=3;
       end
       if f==4
           bckgrnd_multiple=3.5;
       end
       if f==5
           bckgrnd_multiple=6;
       end
       if f==6
           bckgrnd_multiple=3;
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
            if size(nuc_CC.PixelIdxList{obj_count},1)>1500
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
       nuc_edge_stack(:,:,iT)=nuc_edge;                                                 

       
       if iT==51
%            imtool(sum(real_nuc_bw,3))
           imtool(sum(XYZ{nuc_ch},3),[])
%            imtool(real_nuc_pro,[])
           imtool (uint16(sum(XYZ{nuc_ch},3)).*uint16(real_nuc_pro(51:50+size(real_nuc_bw,1),51:50+size(real_nuc_bw,2))),[])
       end
       if iT==451
%            imtool(sum(real_nuc_bw,3))
           imtool(sum(XYZ{nuc_ch},3),[])
%            imtool(real_nuc_pro,[])
           imtool (uint16(sum(XYZ{nuc_ch},3)).*uint16(real_nuc_pro(51:50+size(real_nuc_bw,1),51:50+size(real_nuc_bw,2))),[])
       end

        
       %% PERFORM ANALYSIS ABOVE!!!            
        end
        
        imtool(sum(nuc_edge_stack,3),[])


    end

end

