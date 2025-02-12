%A funtion to analyze ASF vertex Produced interferograms and convert them
%to a value of apparent strain
function Landfast_Ice_App_Strain_calc(region,file)
%If you have mulitple regions thedefault path should get you to the folders
%of the different regions
if ~exist('region', 'var')
    defaultpath = 'G:/My Drive/Stability/RawData';
    spath = uigetdir(defaultpath, ...
                             'Select SLIE data directory');
    splitpath = regexp(spath, filesep, 'split');
    regdir = splitpath{end};
    region = regdir(1:4);
end
%Orgainze the ASF vertex into regions as some regions can have SAR pairs
%with simialr names based on the SAR pairs used
%In each SAR pair output folder the contexts of each folder should be taken out of the
%folder and put into these locations. 
switch region
        case 'EXWR'
            spath = 'G:/My Drive/Stability/RawData/Extra_Wainwright'; 
        case 'UTNQ'
            spath = 'G:/My Drive/Stability/RawData/UTQ-NSQ';
        case 'KZSO'
            spath = 'G:/My Drive/Stability/RawData/Kotzebue_Sound';  
        case 'OTKZ'
            spath = 'G:/My Drive/Stability/RawData/Outer_Kotzebue';                
        case 'WRUT'
            spath = 'G:/My Drive/Stability/RawData/Wainwright-UTQ';    
        case 'PLWR'
            spath = 'G:/My Drive/Stability/RawData/Point_Lay-Wainwright';              
        case 'PBKT'
            spath = 'G:/My Drive/Stability/RawData/Prudhoe_Bay-Katovik';  
        case 'PHPL'
            spath = 'G:/My Drive/Stability/RawData/Pt_Hope-Pt_Lay';                          
        case 'ULSP'
            spath = 'G:/My Drive/Stability/RawData/Uelen-SewardPenn';
        case 'NQKA'
            spath = 'G:/My Drive/Stability/RawData/NSQ-KAK';          
         case 'RUSS'
            spath = 'G:/My Drive/Stability/RawData/Russia';  
         case 'CODT'
            spath = 'G:/My Drive/Stability/RawData/Colville_Delta';
         case 'KAMK'
            spath = 'G:/My Drive/Stability/RawData/Kaktovik_MK_Delta';
         case 'NEKS'
            spath = 'G:/My Drive/Stability/RawData/NE_Kotzebue_Sound';
         case 'RURE'
            spath = 'G:/My Drive/Stability/RawData/Russia_redo';
         case 'AMSS'
            spath = 'G:/My Drive/Stability/RawData/AMSS';
        otherwise
            disp(['Region not recognized: ', region]);
            return
end



%Get a list of all files in the folder with the desired file name pattern.
%It is important that when ordering a sar pair from ASF vertex be sure to
%click the button requesting the wrapped interferograms 
filePattern = fullfile(spath, '*_wrapped_phase.tif'); 
%a list of all wrapped phase interferograms to process 
theFiles = dir(filePattern);

%read each wrapped interferogram as a raster image 
 for k = 1 : length(theFiles)
     baseFileName = theFiles(k).name;
     fullFileName = fullfile(theFiles(k).folder, baseFileName);
     fprintf(1, 'Now reading %s\n', fullFileName);
      %A is the data and R is the geographic information 
     [A,R] = readgeoraster(fullFileName);

%Raster layer infomration represents the phase information 
phase = A;
%Getting the projection information from the file geotiff information 
proj = geotiffinfo(fullFileName);  


%Read GeoTiff of Coherence file 
Corr_root_name = baseFileName(1:60);
suffix1 = '_corr.tif';
Corr = append(Corr_root_name,suffix1);
Corr_file = geotiffread([spath filesep ...
                       Corr]);
%Read in the theta file to account for line of sight motion vs hotizontal motion
theta_root_name = baseFileName(1:60);
suffix2 = '_lv_theta.tif';
theta = append(theta_root_name,suffix2);
theta_file = geotiffread([spath filesep ...
                       theta]);
%Set up the kernals in both directions to calcuate the gradients 
xgradient = [-1, 0, 0, 1];
ygradient = [-1; 0; 0; 1];


%Using principles of Libert et al to calcuate the phase gradient using the complex value of the
%phase difference to account for the phase wrap 
phase_gradient_x = conv2(phase(1:end,:), xgradient, 'same');
cmplx_gradient_x = exp(1j.*phase_gradient_x);
phase_dx = atan2(imag(cmplx_gradient_x),real(cmplx_gradient_x));

phase_gradient_y = conv2(phase(1:end,:), ygradient, 'same');
cmplx_gradient_y = exp(1j.*phase_gradient_y);
phase_dy = atan2(imag(cmplx_gradient_y),real(cmplx_gradient_y));
%Combine the gradient in the x and y direction to get the total magnitude
%of the phase gradient 
phase_gradient = (sqrt(phase_dx.^2 + phase_dy.^2))/3;



%Convert vaule of phase gradient to a real world value based on the length of the sailite wavelength (5.6 cm)
app_strain = phase_gradient*(0.0556/(4*pi));
%80 meter pixel size 


%Create Binary Image from the Coherence to use as a mask later 
%convolution matrix adds two pixels in both x and y direction not allowing
%a direct mask 
%morphologial filtering 
binary_Coher = imbinarize(Corr_file,0.1);
se = strel('disk',6); 
se2 = strel('disk',10);
clean = bwmorph(binary_Coher,'clean');
bigopen = bwareaopen(clean,75);
afterOpening01 = imopen(bigopen,se);
afterOpening02 = imopen(afterOpening01,se);
afterOpening03 = imopen(afterOpening02,se);
afterOpening04 = imopen(afterOpening03,se);
afterOpening05 = imopen(afterOpening04,se);
afterOpening06 = imopen(afterOpening05,se);
afterOpening07 = imopen(afterOpening06,se);
afterOpening08 = imopen(afterOpening07,se);
afterOpening09 = imopen(afterOpening08,se);
afterOpening10 = imopen(afterOpening09,se2);
afterOpening11 = imopen(afterOpening10,se2);
afterOpening12 = imopen(afterOpening11,se2);
afterOpening13 = imopen(afterOpening12,se2);
afterOpening14 = imopen(afterOpening13,se2);
afterOpening15 = imopen(afterOpening14,se2);

afterclosing01 = imclose(afterOpening15,se);
afterclosing02 = imclose(afterclosing01,se);
afterclosing03 = imclose(afterclosing02,se);
afterclosing04 = imclose(afterclosing03,se);
afterclosing05 = imclose(afterclosing04,se);
afterclosing06 = imclose(afterclosing05,se);
afterclosing07 = imclose(afterclosing06,se);
afterclosing08 = imclose(afterclosing07,se);
afterclosing09 = imclose(afterclosing08,se);
afterclosing10 = imclose(afterclosing09,se2);
afterclosing11 = imclose(afterclosing10,se2);
afterclosing12 = imclose(afterclosing11,se2);
afterclosing13 = imclose(afterclosing12,se2);
afterclosing14 = imclose(afterclosing13,se2);
afterclosing15 = imclose(afterclosing14,se2);
final_open = bwareaopen(afterclosing15,100);
%After morphological filtering select the biggest feature remaining to
%remove any not contiguious areas which were big enough to pass through the
%morphological filters
Biggest_Feature = bwareafilt(final_open,1);

%Zero the apparent strain values for all areas which do not meet the coherence requirments 
mask_app_strain = Biggest_Feature.*app_strain;


%set of the output directories 
datadir = fullfile(spath, 'Strain_Files');
if ~exist(datadir, 'dir')
   mkdir(datadir)
end
key = proj.GeoTIFFTags.GeoKeyDirectoryTag;
%retain originl name and add prefix denoting apparent strain
Apparent_strain_name = ['app_strain_' Landmask_root_name];
%Write the geotiff with oringinal name and geogrpahic infomation 
geotiffwrite(Apparent_strain_name,mask_app_strain,R,'GeoKeyDirectoryTag',key);

fprintf(1, 'Saving: %s\n', Apparent_strain_name);
end
