info = niftiinfo('BF_1_t_nii.nii');
v = niftiread('BF_1_t_nii.nii');

for z= 1:91
    for x = 1:91 
        for y = 1:109
            if v(x,y,z)>0 && v(x,y,z-1)>0 && v(x,y,z-2)>0
            v(x,y,z) = (v(x,y,z)+ v(x,y,z-1)+ v(x,y,z-2))/3;
            v(x,y,z)
            end 
        end 
    end 
end 

niftiwrite(v,'BF_smoothed.nii',info)