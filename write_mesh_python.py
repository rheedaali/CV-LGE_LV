#!/usr/bin/env python
from __future__ import print_function


from sys import argv
import numpy as np
import math
import collections


a=np.fromfile(argv[1],int)



pts=np.fromfile(argv[2])
# area_triangles=np.fromfile(argv[3])
pts=np.reshape(pts,(-1,3))
size_pts= pts.shape[0]
pts.astype(int)

print(pts)
b=np.reshape(a, (-1, 3))


# v_vec=np.fromfile(argv[4])
# v_vec=np.reshape(v_vec,(-1,3))
# size_v_vec= v_vec.shape[0]

v_scalars=np.fromfile(argv[3])


# filt=[~np.isinf(v_scalars)]

# v_vec= v_vec[filt]
# v_scalars=v_scalars[filt]

# b=b[filt]
# area_triangles=area_triangles[filt]

c= b.shape[0]
arrayof3s = np.empty(c)
arrayof3s.fill(3)

cells_info=np.column_stack((arrayof3s,b))
cells_info=np.reshape((cells_info), (-1, 4))
cells_info=cells_info.astype(int)
size_cells_info= cells_info.shape[0]



with open('target.vtk','w') as out:
    line1 = "# vtk DataFile Version 3.0 "
    line2 = "vtk output "
    line3 = "ASCII "
    line4 = "DATASET POLYDATA" 
    line5= "POINTS {} float".format(size_pts)
    print("I'm going to write these to the file.")
    out.write('{}\n{}\n{}\n{}\n{}\n'.format(line1,line2,line3,line4,line5))
#    import pdb; pdb.set_trace()
#this is already uplodaed ????
    np.savetxt(out, pts, delimiter=' ',fmt='%.7f')   # X is an array

    line6= "POLYGONS {} {}".format(size_cells_info, size_cells_info*4)
    out.write('{}\n'.format(line6))

    np.savetxt(out, cells_info, delimiter=' ',fmt='%d')   

    # line7="CELL_DATA {}".format(size_cells_info)
    # out.write('{}\n'.format(line7))
    # # line8="VECTORS fibre_orientations float"
    # # out.write('{}\n{}\n'.format(line7,line8))
    # # np.savetxt(out, v_vec, delimiter=' ',fmt='%.7f')   
    
    # # line9="SCALARS area float"
    # # line10="LOOKUP_TABLE default"
    # # out.write('{}\n{}\n'.format(line9,line10))
    # # np.savetxt(out,area_triangles)
    
    # line8="SCALARS CVs float"
    # line9="LOOKUP_TABLE default"
    # out.write('{}\n{}\n'.format(line8,line9))
    # np.savetxt(out,v_scalars)



