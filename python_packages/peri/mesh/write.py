r"""
This is python package for Peridigm Mesh Generation(write).

Parameters:
------------------
outf      : Output data file (mesh) for Peridigm, /path_to_dir/filename
list_l1   : List of data
------------------
"""

import os

def write(outf,list_l1):
    total = 0
    for i in list_l1:
        total += len(i)

    try:
        os.makedirs(os.path.dirname(outf), exist_ok=Tru)
    except:
        pass

    init = 1

    for ind,i in enumerate(list_l1):
        with open(outf+'_list_'+str(ind+1),'w+') as f:
            for j in range(init,init+len(i)):
                       f.write(str(j)+'\n')
            init += len(i)
        f.close

    with open(outf,'w+') as f:
        f.write(str(total)+'\n')
        f.write('Properties=pos:R:3:species:S:1:Volume:R:1'+'\n')
        for i in list_l1:
            for j in i:
                f.write( str(j[0]) + ' ' + str(j[1]) + ' ' + str(j[2]) + ' ' + str(j[3]) + ' ' + str(j[4]) + '\n')
    f.close

    print('Data stored in: ',outf)
