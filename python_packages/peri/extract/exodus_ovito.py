import sys
sys.path.append('/usr/local/trilinos/bin')
import os
from peridigm.extract.exodus_extract import extract
import numpy as np


inpf = sys.argv[1]
mydata = extract(inpf)

try:
     os.makedirs(mydata.mydir+'/ovito_files')
except:
     pass

if not(mydata.skip):
    mydir = mydata.mydir
    coords = mydata.coords
    time = mydata.time
    lines = len(mydata.coords)

    mydata.pos = {}
    mydata.pos['Position'] = []
    for i in range(len(time)):
        str_1 = 'Properties=species:S:1:pos:R:3:'
        print 'Step: ',i
        mydata.pos['Position'].append(coords+np.transpose(np.array((mydata.peri_var[0]['DisplacementX'][:,i],mydata.peri_var[0]['DisplacementY'][:,i],mydata.peri_var[0]['DisplacementZ'][:,i],))))
        out_array = []
        for var_name in sorted(mydata.peri_var[0]):
            try:
                a = mydata.peri_var[0][var_name][:,i]
            except:
                a = mydata.peri_var[0][var_name][i]
            size = 1
            try:
                size = a.shape[1]
            except:
                pass
            type = 'R'
            if var_name == 'Element_Id':
                var_name = 'id'
                type = 'I'
            elif var_name == 'species':
                var_name = 'species'
                type = 'S'
            else:
                pass

            str_1 += ("{}:{}:{}:".format(var_name,type,size))
            out_array.append(a)

        out_array = np.transpose(np.vstack(out_array))
        np.savetxt(mydir+"/ovito_files/"+inpf+"_"+str(i)+".xyz",(np.hstack((mydata.PI[:,None],mydata.pos['Position'][i],out_array))),comments='',header="{}\n{}".format(lines,str_1[:-1]))



    print '\n======================\n\x1b[5mAll data extracted in:\x1b[0m'
    print '\x1b[1m'+mydir+'\x1b[0m'+'\n\n'

    print '======================\n\x1b[5mAll ovito files extracted in:\x1b[0m'
    print '\x1b[1m'+mydir+"/ovito_files/"+'\x1b[0m'+'\n\n'


else:
    print '======================\n\x1b[5;1mEverything Skipped!!!!\x1b[0m'
