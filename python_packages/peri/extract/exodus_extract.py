import sys
sys.path.append('/usr/local/trilinos/bin')
from peridigm.extract.exodus_v2 import exodus
import numpy as np
import os
import ntpath

class pyperi:
    """Documentation: a class"""

    def __init__(self,myname):
        self.name = myname
        self.peri_var = []

def my_check(if_file,check,overwrite):
    if if_file:
        if check:
            print('File already exits. ')
            overwrite = (str(raw_input('Overwrite?(y/n): '))=='y')
            check = (str(raw_input('Repeat for all?(y/n): '))=='n')
        else:
            pass
    return if_file and (overwrite==False),check,overwrite



def extract(inpf):
    skip = (str(raw_input("Skip calculations in case of no overwrite?(y/n)\nNo ovito files will be generated. : "))=='y')

    fpath = os.path.abspath(inpf)
    mydir = fpath + '_data'

    data = exodus(fpath, mode='r', array_type='numpy')
    try:
        os.makedirs(mydir)
    except:
        pass

    mydata = pyperi('Peridigm Data')

    mydata.skip = skip
    mydata.mydir = mydir
    mydata.time=data.get_times()
    time = mydata.time
    print 'Total time steps: ',len(time)
    print 'Extracting in: ',mydir



    overwrite = True
    check = True

    #node_variables
    print 'node_varibale'
    mydata.peri_var.append({})
    list1 = data.get_node_variable_names()

    print('===================================')
    for i in list1:
        print('\n--------------\nWritting '+i)
        if_file = os.path.isfile(mydir+'/'+i)

        bool1,check,overwrite = my_check(if_file,check,overwrite)
        if bool1 and skip:
            print(i+' skipped.\n')
        else:
            df=(data.get_node_variable_values(i,1))
            a = np.ndarray((len(df),len(time)))
            print('Size of array:')
            print(a.shape)
            a[:,0] = df
            for t in range(2,len(time)+1):
                a[:,t-1] = (data.get_node_variable_values(i,t))
            mydata.peri_var[0][i]=a
            if not(bool1):
                np.savetxt(mydir+'/'+i,a)
                print(i+' has been written.\n')
            else:
                print(i+' skipped.\n')


    #global_variables
    print 'global_variables'
    mydata.peri_var.append({})
    list1 = data.get_global_variable_names()
    for i in list1:
        print('\n--------------\nWritting '+i)
        if_file = os.path.isfile(mydir+'/'+i)

        bool1,check,overwrite = my_check(if_file,check,overwrite)
        if bool1 and skip:
            print(i+' skipped.\n')
        else:
            df=(data.get_global_variable_values(i))
            print('Size of array:')
            print(df.shape)
            mydata.peri_var[1][i]=df
            if not(bool1):
                np.savetxt(mydir+'/'+i,df)
                print(i+' has been written.\n')
            else:
                print(i+' skipped.\n')



    #element_property
    print 'element_property'
    list1 = data.get_element_variable_names()
    ids = data.get_elem_blk_ids()
    for i in list1:
        mydata.peri_var[0][i]=[]
        for j in ids:
            print('\n--------------\nWritting '+i+'_blk'+str(j))
            if_file = os.path.isfile(mydir+'/'+i+'_blk'+str(j))

            bool1,check,overwrite = my_check(if_file,check,overwrite)
            if bool1 and skip:
                print(i+' skipped.\n')
            else:
                df=data.get_element_variable_values(j,i,1)
                a = np.ndarray((len(df),len(time)))
                print('Size of array:')
                print(a.shape)
                a[:,0] = df
                for t in range(2,len(time)+1):
                    a[:,t-1] = data.get_element_variable_values(j,i,t)
                mydata.peri_var[0][i].append(a)
                if not(bool1):
                    np.savetxt(mydir+'/'+i+'_blk'+str(j),a)
                    print(i+'_blk'+str(j)+' has been written.\n')
                else:
                    print(i+' skipped.\n')

        if not(bool1 and skip):
            mydata.peri_var[0][i] = np.vstack(mydata.peri_var[0][i])


    #Coords
    i = 'Coords'
    print('\n--------------\nWritting '+i)
    if_file = os.path.isfile(mydir+'/'+i)
    bool1,check,overwrite = my_check(if_file,check,overwrite)

    if bool1 and skip:
        print(i+' skipped.\n')
    else:
        ids = data.get_elem_blk_ids()
        PI = []
        list1 = data.get_element_variable_names()
        for id in ids:
            PI += [id for j in range(len(data.get_element_variable_values(id,list1[0],1)))]
        mydata.coords = np.transpose(data.get_coords())
        mydata.PI = np.array(PI)
        if not(bool1):
            np.savetxt(mydir+'/'+i,mydata.coords)
            np.savetxt(mydir+'/'+'PI',mydata.PI)
            print(i+'and PI has been written.\n')
        else:
            print(i+' skipped.\n')


    np.save(mydir+'/'+'all_data',mydata)
    print '======================\n\x1b[5mAll data extracted in:\x1b[0m'
    print '\x1b[1m'+mydir+'\x1b[0m'+'\n\n'

    return mydata


#===========================================
if __name__ == '__main__':
    extract(sys.argv[1])
