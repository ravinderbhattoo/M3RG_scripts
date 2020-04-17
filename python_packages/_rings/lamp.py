r"""
This is python package for R.I.N.G.S software.

"""


def lamp2ring(inpf,outf,sfile,lfile):
    r"""Covert lammps coordinates output file to R.I.N.G.S software compatible.
    Replace atom number with symbols. (for e.g 1 -> C).
    Scale xs,ys,zs if neccessory.

    -------------
    Parametrs:
    inf    : lamps output file
    sfile  : scale parameter file
    lfile  : label file
    -------------

    Typical file format
    sfile:

    	#scale
	200         #x
	200         #y
	2           #z
	#end of file

    lfile:

     	#label add labels from line no. 2. Make a single column only.  #-> initiate comments
        C            #type1
        endoffile    #do not edit this line


        """
    import numpy as np

    x,y,z=np.loadtxt(inpf,skiprows=1,usecols=(2,3,4),unpack=True)
    tp=np.loadtxt(inpf,skiprows=1,usecols=(1))
    scale=np.loadtxt(sfile,skiprows=1)
    label=np.loadtxt(lfile,skiprows=1,dtype=str)

    tps=[]
    for ind,item in enumerate(tp):
        try:
            tps.append(label[int(item)-1])
        except:
            print('check label file!!!')

    x=x*scale[0]
    y=y*scale[1]
    z=z*scale[2]

    with open(outf,'w+') as f:
        for i,item in enumerate(tps):
            a=(tps[i])
            b=('%.3f' % x[i])
            c=('%.3f' % y[i])
            d=('%.3f' % z[i])
            f.write(a+' '+b+' '+c+' '+d+'\n')
        print('Data stored in : ',outf)
    f.close
