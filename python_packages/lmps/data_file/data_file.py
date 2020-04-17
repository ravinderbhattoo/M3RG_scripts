import math
from random import shuffle
import sys
import logo

def main():
    cs()
    
def cs():
    print('\nData file for core-shell type potential.\n')
    CF = input('Chemical Formula: ')
    first_line = 'Coreshell: '+CF+'\n\n'

    cs_nat = int(input('Number of CS atoms type?: '))
    cs_symb=[]
    cs_mass=[]
    cs_charge=[]
    cs_n=[]
    cs_m=[]
    cs_fn=[]
    for i in range(0,int(cs_nat)):
        s = input('Symbol for '+str(i+1)+'-CS atom: ')
        cs_symb.append(s)
        m1 = input('Mass for core of '+s+' atom: ')
        m2 = input('Mass for shell of '+s+' atom: ')
        c1 = input('Charge for core of '+s+' atom: ')
        c2 = input('Charge for shell of '+s+' atom: ')
        cs_mass.append([m1,m2])
        cs_charge.append([c1,c2])
        cn = input('Atom number for core of '+s+': ')
        sn = input('Atom number for shell of '+s+': ')
        cs_n.append([cn,sn])
        n = input('Molecule number of '+s+': ')
        cs_m.append(n)
        n = input('Number of atoms for '+s+' in complete formula: ')
        cs_fn.append(int(n))

    ncs_nat = int(input('Number of non-CS atoms type?: '))
    ncs_symb=[]
    ncs_mass=[]
    ncs_charge=[]
    ncs_n=[]
    ncs_m=[]
    ncs_fn=[]
    for i in range(0,int(ncs_nat)):
        s = input('Symbol for '+str(i+1)+'-non-CS atom: ')
        ncs_symb.append(s)
        m1 = input('Mass of '+s+' atom: ')
        c1 = input('Charge of '+s+' atom: ')
        ncs_mass.append(m1)
        ncs_charge.append(c1)
        n = input('Atom number of '+s+': ')
        ncs_n.append(n)
        n = input('Molecule number of '+s+': ')
        ncs_m.append(n)
        n = input('Number of atoms for '+s+' in complete formula: ')
        ncs_fn.append(int(n))

    nr = int(input('Number of replecation of base formula: '))
    num_str = str(int(nr)*(2*sum(cs_fn)+sum(ncs_fn)))+' atoms\n'+str(int(nr)*sum(cs_fn))+' bonds\n0 angles\n0 dihedrals\n\n' 
    typ_str = str(2*cs_nat+ncs_nat)+' atoms types\n'+str(cs_nat)+' bonds types \n0 angles types\n0 dihedrals types\n\n' 


    #add_masses
    ad_mass = 'Masses\n\n'
    str4=''
    for ind,i in enumerate(cs_n):
        str4 += i[0]+'\t'+cs_mass[ind][0]+'\n'
        str4 += i[1]+'\t'+cs_mass[ind][1]+'\n'
    for ind,i in enumerate(ncs_n):
        str4 += i+'\t'+ncs_mass[ind]+'\n'
        



    #add_atoms
    ad_atoms = '\nAtoms\n\n'
    total_cs_atoms = nr*cs_fn
    total_ncs_atoms = nr*ncs_fn
    total = sum(total_cs_atoms) + sum(total_ncs_atoms) 
    exp_bl = 1.6

    side = int(math.pow(total,1/3)+1)

    coords=[]
    for i in range(side):
        for j in range(side):
            for k in range(side):
                coords.append([str(exp_bl*k),str(exp_bl*j),str(exp_bl*i)])    


    xyz=[float(min(coords[:][0]))-0.1,float(max(coords[:][0]))+0.1,
         float(min(coords[:][1]))-0.1,float(max(coords[:][1]))+0.1,
         float(min(coords[:][2]))-0.1,float(max(coords[:][2]))+0.1]
    s=['xlo','xhi','ylo','yhi','zlo','zhi']
    dim_str = str(xyz[0])+' '+str(xyz[1])+' '+'xlo xhi\n'+str(xyz[2])+' '+str(xyz[3])+' '+'ylo yhi\n'+str(xyz[4])+' '+str(xyz[5])+' '+'zlo zhi\n\n'


    ind1 = [i for i in range(total)]
    shuffle(ind1)
    shift=0
    bshift=1
    str5=''
    str3=''
    str1=''
    for i in range(nr):
        for ind,atoms in enumerate(cs_fn):
            for j in range(atoms):
                str1 += (str(2*shift+1)+'\t'+cs_m[ind]+'\t'+cs_n[ind][0]+'\t'+cs_charge[ind][0]+'\t'+coords[ind1[shift]][0]+'\t'+coords[ind1[shift]][1]+'\t'+coords[ind1[shift]][2]+'\n')
                str1 += (str(2*shift+2)+'\t'+cs_m[ind]+'\t'+cs_n[ind][1]+'\t'+cs_charge[ind][1]+'\t'+coords[ind1[shift]][0]+'\t'+coords[ind1[shift]][1]+'\t'+coords[ind1[shift]][2]+'\n')
                str3 += (str(bshift)+'\t'+str(ind+1)+'\t'+str(2*shift+1)+'\t'+str(2*shift+2)+'\n')              
                str5 += (str(2*shift+1)+'\t'+str(shift+1)+'\n'+
                         str(2*shift+2)+'\t'+str(shift+1)+'\n')
                shift += 1
                bshift += 1
    str2=''
    nshift = 1*shift
    for i in range(nr):
        for ind,atoms in enumerate(ncs_fn):
            for j in range(atoms):
                str2 += (str(nshift+shift+1)+'\t'+ncs_m[ind]+'\t'+ncs_n[ind]+'\t'+ncs_charge[ind]+'\t'+coords[ind1[shift]][0]+'\t'+coords[ind1[shift]][1]+'\t'+coords[ind1[shift]][2]+'\n')
                shift += 1


    #add_bond
    ad_bonds = '\nBonds\n\n'
    #str3

    #add_info
    ad_info = '\nCS-info\n\n'
    #str5



    all_data = first_line+num_str+typ_str+dim_str+ad_mass+str4+ad_atoms+str1+str2+ad_bonds+str3+ad_info+str5

    f = open('./CS_datafile','w+')
    f.write(all_data)
    f.close()
    print('\n\n'+'\033[91m'+'\033[1m'+'Data stored in CS_datafile.'+'\033[0m'+'\n\n')


    logo.print_logo()
    









if __name__ == "__main__": main()

