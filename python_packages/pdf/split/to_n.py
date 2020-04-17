import sys
import fnmatch
import PyPDF2 
import os

def to_n(n,inp):
    pdfs = []
    if len(inp)==1:
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, inp[0]):
                pdfs.append(file)
    else:
        pdfs = inp
        
    for pdf in pdfs:
        mypdf = PyPDF2.PdfFileReader(pdf)
        tnp = mypdf.getNumPages()

        step = int(tnp/n) + 1*(tnp/n > int(tnp/n) + 0.0000000000001)

        for j in range(n):
            out_pdf = PyPDF2.PdfFileWriter()
            for i in range(step):
                try:
                    out_pdf.addPage(mypdf.getPage(j*step+i))
                except:
                    pass
        
            with open(str(j)+'_'+pdf, 'wb') as fout:
                out_pdf.write(fout)
            out_pdf = None
            
            fout.close()    
    
if __name__ == '__main__':
    to_n(int(sys.argv[1]),sys.argv[2:])
