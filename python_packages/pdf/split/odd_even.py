import sys
import fnmatch
import PyPDF2 
import os

def to_n(inp):
    pdfs = []
    if len(inp)==1:
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, inp[0]):
                pdfs.append(file)
    else:
        pdfs = inp

    label = ['odd','even']    
    for pdf in pdfs:
        mypdf = PyPDF2.PdfFileReader(pdf)
        tnp = mypdf.getNumPages()

        p = int(tnp/2)    
        
        for j in range(2):
            out_pdf = PyPDF2.PdfFileWriter()
            pages = [2*n+j for n in range(p+1-j)]
            for i in pages:
                try:
                    out_pdf.addPage(mypdf.getPage(i))
                except:
                    pass
        
            with open(label[j]+'_'+pdf, 'wb') as fout:
                out_pdf.write(fout)
            out_pdf = None
            
            fout.close()    
    
if __name__ == '__main__':
    to_n(sys.argv[1:])
