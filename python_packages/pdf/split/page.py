import sys
import fnmatch
import PyPDF2 
import os

def page(str1,inp):
    pdfs = []
    if len(inp)==1:
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, inp[0]):
                pdfs.append(file)
    else:
        pdfs = inp

    for pdf in pdfs:
        mypdf = PyPDF2.PdfFileReader(pdf)
        
       
        out_pdf = PyPDF2.PdfFileWriter()
        list1 = str1.split(',')
        pages =[]
        for j in list1:
            if ':' in j:
                list2 = j.split(':')
                pages += list( range( int(list2[0])-1 , int(list2[1]) ) )
            else:
                pages.append(int(j)-1)
                

        for i in pages:
            try:
                out_pdf.addPage(mypdf.getPage(i))
            except:
                pass
        
        with open('Pages_'+pdf, 'wb') as fout:
            out_pdf.write(fout)
        out_pdf = None
            
        fout.close()    
    
if __name__ == '__main__':
    page(sys.argv[1],sys.argv[2:])
