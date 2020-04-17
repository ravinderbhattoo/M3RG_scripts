import sys
import fnmatch
from PyPDF2 import PdfFileMerger
import os

def merge(inp):
    if len(inp)==1:
        pdfs = []    
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, inp):
                pdfs.append(file)
    else:
        pdfs = inp
        
    merger = PdfFileMerger()

    for pdf in pdfs:
        merger.append(open(pdf, 'rb'))

    with open('merged_document.pdf', 'wb') as fout:
        merger.write(fout)
    
if __name__ == '__main__':
    merge(sys.argv[1:])
