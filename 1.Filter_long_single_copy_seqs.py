#!/usr/python
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import glob,sys,getopt
import os,shutil
import subprocess


##python 1.Filter_long_single_copy_seqs.py  -Q Tribolium_castaneum_exon_0414_bu.fas -D Tribolium_castaneum.Tcas3.21.dna.toplevel.fa -e 1E-10 -L 600 -C 0.30 -S 0.50######
##out_file='single_copy_long_result.fas'
##xml_out_file='blastn_results.xml'

def blastn(query_file,DB_file,E_value,xml_file):
    command1 = 'makeblastdb.exe -in '+DB_file+' -dbtype nucl\n'
    
    command2 = 'blastn.exe -task blastn -db  '+DB_file+' -query '+\
              query_file+'  -outfmt 5  -evalue 1E-10 -out '+xml_file+'\n'
    os.system(command1)
    os.system(command2)
    
    
def single_copy_long_blast_filter(xml_file,Length,Cover,Similarity,out_file,dic):
    xml_handle=open(xml_file)
    blast_records = NCBIXML.parse(xml_handle)
    print xml_file
    out_handle=open(out_file,'w')
    for blast_record  in blast_records:        
        query = str(blast_record.query)
        query_len = str(blast_record.query_length)
        pointer =0
        q=0
        
        for alignment in blast_record.alignments:
            sbjct_name=alignment.title
            ali_len= str(alignment.length)
            pointer=pointer+1
            for HSP in alignment.hsps:
                seq_name=query
                if pointer == 1:
                   first_identify=float(HSP.identities)/HSP.align_length
                   first_cover=float(HSP.align_length)/float(query_len)
                   #print first_identify,first_cover
                   break
                if pointer == 2:
                   second_identify=float(HSP.identities)/HSP.align_length
                   second_cover=float(HSP.align_length)/float(query_len)
                   break
        if pointer ==1 :
           record='>'+str(query)+'\n'+dic[query]+'\n'
           q=1
        elif pointer >1:
            if second_identify >=float(Similarity) and second_cover>=float(Cover):
               pass
            else:
               record='>'+str(query)+'\n'+dic[query]+'\n'
               q=1
        else:
            pass
        print q,Length,len(dic[query])
        if q == 1  and len(dic[query]) >=int(Length):
            out_handle.write(record)
            q=q+1 
    #print q            
    xml_handle.close()
    out_handle.close()
    print 'single-copy and long filter done!'    



def run (query_file,DB_file,E_value,xml_file,Length,Cover,Similarity,out_file):
    dic={}
    records = SeqIO.parse(query_file,'fasta')
    for record in records:
        dic.update({str(record.id):str(record.seq)})
    print 'dic Done!'

    blastn(query_file,DB_file,E_value,xml_file)
    single_copy_long_blast_filter(xml_file,Length,Cover,Similarity,out_file,dic)

    
if __name__ == "__main__":
    print "Filter single-copy long seqs..."
    try:
        print sys.argv[1:]
        opts, args = getopt.getopt(sys.argv[1:], "Q:D:e:L:C:S:")
    except:
        sys.exit()
    for opt, arg in opts:
        if opt == "-Q":
            query_file = os.path.abspath(arg)
        if opt == "-D":
            DB_file = os.path.abspath(arg)
        if opt == "-e":
            E_value = arg
        elif opt == "-L":
            Length = arg
        elif opt == "-C":
            Cover= arg
        elif opt == "-S":
            Similarity= arg
        else:
            pass
    out_file='single_copy_long_result.fas'
    xml_file='blastn_results.xml'
    run(query_file,DB_file,E_value,xml_file,Length,Cover,Similarity,out_file)
    print "OK"
