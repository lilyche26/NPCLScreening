from Bio import SeqIO
import glob,sys,getopt
import os,shutil



#####python 2.Auto_filter_alignments.py  -f Input_MSA -o Output_Filter_MSA -L 7 -I 0.80 -A 100

################for alignment conservation############
def ComputeIdentity(s):
    s=s.upper()
    temp_s = s
    majority = 0
    gap_num=s.count('-')
    if gap_num>2:
        return 0
    else:
        for i in range (0,len(s)):
            n = temp_s.count(s[i])
            temp_s = temp_s.replace(s[i],'')
            if majority < n:
                majority = n
            if temp_s == '':
               return float(majority)*1.0/len(s)
               break

###############abandon bad taxa sequence ,build new_records#########
def Auto_filter_alignment(MSA_list,f1,Block_Length,Expect_identity,Min_Apart):
    
    WindowWidth= int(Block_Length)
    Expect_identity= float(Expect_identity)
    MinLen= int(Min_Apart)

    flhandle=open(f1,'w')
    li=['Gene','No. of Taxa','Filter Results','Primer_Start','Primer_End','Apart_Length']
    line1='\t'.join(li)+'\n'
    flhandle.write(line1)

    for filename in MSA_list:
        pp=0
        records = list(SeqIO.parse(filename,'fasta'))
        new_records=[]
        new_records=records
        #print "Found %i records" % len(records)
                
        pointer = 0
        column =[]##read by column
        
        for new_record in new_records:
            #print str(new_record.id)        
            seq = str(new_record.seq)
            pointer = pointer +1
            for i in range (0,len(seq)):
                #if seq[i]!='-':
                    if pointer == 1 :
                       column.append(seq[i])
                    else:
                       column[i]=column[i]+seq[i]

        cut_start = 0
        cut_end = 1       
        #print len(column[0]),column[0]            
        ##########slide  window 
        for i in range (0,len(column)-WindowWidth):
            total_identity = 0
            for j in range (0,WindowWidth):               
                total_identity = total_identity + ComputeIdentity(column[i+j])
                #print ComputeIdentity(column[i+j]),i,'aaaaa'
            if total_identity/WindowWidth >= Expect_identity :
                cut_start = i+1
                break
        #print cut_start,float(total_identity)/WindowWidth,filename

        for i in range (len(column)-1,WindowWidth,-1):
            total_identity = 0
            for j in range (0,WindowWidth):           
                total_identity = total_identity + ComputeIdentity(column[i-j])
            if total_identity/WindowWidth >= Expect_identity:
                cut_end = i+1
                break
        length=cut_end - cut_start 
        identity=0
        if cut_end - cut_start >= MinLen:
            pp=1
            #print 'Yes!',length,cut_start,cut_end 
        else:
            pass
            #print 'NO',length,cut_start,cut_end

        if pp==1:
           #print filename,filename.split('\\')[-1],'dsagftty'
           out=open(outDir+'\\'+os.path.basename(filename),'w')
           SeqIO.write(new_records,out,'fasta')
           line=os.path.basename(filename)+'\t'+str(len(records))+'\tSuitable\t'+str(cut_start)+'\t'+str(cut_end)+'\t'+str(length)+'\n'
        else:
           line=os.path.basename(filename)+'\t'+str(len(records))+'\tNot Suitable\n'
           pass
           #print 'NO',filename
        flhandle.write(line)
    flhandle.close()

def run (Inputdir,outDir,f1,Block_Length,Expect_identity,Min_Apart):
    InputFilesLists=glob.glob(os.path.join(Inputdir,'*.fasta'))
    Auto_filter_alignment(InputFilesLists,f1,Block_Length ,Expect_identity, Min_Apart)

   
if __name__ == "__main__":
    print "Filter_alignments..."
    try:
        print sys.argv[1:]
        opts, args = getopt.getopt(sys.argv[1:], "f:o:L:I:A:")
    except:
        sys.exit()
    for opt, arg in opts:
        if opt == "-f":
            InputFilesLists = os.path.abspath(arg)
        elif opt == "-o":
            outDir = os.path.abspath(arg)
            if not os.path.isdir(outDir):
                os.makedirs(outDir)
        elif opt == "-L":
            Block_Length = arg
        elif opt == "-I":
            Expect_identity= arg
        elif opt == "-A":
            Min_Apart= arg
        else:
            pass
    f1='Filter_alignments_results.xls'
    run(InputFilesLists, outDir, f1,Block_Length , Expect_identity, Min_Apart)
    print "OK"
    
print "Complete"
