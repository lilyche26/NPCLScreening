#!/usr/python
# -*- coding: UTF-8 -*-
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import glob
import os ,sys
import score_cal_rev
import glob,sys,getopt
import os,shutil
import subprocess

### python 3.Auto_exclusion_introns.py -q CAD_contigs.fasta

def query_dic_make(query_file):
    query_dic = {}
    records = SeqIO.parse(query_file,'fasta')
    for record in records:
        query_dic[str(record.description)]=str(record.seq)
    return query_dic
    print 'query_dic Done!'

def ref_dic_make(DB_file):
    ref_dic = {}
    records2 = SeqIO.parse(DB_file,'fasta')
    for record2 in records2:
        ref_dic[str(record2.description)]=str(record2.seq)
    return ref_dic
    print 'ref_dic Done!'


def blastx (query_file,DB_file,xml_file):
    command1 = 'makeblastdb.exe -in '+DB_file+' -dbtype prot\n'
    command2 = 'blastx.exe  -db  '+DB_file+' -query '+\
              query_file+'  -outfmt 5  -evalue 1E-10 -out '+xml_file+'\n'
    os.system(command1)
    os.system(command2)
    print 'blastx Done!'


def pair_identity(a,b):
    count_over=0
    for ii in zip(a,b):
        if ii[0]==ii[1]:
           count_over=count_over+1
    if len(a)!=0:
        over_identity=float(count_over)/len(a)
    elif len(a)==0 and count_over==0:
        over_identity=float(count_over)/1
    elif len(b)==0 and count_over==0:
        over_identity=float(count_over)/1
    return over_identity

def exclusion_introns(xml_file,out_file):
    
    query_dic = query_dic_make(query_file)
    ref_dic = ref_dic_make(DB_file)
    
    xml_handle=open(xml_file)
    blast_records = NCBIXML.parse(xml_handle)
    fff=open(out_file,'a')
    for blast_record  in blast_records:
        hit_num=0
        query = str(blast_record.query)
        query_len = str(blast_record.query_length)                
        for alignment in blast_record.alignments:
            hit_num=hit_num+1
            query_geneName=query.split('_')[0]

            sbjct_name=str(alignment.title.split(' ')[1])
            ali_len= alignment.length
            ref_cover_list=[]
            identity_list=[]
            Hsp_start_list=[]
            Ref_start_list=[]
        
            Hsp_end_list=[]
            Ref_end_list=[]

            Que_site_list=[]
            Ref_site_list=[]
        
            Hsp_sort_end_list=[]
            #Ref_sort_end_list=[]
            new_seq_list=[] ###!!!!!
            Intron_site_list = []

            Hsp_ref_list=[]
            Hsp_que_list=[]

            Query_list=[]
            Ref_list=[]
            
            temp=0
            #temp1=1
            pointer =0
            K=0.9##differ
            count_Iden=0
            for HSP in alignment.hsps:
                pointer=pointer+1
                start1=int(HSP.query_start )
                end1=int(HSP.query_end)
                start2=int(HSP.sbjct_start )
                end2=int(HSP.sbjct_end)
                Q= str(HSP.query)
                S=str(HSP.sbjct)
                seq_name=query
                identity=float(HSP.identities)/HSP.align_length
                identity_list.append(identity)
                
                ref_cover=(end2-start2)/float(ali_len)
                ref_cover_list.append(ref_cover)
                #print identity,ref_cover
                frame=int(HSP.frame[0])
                ppo=0

                if  frame>0 :
                    if  pointer==1:
                        ppo=1
                        S1=start2
                        S2=end2
                        Q1=start1
                        Q2=end1
                    else:
                        ##left
                        if  S1-start2>1 and S2-end2>0:
                            ds=S1-end2
                            dq=Q1-end1
                            if  3*ds<=dq*K and dq>0 and ds>0:##ref & query no overlap
                                ppo=1
                            elif dq<=0 and ds<=0 and 3*abs(ds)-abs(dq)*K>=0:##########query & ref overlap###intron=3*ref_overlap-query_overlap;
                                ppo=1
                                print 'query and ref overlap',query,xml_file
                            elif dq>=0 and ds<=0 and 3*abs(ds)+dq*K>=0:###intron=3*abs(ds)+dp#####intron=3*ref_overlap+query_overlap;
                                ppo=1
                        ##right
                        elif end2-S2>1 and start2-S1>0: ###right
                            #if start2>S2:###no overlap
                            ds=start2-S2
                            #elif start2<S2:###overlap
                            #ds=end2-S2
                            dq=start1-Q2
                            if 3*ds<=dq*K and dq>0 and ds>0:###ref & query no overlap
                               ppo=1
                            elif dq<=0 and ds<=0 and 3*abs(ds)>=abs(dq)*K:##########query & ref overlap######intron=3*ref_overlap-query_overlap
                                ppo=1
                                print 'query and ref overlap ',query,xml_file
                            elif dq>=0 and ds<=0 and 3*abs(ds)+dq*K>=0:###intron=3*abs(ds)+dp#####intron=3*ref_overlap+query_overlap;
                                ppo=1
                if  ppo==1:
                    if identity<0.5:
                       count_Iden=count_Iden+1
                    Hsp_start_list.append(start1)
                    Hsp_end_list.append(end1)
                
                    Ref_start_list.append(start2)
                    Ref_end_list.append(end2)
 
                        
                    Ref_site_list.append(start2)
                    Ref_site_list.append(end2)
                        
                    Que_site_list.append(start1)
                    Que_site_list.append(end1)

                    Hsp_ref_list.append([start2,end2])
                        
                    Hsp_que_list.append([start1,end1])

                    Query_list.append([(start1,end1,start2,end2),Q])
                    Ref_list.append([(start2,end2),S])

                    
            pin_count=str(len(Hsp_ref_list))
            print Hsp_ref_list,query.split('_')[1],pin_count,xml_file
            Hsp_sort_start_list=sorted(Hsp_start_list)
            Ref_sort_start_list=sorted(Ref_start_list)
            sorted_Ref_list=sorted(Ref_list)

            if pointer!= 1:###many HSPs
               many_HSP(Ref_sort_start_list,Ref_start_list,Ref_end_list,Ref_site_list,Que_site_list,Hsp_start_list,Hsp_end_list,Hsp_ref_list,\
                        Hsp_que_list,sorted_Ref_list,Query_list,new_seq_list,query_dic,ref_dic,seq_name,fff,sbjct_name,ali_len,count_Iden,pin_count,Intron_site_list) 
        #############a correct of frame translation，have intron,have gap################
            elif pointer == 1  and '*'in Q :
                 One_HSP_Stop_Codon (seq_name,start1,S,Q,query_dic,new_seq_list,fff,ali_len,pin_count,Intron_site_list)
            elif pointer==1 and  '*' not in Q and Ref_start_list!=[] :######one HSP
                 One_HSP_NO_Stop_Codon (query_dic,seq_name,Hsp_start_list,Hsp_end_list,fff,ali_len,pin_count)
                 
            break
    fff.close()
    print "Done!"




def many_HSP(Ref_sort_start_list,Ref_start_list,Ref_end_list,Ref_site_list,Que_site_list,Hsp_start_list,Hsp_end_list,Hsp_ref_list,\
             Hsp_que_list,sorted_Ref_list,Query_list,new_seq_list,query_dic,ref_dic,seq_name,fff,sbjct_name,ali_len,count_Iden,pin_count,Intron_site_list): 
   for key in query_dic:
      if key == seq_name:
         point=0
         overlap_list = []
         for sort in Ref_sort_start_list:
             for r in zip(Ref_start_list,Ref_end_list):
                     if sort==r[0]: 
                       
                        if point==0 :
                           sorted_1st_HSP(sort,key,Ref_site_list,Que_site_list,Hsp_start_list,Hsp_end_list,new_seq_list,query_dic,Intron_site_list)

                        if   point >0 and overlap_list!=[] and overlap_list[-1]<r[0] :##########no overlap case 1 2 
                             print seq_name.split(' ')[0],'no overlap',overlap_list[-1],r[0],sort
                             case1_no_overlap(sort,key,seq_name,r[0],overlap_list[-1],Ref_site_list,Que_site_list,Hsp_start_list,Hsp_end_list,new_seq_list,query_dic,Intron_site_list)
                             print 'case 1  :no_overlap_maybe_refine',new_seq_list,'new_seq_list'  
                        elif point>0 and overlap_list!=[] and overlap_list[-1]>=r[0] :###case 3 :temp>=r[0],overlap
                             print seq_name.split(' ')[0],'overlap',overlap_list[-1],r[0],sort
                             case3_overlap (sort,key,sbjct_name,overlap_list[-1],r[0],Hsp_ref_list,Hsp_que_list,sorted_Ref_list,Query_list,new_seq_list,query_dic,ref_dic,Intron_site_list)
                             

                        temp=r[1]
                        overlap_list.append(temp)          
                        point=point+1
                        break
         new_seq=('').join(new_seq_list)
         stop_count=str(str(Seq(new_seq).translate()).count('*'))
         missing=str(Seq(new_seq).translate()).count('X')
         cover=str(int((len(Seq(new_seq).translate())-missing)/float(ali_len)*100))
         if count_Iden<=0:
            seq_name=seq_name.split(' ')[0]+'_'+pin_count+'pin_'+stop_count+'*_'+cover+'cover'
         else:
            seq_name=seq_name.split(' ')[0]+'_'+pin_count+'pin_'+stop_count+'*_'+cover+'cover_'+str(count_Iden)+'_Ide0.5' 
         line1='>'+seq_name+'\n'+new_seq+'\n'
         fff.write(line1)

    
         

def sorted_1st_HSP(sort,key,Ref_site_list,Que_site_list,Hsp_start_list,Hsp_end_list,new_seq_list,query_dic,Intron_site_list):
    for rr in zip(Ref_site_list,Que_site_list):
        if sort ==rr[0]:
           x1=rr[1]
    for rrr in zip (Hsp_start_list,Hsp_end_list):
        if x1==rrr[0]:
           y=rrr[1]
           x=x1
           new_seq_list.append(query_dic[key][x-1:y])
           Intron_site_list.append(x-1)
           Intron_site_list.append(y)



def case1_no_overlap(sort,key,seq_name,START,END,Ref_site_list,Que_site_list,Hsp_start_list,Hsp_end_list,new_seq_list,query_dic,Intron_site_list):
    f1=open('no_overlap_may_refine_by_mannal.txt','a')
    cut=START-END
    if cut != 1:
       line=seq_name+'\t'+str(cut)+'AA\n'
       f1.write(line)
    f1.close()
    for rr in zip(Ref_site_list,Que_site_list):
           if sort ==rr[0]:
              x=rr[1]
    for rrr in zip (Hsp_start_list,Hsp_end_list):
           if x==rrr[0]:
              y=rrr[1]
              bu_gap=(cut-1)*3*'N'
              new_seq_list.append(bu_gap)
              new_seq_list.append(query_dic[key][x-1:y])
              print query_dic[key][x-1:y],x,y,'case1_no_overlap'
              Intron_site_list.append(x-1)
              Intron_site_list.append(y)

def case3_overlap (sort,key,sbjct_name,END,START,Hsp_ref_list,Hsp_que_list,sorted_Ref_list,Query_list,new_seq_list,query_dic,ref_dic,Intron_site_list):
    overlap=END-START
                                          
    if overlap>=0:
       for rr in zip(Hsp_ref_list,Hsp_que_list):
           if sort==rr[0][0]:
              ref_yb0=rr[0][0]
              yb0=rr[1][0]
              yb2=rr[1][1]##end
              yb1=yb0+(overlap+1)*3
           
       l=len(new_seq_list[-1])
       nuc_q1=str(new_seq_list[-1][l-(overlap+1)*3:])
       over_q1=Seq(str(new_seq_list[-1][l-(overlap+1)*3:])).translate()


       over_q2=Seq(query_dic[key][yb0-1:yb1-1]).translate()
       nuc_q2=query_dic[key][yb0-1:yb1-1]
       

       for ee in sorted_Ref_list:
           #print ee[0][1] 
           if START==int(ee[0][0]):
              if  overlap==0:
                  over_ref2=ee[1][0]
              else:
                  over_ref2=ee[1][0:overlap+1]
              for eee  in Query_list:
                  if  START==int(eee[0][2]):
                      if overlap==0:
                         over_qq2=eee[1][0]
                      else:
                         over_qq2=eee[1][0:(overlap+1)]
           if END==int(ee[0][1]):
              over_ref1=ee[1][-overlap-1:]
              for eee1  in Query_list:
                  if  END==int(eee1[0][3]):
                      over_qq1=eee1[1][-overlap-1:]

       b1=Seq(query_dic[key][yb1-1:yb2]).translate()

       a1=new_seq_list[-1][0:l-(overlap+1)*3]
       new_seq_list.pop()
       new_seq_list.append(a1)


       
       for key2 in ref_dic:
           if key2==sbjct_name:
              over_ref=ref_dic[key2][START-1:END]
       dic_identity={}
       tempp=0

       
       for ii in range(0,len(over_q1)+1):
           score_a=score_cal_rev.HSP_score(over_q1[0:ii],over_ref1[0:ii])
           score_b=score_cal_rev.HSP_score(over_q2[ii:],over_ref2[ii:])
           identity_all=float(score_a+score_b)/2
           identity_all=round(identity_all,1)
           if identity_all>=tempp:
              tempp=identity_all
           dic_identity.update({ii:identity_all})
            
       #print dic_identity,'dic_identity',tempp,'MAX'
       for key22 in dic_identity:
           if dic_identity[key22]==tempp:
                  j=int(key22)

       #############caculate two identity##########

                  
       list_true=[nuc_q1[0:j*3],nuc_q2[j*3:]]
       
       Intron_site_list.pop()
       xx=Intron_site_list[-1]+len(a1)+j*3
       Intron_site_list.append(xx)
       Intron_site_list.append(yb0-1+j*3)
       Intron_site_list.append(yb2)
       over_true=('').join(list_true)
    
       print Seq(over_true).translate(),'over_true',len(over_true),str(j)+'j',str(len(over_q1))+'over_len'
       new_seq_list.append(over_true)
       new_seq_list.append(query_dic[key][yb1-1:yb2])###b1


def One_HSP_Stop_Codon (seq_name,start1,S,Q,query_dic,new_seq_list,fff,ali_len,pin_count,Intron_site_list):
    list_site=[]
    list_site.append(0)
    for i in range(0,len(S)-1):
        if S[i]!='-' and S[i+1]=='-':
              tempp=i+1
              list_site.append(tempp)
        elif S[i]=='-' and S[i+1]!='-':
              tempp=i+1
              list_site.append(tempp)
        else:
            pass
    list_site.append(len(S))
    #print list_site
    #print S
    print seq_name,'right codon ;or seq quality is bad;or the relationship of phylogeny is far;'
    for key in query_dic:
              if key==seq_name:
                 for pp in range(0,len(list_site)-1,2):
                     m=start1+(list_site[pp]-0)*3
                     n=m+(list_site[pp+1]-list_site[pp])*3
                     new_seq_list.append(query_dic[key][m-1:n-1])
                     Intron_site_list.append(m-1)
                     Intron_site_list.append(n-1)
                     

                 new_seq=('').join(new_seq_list)         
                 stop_count=str(str(Seq(new_seq).translate()).count('*'))
                 missing=str(str(Seq(new_seq).translate()).count('?')) 
                 cover=str(int(len(Seq(new_seq).translate())/float(ali_len)*100))
                 seq_name=seq_name+'_'+pin_count+'pin_'+stop_count+'*_'+cover+'cover'
                 print seq_name
                 line2='>'+seq_name+'\n'+new_seq+'\n'
                 fff.write(line2)
        


def One_HSP_NO_Stop_Codon (query_dic,seq_name,Hsp_start_list,Hsp_end_list,fff,ali_len,pin_count):
    for key in query_dic:
        if key==seq_name:
           x=Hsp_start_list[0] 
           y=Hsp_end_list[0]
           new_seq=query_dic[key][x-1:y]
           stop_count=str(str(Seq(new_seq).translate()).count('*'))
           missing=str(Seq(new_seq).translate()).count('X')
           cover=str(int((len(Seq(new_seq).translate())-missing)/float(ali_len)*100))
           
           seq_name=seq_name.split(' ')[0]+'_'+pin_count+'pin_'+stop_count+'*_'+cover+'cover'
           print seq_name
           line='>'+seq_name+'\n'+query_dic[key][x-1:y]+'\n'
           fff.write(line)

                        

def run (query_file,DB_file,xml_file,out_file):
    blastx(query_file,DB_file,xml_file)
    exclusion_introns(xml_file,out_file)
    print 'All Done!'

    
if __name__ == "__main__":
    print "Exclusion introns..."
    try:
        print sys.argv[1:]
        opts, args = getopt.getopt(sys.argv[1:], "q:")
    except:
        sys.exit()
    for opt, arg in opts:
        if opt == "-q":
            query_file = os.path.abspath(arg)           
        else:
            pass
    xml_file='blastx_results.xml'
    DB_file='Tribolium_ref_pro_DB_beetles_final.fasta'
    out_file='results_exclusion_introns.fas'
    run(query_file,DB_file,xml_file,out_file)
    print "OK"



