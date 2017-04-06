from __future__ import division
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import urllib2
import re
import collections
import sys
import requests

with open("rosalind_dna.txt","r") as file:
    s=file.read()

Rproteins ={"UUU":"F", "UUC":"F", "UUA":"L","UUG":"L",
	    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
	    "UAU":"Y","UAC":"Y","UGU":"C","UGC":"C",
	    "UGG":"W","CUU":"L","CUC":"L","CUA":"L",
	    "CUG":"L","CCU":"P","CCC":"P","CCA":"P",
	    "CCG":"P","CAU":"H","CAC":"H","CAA":"Q",
	    "CAG":"Q","CGU":"R","CGC":"R","CGA":"R",
	    "CGG":"R","AUU":"I","AUC":"I","AUA":"I",
	    "AUG":"M","ACU":"T","ACC":"T","ACA":"T",
	    "ACG":"T","AAU":"N","AAC":"N","AAA":"K",
	    "AAG":"K","AGU":"S","AGC":"S","AGA":"R",
	    "AGG":"R","GUU":"V","GUC":"V","GUA":"V",
	    "GUG":"V","GCU":"A","GCC":"A","GCA":"A",
	    "GCG":"A","GAU":"D","GAC":"D","GAA":"E",
	    "GAG":"E","GGU":"G","GGC":"G","GGA":"G",
	    "GGG":"G"}
stop_list=["UGA","UAG","UAA"]

def countides():
  global s

  bases ={'A':0,'C':0,'G':0,'T':0}
  
  for k in s:
    if k=='A':
      bases['A']+=1
    elif k=='C':
      bases['C']+=1
    elif k=='G':
      bases['G']+=1
    elif k =='T':
      bases['T']+=1

  return bases


#print countides()




def DNA2RNA (s):
  

  s1=""

  for k in s:
    if k=='T':
      s1 = s1 + "U"
    else: s1 = s1 + k
  
  return s1


def CompDNA():
  global s
  
  l2=[]
  
  for k in s:
    if k=='T':
      l2.append('A')
    elif k=='C':
      l2.append('G')
    elif k=='G':
      l2.append('C')
    elif k=='A':
      l2.append('T')
    
  l2=l2[::-1]
  
  
  return ''.join(l2)


#print CompDNA()


def RNA2Prot(s):
  global Rproteins
  s2=""
  #s="AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
  
  for i in range(0,len(s)//3):
    if s[3*i:3*i+3] in stop_list:
      break
    s2 = s2 + Rproteins[s[3*i:3*i+3]]
  
  return s2

#print Rproteins.values()

#print RNA2Prot()
  
def mRNA_from_Prot():
  global s, Rproteins
  no_of_mrna=1
  for k in s:
    if k in Rproteins.values():
      no_of_mrna= no_of_mrna* Rproteins.values().count(k)
#print k
  return (no_of_mrna*3)%1000000
#because of the stop codon x3 modulo 1M
#print mRNA_from_Prot()
  



def Mendel_for_Dom():
  HomD=21
  Het=21
  Rec=21
  Tot= HomD+Het+Rec
  
  Pr_het= Het*(Het-1)/(Tot*(Tot-1)) 
  Pr_rec= Rec*(Rec-1)/(Tot*(Tot-1))
  Pr_rec_het = (2*Rec*Het/(Tot*(Tot-1)))
  
  Pr= Pr_het/4 +Pr_rec+ Pr_rec_het/2
  
  
  P= 1 -Pr

  return P

#print Mendel_for_Dom()
  
  
def DNAmotif(s1, s2):
  #ss="TGTATTTAGTACTTGAGTATTTAGTATTTATCACGCGTATTTACCTAGTATTTAAGTATTTAGCCTGTATTTACCGTATTTAAGCATCACTCAAATAGGGTATTTAGAGTGTATTTAGTATTTATAGTATTTAGTATTTAGTATTTACGTATTTATGTACGGTATTTAATGGTATTTAATATTACAGACGGAAGTATTTAGTGTCTGTATTTATGTATTTACAATTACCGGTATTTAGGTATTTACAGGAGTCTGGTATTTAAGGTATTTAGTATTTAAGGTAGTATTTAGTATTTAACGTATTTACGCTTTGTATTTATGTATTTACGTATTTACCAGGGTATTTAGTATTTATCACATGTATTTAGTATTTAAGCCGACGTATTTACGAACGCCAGTATTTATCGTATTTACCGGGTATTTAGACCCTGTATTTACCACCGTATTTAAGTATTTAGTATTTAGGTATTTATGTAGTATTTAGTGTATTTAAGTTAGGCATCGTATTTAGTATTTACCGTATTTATCGTATTTAGCGTATTTAGTATTTAAACTGTATTTAGTATTTAGTAGTTCAAGGTATTTACATCGAACCGGTATTTATGTATTTACGTATTTAGTATTTAGTGTACCTGTATTTAGTATTTAGACCTGGGTATTTAGTATTTAGGGAGGTATTTAGGTATTTACCTGGTATTTAATGTATTTACGTATTTAGTACGCGGCAAGTATTTAGTATTTATCGATTGCTGTATTTAGGTGTATTTATGTATTTAATAGTATTTAACGCATCCTCGTATTTAGAGTATTTAGTATTTAGTATTTACCGGTATTTATTGTGCGTGTATTTACGCGCCCTGTATTTAATAGCAGTATTTACGGTATTTAGTAGTATTTAGGCCGCGGGAGTATTTACGTATTTAGTATTTAGTATTTAAGTGTATTTA GTATTTAGT"
  
  #s1=ss.split()[0]
  #s2=ss.split()[1]
  n=len(s1)
  l1= []
  j=0
  
  for k in s1:
    
    j=s1.find(s2,j+1)
    if j<0:
      return l1
    l1.append(j+1)
      
  

def ProtMotif(filenm):
  f = open(filenm,'r')
  s =  f.read().strip().split('\n')
 
#print entradas
 
  for url in s:
      r = requests.get('http://www.uniprot.org/uniprot/%s.fasta' % url)
      prot = r.text
      newline_ind = prot.find('\n')
      prot = prot[(newline_ind+1):]
      prot = prot.replace('\n','')
      pos = re.finditer('(?=(N[^P][ST][^P]))',prot)
      p2=[]
      for i in pos:
	  p2.append(i.start()+1)
      if len(p2)>0:
	  print url
	  for i in p2:
	      print i,
	  print
  #with open("rosalind_mprt.txt","r") as file:
    #urls=file.read().strip().split('\n')
  
  #print urls
  
  #zeros=15*[0]
  #urls.sort()
  #result = dict(zip(urls,zeros))
  
  #for u in urls:
    
    #link="http://www.uniprot.org/uniprot/"+u+".fasta"
    #f = urllib2.urlopen(link)
    #myfile = f.read()
    
    #beg=myfile.find('\n')
    #s1=myfile[beg+1:]
    #print s1
    #p2=[]
    #protseq=''.join(s1)
    #p=re.compile("(?=(N[^P][ST][^P]))")
    
    #for m in p.finditer(protseq):
      #p2.append(m.start())
    #result[u]=p2
  
  #rs2 =collections.OrderedDict(sorted(result.items()))
  #return rs2

#rec=ProtMotif()
#for r in rec:
  #if rec[r]!=[]:
    #print r
    #for k in rec[r]:
      #print k
    


def CountPointMut():
  global s
  sl=s.split('\n')[:2]
  s1=sl[0]
  s2=sl[1]
  r=0
  
  for i in range(0,len(s1)):
    if s1[i]!=s2[i]:
      r=r+1
      
  return r

#print CountPointMut()
  

def SharedMotif():
  global s
  sl=re.split('>Rosalind_\d\d\d\d',s)[:len(s)]
  
  for st in sl:
    sl[sl.index(st)]= st.replace('\n','')
  
  sl.sort(key = len)
  
  #print sl
  
  check=False
  
  sls=sl[1:]
  s1=sls[0]
  for i in range(len(s1),1,-1):
    for l in range(0,len(s1)-i):
      sample= s1[l:l+i]
      #print sample
      for sk in sls[1:]:
	if sk.find(sample)!=-1:
	  check=True
	elif sk.find(sample)==-1: 
	  check = False
	  break
      if (check==True):
	return sample
	
#print SharedMotif()
  
def RNASplice():
  global s
  sl=re.split('>Rosalind_\d\d\d\d',s)[1:len(s)]
  for st in sl:
    sl[sl.index(st)]= st.replace('\n','')
  #print sl
  dna=sl[0]
  introns=sl[1:]
  
  for i in introns:
    dna=dna.replace(i,'')
  rna=DNA2RNA(dna)

  #print rna
  prot=RNA2Prot(rna)
  
  return prot
 
#print RNASplice()



ProtMotif("rosalind_mprt.txt")