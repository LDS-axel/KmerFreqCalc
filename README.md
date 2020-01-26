KmerFreqCalc (KFC)    
==================================    
DATE:2020-1-24    
version:0.2.1


Author Information:    
==================     
Name:  Liu Shuai     
Department:  Anhui University and the Institute of Zoology, Chinese Academy of Sciences jointly cultivation     
E-mail:  ls2659614061@126.com


Supervisor Information:     
======================     
Name:  Wu Qi     
Department:  Institute of Microbiology, Chinese Academy of Sciences    
E-mail:  wuq@im.ac.cn


General Procedure     
=================     
* kmer frequency vector
    * 1.Read the sequences in the filement into a dictionary.
    * 2.Count the number of kmer in sequences.
    * 3.Make a summary statistic of all kmers' number of all subsequences.   
    * 4.Calculate the frequencies of kmers.
    * 5.Make a cv.txt formated output.  
    #example.cv.txt  
    K   
    Significant kmer number   
    Inner product    
    1 0.044418402253529    
    2 -0.13694127815637536    
    3 0.02302543931397616    
    4 0.08209917560394331     
    ......
* distance matrix
    * 1.Employee Cosine dissimilarity or Euclidean distance to calculate the difference among species
    * 2.Make a .meg format output opened up by MEGA.


Manual
========
kmer frequency vector
* python3 kmer_frequence_vector.py -K -fna.path -speciesname -CPU.NUM    
-K      -num    the length of kmer      K>=1    
-fna.path       -str    the path of dna sequence as an input,whose suffix is fna,fa,ffn,or fasta    
-tag.name    -str    the name of sequence     
-CPU.NUM        -num    the number of CPU be used

distance matrix
* python3 distance_matrix.py cvtxt.directory CPU.NUM distance-type(COS or EUC)     
cvtxt.directory     -str    The path of cv.txt of species stored      other suffix files must nonexistent     
CPU.NUM     -num            The number of CPU be used    
distance-type       -str    Euclidean distance or Cosine disimilarity
