#################################################################################################################
#                                                                                                               #
#   Software     :  DNSS (A Deep Learning Network Approach to ab initio Protein Secondary Structure Prediction) #
#   Release      :  1.1  (October 2016)                                                                         #
#                                                                                                               #
#   Author(s)    :  Matt Spencer, Jesse Eickholt, and Jianlin Cheng                                             #
#   Maintainance :  Jie Hou(jh7x3@mail.missouri.edu ), Badri Adhikari (bap54@mail.missouri.edu)                 #
#   Copyright    :  Bioinformatics, Data Mining, and Machine Learning Lab (BDM)                                 #
#                    Department of Computer Science                                                             #
#                    University of Missouri, Columbia                                                           #
#                                                                                                               #
#################################################################################################################

Note that probably every run will result in a warning from the blast 
software. This is not a problem, and SS predictions should work fine
despite this error.



----------------------------------------------------------------------------

1. Before installing DNSS, the following tools must be downloaded and installed

    1). Python2: DNSS was fully developed and tested under Python 2.7.6
        
        * Note: Python is recommended to install in system, in our program, will call 'python2 *py' directly.

    2). ncbi-blast-2.2.25+:  blast programs should be in folder './programs', for example, the binary file should be in 'programs/ncbi-blast-2.2.25+/bin/'

        * Note: We had hard code of 'ncbi-blast-2.2.25+' in 'scripts/gen-pssm-less-stringent.sh' and 'scripts/generate-pssm.sh'
                So if other version is used, please make sure to change the version in two scripts.
                Executable blast can be downloaed from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.25/
                For example: 
                
                $ cd ./programs
                $ wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.25/ncbi-blast-2.2.25+-x64-linux.tar.gz
                $ tar -zxvf ncbi-blast-2.2.25+-x64-linux.tar.gz

    3). Non-redundency database(nr90): Please put nr90 files that were formated by blast into folder './nr_database' 
        
        *Note: Please make sure the prefix of nr database is 'nr90', for example, nr90.00.phr, nr90.00.pin, etc
               If other version of nr database is used, please change the nr version in 'scripts/gen-pssm-less-stringent.sh' and 'scripts/generate-pssm.sh'
               Our lab provided nr70, nr90, and nr for users to download, please try 
               
               $ cd ./nr_database
               
                a) nr90: $ wget http://sysbio.rnet.missouri.edu/bdm_download/nr_database/nr90.tar.gz
                         $ tar -zxvf nr90.tar.gz
                         $ mv nr90/* ./
                         $ rm -rf nr90 nr90.tar.gz

                b) nr70(optional): 
                         $ wget http://sysbio.rnet.missouri.edu/bdm_download/nr_database/nr70.tar.gz
                         $ tar -zxvf nr70.tar.gz
                         $ mv nr70/* ./
                         $ rm -rf nr70 nr70.tar.gz
                         

                c) nr(optional):   
                         $ wget http://sysbio.rnet.missouri.edu/bdm_download/nr_database/nr.tar.gz
                         $ tar -zxvf nr.tar.gz
                         $ mv nr/* ./
                         $ rm -rf nr nr.tar.gz



----------------------------------------------------------------------------

2. Configuation

    Usage:
    $ cd DNSS
    $ perl configure.pl <DNSS folder>

    Example:
    $ perl configure.pl /home/jh7x3/DNSS/


----------------------------------------------------------------------------

3. Instructions for predicting secondary structure:
    * Enter the scripts/ directory.
    * Use the PredictSS.pl script to predict SS.
    * This tool also support both GPU and CPU to run prediction, the default is running on CPU. 
      Before running CPU, please check if GPU can be called for cudamat
      
      $ cd scripts/
      $ perl GPU_detect.pl  

    
    * There are two ways to indicate the protein to predict:
    ----------------------------------------------------------------------------
    (1) Directly give a sequence:
          You can indicate a residue sequence in the command line directly,
          and name the sequence to make resulting files easier to find.
          To use this procedure, it is reccomended that the -name tag is used
          in addition to the -seq tag. -name is not required, but if it is
          not used, the program will automatically assign a less useful name
          to the protein (such as 0000). 
          
          Usage:
          $  perl PredictSS.pl  -seq <AA sequence> -name <Protein Name> -device <GPU/CPU> -out   -out <output folder>
          
          Example:
          a). CPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -seq GNVVIEVDMANGWRGNASGSTSHSGITYSADGVTFAALGDGVGAVFDIARPTTLEDAVIAMVVNVSAEFKASEANLQIFAQLKEDWSKGEWDALAGSSELTADTDLTLTATIDEDDDKFNQTARDVQVGIQAKGTPAGTITIKSVTITLAQEA -name Prot1  -out ./output/
          b). GPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -device GPU -seq GNVVIEVDMANGWRGNASGSTSHSGITYSADGVTFAALGDGVGAVFDIARPTTLEDAVIAMVVNVSAEFKASEANLQIFAQLKEDWSKGEWDALAGSSELTADTDLTLTATIDEDDDKFNQTARDVQVGIQAKGTPAGTITIKSVTITLAQEA -name Prot1  -out ./output/

    ----------------------------------------------------------------------------
    (2) Predict from protein file:
          This is the suggested method for predicting single proteins, as it is
          less confusing. Create a fasta file containing the residue sequence
          and name the file <prot-name>.fasta. The program will use the name of 
          the file to indicate the protein that was predicted, and name output
          files accordingly.
          
          Usage:
          $ ./PredictSS.pl -seq <file name>.fasta -file -device <GPU/CPU>  -out <output folder>
          
          Example:
          a) CPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -seq test/1GNY-A.fasta -file -device CPU -out output/
          b) Gpu: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -seq test/1GNY-A.fasta -file -device GPU -out output/
          
          
          * Note that the 1GNY-A.fasta file is available at that location to use as
            a test protein.

          
          * Note that using the -out tag to designate an output directory is recommended,
            or else the scripts folder will get cluttered by predictions. (In fact, I
            recommend changing the default output directory to something else to avoid
            this fate.)
            
          * More details about the parameter, can simply type 
            $ perl scripts/PredictSS.pl

----------------------------------------------------------------------------

4. Predicting multiple proteins:
      It is quite easy to predict multiple proteins in sequence using the -indir
      tag. Just save fasta files of all the proteins you want to predict in a
      directory, use the -indir tag, and all fasta files in that directory will
      be predicted.
      
      Usage:
      $ ./PredictSS.pl -indir <input directory> -out <output directory> -device <GPU/CPU>
      
      Example:
      a) CPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -indir ./test/ -out ./output/ -device CPU
      b) GPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -indir ./test/ -out ./output/ -device GPU

----------------------------------------------------------------------------

5. Recommended improvements:
  As stated above, I recommend changing the default directory for output
    to avoid cluttering the scripts/ directory.
  Currently the program saves all intermediate files in the temp/ directory
    and its subdirectories. I recommend revising the script to delete these
    intermediate files, as they will build up to be quite large if left
    unattended.


----------------------------------------------------------------------------

If you have questions or suggestions, please contact:        
        Jie Hou(jh7x3@mail.missouri.edu), Jianlin Cheng, PhD
        Bioinformatics, Data Mining, and Machine Learning Lab (BDM)
        Department of Computer Science
        University of Missouri, Columbia
        Email: chengji@missouri.edu

----------------------------------------------------------------------------

Citation: Spencer, Matt, Jesse Eickholt, and Jianlin Cheng. "A deep learning network approach to ab initio protein secondary structure prediction." IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB) 12.1 (2015): 103-112.
