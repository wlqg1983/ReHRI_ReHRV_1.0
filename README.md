    In organelle genomes, especially mitochondrial genomes, there is a widespread phenomenon of genome recombination mediated by repetitive sequences (direct repeats and inverted repeats), which leads to the emergence of subconfigurations in the genomes. The ReHRI software can utilize second- and third-generation sequencing data (reads) to explore the potentially existing subconfigurations within the main configuration of the genome. Meanwhile, it calculates the probabilities of genome recombination mediated by different repetitive sequences based on the number of reads. The ReHRV software can draw genome maps before and after genome recombination to display the genomic structural variations.
    ReHRI and ReHRV need to run in an environment with at least Ubuntu (v20.04) and conda (v23.5.2). The tutorials for software installation and usage are as follows:

**1. Install the software**

    (1) Download and install the software from https://github.com/wlqg1983/ReHRI_ReHRV_1.0/archive/refs/heads/main.zip

    unzip ReHRI_ReHRV_1.0-main.zip
    cd  ReHRI_ReHRV_1.0-main
    conda env create -f  ReHRI_ReHRV_1.0.yml
    conda activate ReHRI_ReHRV_1.0

    
    (2) To speed up the installation process, you can modify the mirror URLs with a trusted source in the "channels" section of the .yml file. 

    (3) Validate install

    python bin/ReHRI.py -v
    ReHRI 1.0

    python bin/ReHRV.py -v
    ReHRV 1.0


**2. Script of searching subconfigurations**

    python bin/ReHRI.py -c ReHRI.config.ini
    
    ReHRI.py: The main program for subconfiguration search.    
    ReHRI.config.ini: The configuration file for subconfiguration search.


**3. Draw the map of recombination organelle genome**

    python bin/ReHRV.py -c ReHRV.config.ini
    
    ReHRV.py: The main program for drawing genome recombination maps.
    ReHRV.config.ini: The configuration file for drawing genome recombination maps.
    

**For more detailed usage instructions, please refer to the tutorial document！！**
  
   
