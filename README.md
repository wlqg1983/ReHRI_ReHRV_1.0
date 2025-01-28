    In organelle genomes, especially mitochondrial genomes, there is a widespread phenomenon of genome recombination mediated by repetitive sequences (direct repeats and inverted repeats), which leads to the emergence of subconfigurations in the genomes. The MiRI software can utilize second- and third-generation sequencing data (reads) to explore the potentially existing subconfigurations within the main configuration of the genome. Meanwhile, it calculates the probabilities of genome recombination mediated by different repetitive sequences based on the number of reads. The MiRIV software can draw genome maps before and after genome recombination to display the genomic structural variations.
    MiRI and MiRIV need to run in an environment with at least Ubuntu (v20.04) and conda (v23.5.2). The tutorials for software installation and usage are as follows:

**1. Install the software**

    Download the software from https://github.com/wlqg1983/MiRI_MiRIV_1.0/archive/refs/heads/main.zip

    unzip MiRI_MiRIV_1.0-main.zip

    cd  MiRI_MiRIV_1.0-main

    conda env create -f  MiRI_MiRIV_1.0.yml

    conda activate MiRI_MiRIV_1.0

    sh Install.sh

    rm Install.sh


**2. Script of searching subconfigurations**

    python bin/MiRI.py -c MiRI.config.ini

    MiRI.py: The main program for subconfiguration search.
    
    MiRI.config.ini: The configuration file for subconfiguration search.


**3. Draw the map of recombination organelle genome**

    python bin/MiRIV.py -c MiRIV.config.ini

    MiRIV.py: The main program for drawing genome recombination maps.

    MiRIV.config.ini: The configuration file for drawing genome recombination maps.


**4. Results**

    The results of genome recombination mediated by different repetitive sequences are stored in a folder named after the project_id.
   
    (1) Results of MiRI
    

    (2) Results of MiRIV

    

**5. Core configuration parameters in the configuration file**
   




