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

    The results of genome recombination mediated by different repeat sequences are stored in a folder named after the project_id.
    
(1) Core results of MiRI (final_repeat-spanning_results_AR/paired_repeats_recomb-supporting_ratio.tsv)
![10](https://github.com/user-attachments/assets/3e620ae4-5afd-47bc-91b2-6398874ddc0f)

① fragment_id: The ID of repeat sequence. 
② length: The length of repeat sequence.
③ start: The start positin of repeat sequence in the genome.
④ end: The end positin of repeat sequence in the genome.
⑤ direction: The postion in plus or minus strand of DNA.
⑥ chromosome: The name of chromosome in genome.
⑦ plus_ratio(s/m): The probability of repeat sequence-mediated genome recombination on the plus strand of DNA.
⑧ minus_ratio(s/m): The probability of repeat sequence-mediated genome recombination on the minus strand of DNA.
⑨ combined_ratio: The overall ratio of repeat sequence-mediated genome recombination on the two strands of DNA.
⑩ type: The type of repeat (direct or inverse repeat).
○11 paired: Represents the other repeat sequence in a pair of repeat sequences that mediate genome recombination.


    (2) Results of MiRIV
![13](https://github.com/user-attachments/assets/06688f10-42a7-4e49-8e91-6d97ed34acce)
![14](https://github.com/user-attachments/assets/7c131b36-61fd-4fbc-a7a7-4100dd7dcc81)
![15](https://github.com/user-attachments/assets/48fbbf6d-1c14-491e-bdc6-d617fd68ac81)
![nine_squares_VV](https://github.com/user-attachments/assets/b12e443b-68a5-4512-b64b-a00874525c67)


**5. Core configuration parameters in the configuration file**
  
    (1) Core configuration file of MiRI
![1 2 3](https://github.com/user-attachments/assets/1b8531bb-2afd-4f75-ae68-b9abf7bbb8d2)
![4](https://github.com/user-attachments/assets/cbee84a6-757f-40c1-8be3-319695dff202)
![5 6 7](https://github.com/user-attachments/assets/1c3031b7-8d18-486c-8c2d-259a46075ae9)

① Set project ID.  
② Organelle genome sequence.  
③ Circular or linear genome.  
④ Set the length of repeat sequence.  
⑤ Set the NGS reads (single or paired).  
⑥ Set thr TGS reads.  
⑦ The type of sequencing platform (pacbio or ont) of TGS reads.


    (2) Core configuration file of MiRIV
    
① Set the project ID.

![1](https://github.com/user-attachments/assets/458090df-cd1f-49ea-8925-b674e7924801)

② If draw the map of mainconfiguration.

![2](https://github.com/user-attachments/assets/31c126fc-ad1a-4673-944d-688324011518)

③ If draw the map of subconfiguration mediated by IR.

![3](https://github.com/user-attachments/assets/95d695cd-406f-4022-9cc4-e6d9901d573d)

④ If draw the map of subconfiguration mediated by DR (1to2).

![4-4](https://github.com/user-attachments/assets/953aac22-3942-4667-ae16-fc6e760c203f)

⑤ If draw the map of subconfiguration mediated by DR (2to1).

![5](https://github.com/user-attachments/assets/fd8888c0-0cbf-461c-bbc7-2d8e40ca2f94)

⑥ If draw the map of subconfiguration mediated by DR (2to2).

![6](https://github.com/user-attachments/assets/4893bb2f-a5d1-4273-835c-484f9904c6fc)

⑦ If arrange the drawed maps in a grid of nine squares.

![7](https://github.com/user-attachments/assets/a4d6e947-4562-431c-a993-13268a3d1b97)


