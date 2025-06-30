    In organelle genomes, especially mitochondrial genomes, there is a widespread phenomenon of genome recombination mediated by repetitive sequences (direct repeats and inverted repeats), which leads to the emergence of subconfigurations in the genomes. The ReHRI software can utilize second- and third-generation sequencing data (reads) to explore the potentially existing subconfigurations within the main configuration of the genome. Meanwhile, it calculates the probabilities of genome recombination mediated by different repetitive sequences based on the number of reads. The ReHRV software can draw genome maps before and after genome recombination to display the genomic structural variations.
    ReHRI and ReHRV need to run in an environment with at least Ubuntu (v20.04) and conda (v23.5.2). The tutorials for software installation and usage are as follows:

**1. Install the software**

    (1) Download and install the software from https://github.com/wlqg1983/ReHRI_ReHRV_1.0/archive/refs/heads/main.zip

    unzip ReHRI_ReHRV_1.0-main.zip
    cd  ReHRI_ReHRV_1.0-main
    conda env create -f  ReHRI_ReHRV_1.0.yml
    conda activate ReHRI_ReHRV_1.0
    sh Install.sh
    rm Install.sh
    
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
    

**4. Core configuration parameters in the configuration file**
  
    (1) Core configuration file of ReHRI
![1 2 3](https://github.com/user-attachments/assets/1b8531bb-2afd-4f75-ae68-b9abf7bbb8d2)
![4](https://github.com/user-attachments/assets/cbee84a6-757f-40c1-8be3-319695dff202)
![5 6 7](https://github.com/user-attachments/assets/1c3031b7-8d18-486c-8c2d-259a46075ae9)

① Set project ID.  
② Organelle genome sequence.  
③ Circular or linear when one chromosome is present.  
④ Set the length of repeat sequence.  
⑤ Set the NGS reads (single or paired).  
⑥ Set the TGS reads.  
⑦ The type of sequencing platform (pacbio or ont) of TGS reads.


    (2) Core configuration file of ReHRV
    
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


**5. The core input tsv file of ReHRV**

Provide the information of each repeat pair (.tsv).

![16](https://github.com/user-attachments/assets/f0cdbf80-9173-4a3c-b016-86f9e9981574)


**6. Results**

    The results of genome recombination mediated by different repeat units are stored in a folder named after the project_id.
    
(1) Core results of ReHRI (final_repeat-spanning_results_AR/paired_repeats_recomb-supporting_ratio.tsv)
![10](https://github.com/user-attachments/assets/3e620ae4-5afd-47bc-91b2-6398874ddc0f)

① fragment_id: The ID of repeat unit.

② length: The length of repeat unit.

③ start: The start positin of repeat unit in the genome.

④ end: The end positin of repeat unit in the genome.

⑤ direction: The postion in plus or minus strand of DNA.

⑥ chromosome: The name of chromosome in genome.

⑦ plus_ratio(s/m): The probability of repeat-mediated genome recombination on the plus strand of DNA.

⑧ minus_ratio(s/m): The probability of repeat-mediated genome recombination on the minus strand of DNA.

⑨ combined_ratio: The overall ratio of repeat-mediated genome recombination on the two strands of DNA.

⑩ type: The type of repeat (direct or inverse repeat).

○11 paired: Represents the other repeat unit in a pair of repeat units that mediate genome recombination.

    (2) Results of ReHRV

![13](https://github.com/user-attachments/assets/06688f10-42a7-4e49-8e91-6d97ed34acce)
![14](https://github.com/user-attachments/assets/7c131b36-61fd-4fbc-a7a7-4100dd7dcc81)

The maps of Direct repeat-mediated recombination genome.

Linear genome.
![15](https://github.com/user-attachments/assets/48fbbf6d-1c14-491e-bdc6-d617fd68ac81)

The grid of nine squares.
![nine_squares_VV](https://github.com/user-attachments/assets/b12e443b-68a5-4512-b64b-a00874525c67)


**7. Advanced features of ReHRI**

(1) Detect custom repeats for genome recombination
![Snipaste_2025-02-14_23-37-59](https://github.com/user-attachments/assets/30bea6e6-e8c8-4692-83f5-1b05909f2f58)

① Set mode=C to enter the mode where the user provides repeated sequences.

② The complete genome sequence is also required.

③ Provide user-defined repeats that may mediate genome recombination.

④ The data format (.tsv, 5 or 6 colume tables) of user-defined repeat information.

![Snipaste_2025-02-15_00-10-09](https://github.com/user-attachments/assets/54e3e909-62a2-4598-ac82-730597730f8a)

The prediction of repeats is a challenging task, and different algorithms may obtain different repeats. This model enables users to study the impact of specific repeats on genome recombination.

(2) Refilter the repeat that can mediate the recombination of organelle genome
![20250214232935](https://github.com/user-attachments/assets/b52d3cae-b4d0-4a14-8702-ad8de85a77f0)

① Set refiter_made=Y to enter the mode of re-filtering results.

② Set a project ID.

③ Reset value of spanning_read_flanking_repeat_length.

④ Reset value of spanning_read_number.

The refiltering operation can filter multiple times based on previous results.

