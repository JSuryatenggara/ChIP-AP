<p align="center">
<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/ChIP-AP_logo.png>
</p>
<p align="center">
	<b>Ch</b>romatin <b>I</b>mmuno<b>P</b>recipitation Sequencing <b>A</b>nalysis <b>P</b>ipeline
</p>
<p align="center">
	<b>Inherently Simple ChIP-Seq Analysis</b>
</p>
<p align="center">
	Developed by Jeremiah Suryatenggara and Mahmoud A. Bassal
</p>


## Introduction
ChIP-Seq is a technique used to analyse protein-DNA interactions and identify the binding sites of said protein. A key step in the bioinformatics analysis is the calling of “peaks,” or regions of enrichment, corresponding to binding sites of the protein to the DNA. 


It has been widely reported that of all the bioinformatics analysis steps, the step most likely to affect study outcome is peak calling (Chen et al., 2012; Johnson et al., 2007). Numerous reviews (Koohy et al., 2014; Laajala et al., 2009; Wilbanks and Facciotti, 2010) have been published benchmarking peak callers to ascertain superiority of one caller over others. Consistently however, these benchmarks report mixed, potentially conflicting, findings with no single peak caller out-performing others. Rather, each peak-caller excels in specific types of datasets and one cannot know the overall best performance of a single peak-caller in advance. Using a sub-optimal peak caller will result in few peaks called that do not accurately describe the binding landscape of the investigated protein, or worse, report a large number of false-positive peaks which, short of validating every peak manually, will be indiscernable from true-positive peaks. Compounding issues further, benchmarking reviews have consistently shown little overlap between called peak sets of different peak callers. Therefore, we can conclude that each peak caller has distinct selectivity and specificity characteristics which are often not additive and seldom completely overlap in many scenarios. As such, results obtained from individual peak callers, even sub-optimal callers, will contain sub-sets of true and false-positive peaks that are exclusive to that caller. 


This leads the user to the peak-caller dilemma, if there is no universal peak caller to use, which peak caller _does_ one use for their dataset? We therefore rationalized that without a perfect peak caller (which is near impossible owing to the significant variability introduced from 
wet-lab experiments), the only option is to leverage the best performing, currently available peak callers and identify unique and/or overlapping peaks to better discern the confidence of peaks without manual validation, and identify overlapping, high-confidence peak sub-sets.


To that end, our <b>ChIP-Seq Analysis Pipeline (ChIP-AP)</b> utilizes multiple peak callers, which enables peak calling and merging of results from multiple peak callers (MACS2, GEM/SICER2, HOMER, and Genrich) to: 

1. Selectively obtain high confidence peaks based on overlaps by different peak callers

2. Avoid the caveats of depending on a single peak caller in datasets with different peak characteristics and signal-to-noise ratios

3. Freely choose the “sweet spot” between the consensus peaks set (intersection of detected peaks from all peak callers) and the total peaks set (the union of all detected peaks), that answers one’s biological questions which may be supported by additional evidence from other experiments. This can be done from the output without re-processing data.

ChIP-AP is a fully automated ChIP-seq data processing and analysis pipeline:

1. For input, it takes unaligned sequencing reads in fastq format (extension: .fastq / .fq / .fastq.gz / .fq.gz) or previously aligned reads in bam format (extension: .bam) from any aligner

2. ChIP-AP is capable of processing and analyze multiple sample replicates

3. ChIP-AP is a complete, integrated workflow which performs all analysis steps (QC, cleanup, alignment, peak-calling, pathway analysis) and outputs a final integrated peak set

4. The output of ChIP-AP is a detected peaks list, annotated with information on associated gene names, IDs, functions, gene ontology, and related pathways

<br>

## Computational requirements
<b>OS</b> – Linux (Ubuntu-variants 16.04+ tested), MacOS (10.13+), Windows 10 (v1903+), Windows 11

<b>CPU</b> – (minimum) Quad-Core Intel/AMD CPU, (recommended) Octa-Core Intel/AMD CPU. Apple Silicon Macs are compatible though Rosetta 2.

<b>RAM</b> – (minimum) 8Gb, (recommended) 16Gb+

<b>Storage (SSD/HDD)</b> – This varies widely depending on the samples being processed and the sequencing depth. Rough estimate, (minimum) 30Gb, (recommended) 100Gb+

<b>Screen Resolution</b> – A minimum resolution of 1920*1080 is required for the dashboard interface. If your screen resolution is less than this, you will be limited to only using the wizard
<br>

### Pre-Configured Image
For users who do not wish to try to configure things on their own, we provide a preconfigured virtual machine image for use. Users can download a pre-configured virtual-machine image to run in VirtualBox (Oracle VM VirtualBox). A detailed tutorial can be found on our wiki [here](https://github.com/JSuryatenggara/ChIP-AP/wiki/Virtual-Machine-Installation-Guide). 

We do recommend installing ChIP-AP locally though.
<br>

## Software Installation
ChIP-AP has been designed to be as simple as possible to run for end-users, be they bioinformaticians or wet-lab biologists with no coding experience. For non-experienced users a little command-line dabbling is required in order to get everything set up just right however. We have prepared full installation guides and wiki’s located here on our github for your operating system (OS) of choice (look at the top you will see a “Wiki” tab). These are complete step-by-step guides with accompanying screenshots to walk users through every step of the installation processes. For fresh installations, we have not encountered issues when following our installation instructions accurately. If you have previous instalaltions of configurations, some tweaking may be necessary to get everything to work right – just be aware.

For advanced users, the installation instructions are as follows:

1 - Download the _chipap_installation_ zip in our github repo and unzip in the folder of your choosing.
	
	wget https://github.com/JSuryatenggara/ChIP-AP/raw/main/chipap_installation.zip

2 - Download the pre-generated, required genome files from our [dropbox](https://www.dropbox.com/s/ioqd3hwdahh9xon/genomes.zip) and unzip in the same folder as the ChIP-AP scripts – ie 

	<path>…/chipap_installation/chipap_scripts/
	
3 - Anaconda3 is required for ChIP-AP to install and run. If you have already installed Anaconda3 on your machine, proceed to the next step. If you need to setup Anaconda3, download the required command-line installer from [Anaconda](https://www.anaconda.com/products/individual) using wget. We recommend isntalling anaconda with default settings and when prompted at the end of the installation, initialize the base conda environment.  

Before continuing, you must close and re-open the terminal to initialize conda.

4 - In a new terminal window, you should now see (base) written to the left of the command prompt, this tells you conda has been properly initialized.  If it hasn’t been initialized, refer to the [conda documentation](https://docs.anaconda.com/anaconda/user-guide/faq/#:~:text=In%20order%20to%20initialize%20after,and%20then%20run%20conda%20init%20.&text=Replace%20%3Cpath%2Dto%2Danaconda,of%20your%20installed%20Anaconda%20file.), section _"Should I add Anaconda to the macOS or Linux PATH?"_

5 - Next navigate to the chipap_installation folder. Check that the chipap_installer.py script is executable (if not, run the command)

	chmod +x ./chipap_installer.py

and then execute it. The installer will then ask you a couple of questions regarding installing in its own environment (highly recommended) and then proceed to install everything and all dependencies. Depending on your CPU and internet speeds this will take a while as it needs to download ~10-15Gb of genome data for HOMER.

6 - Once done, to use ChIP-AP, if you installed in its own conda environment you must activate the conda environment before each run. Otherwise, you can simply run 

	chipap_vxx.py
	chipap_wizard.py 
	chipap_dashboard.py

depending on your preference of running.
<br>

## Quick start – For Command line User Only
ChIP-AP is capable of handling multiple sample replicates in a single run. It is also capable of hanlding an un-balanced number of sample replicates (ie 3 ChIP, 2 Controls). It does so by merging each corresponding sample type (so merge all ChIP samples together and merge all control samples together) following read alignment and these merged files are used for all down-stream processing.

For controls, input (not IgG) is recommended and is consensually considered best practice.

For peak calling, peaks are called as ChIP over control.

### Example : To process single-end unaligned reads with default settings:
    chipap  --mode single 
            --chipR1 [chip fastq replicate1] [chip fastq replicate2] … 
            --ctrlR1 [control fastq replicate1] [control fastq replicate2] …
            --genome [path to genome folder]
            --output [path to output save folder]
            --setname [dataset name]

### Example : To process paired-end unaligned reads with default settings:
    chipap  --mode paired 
            --chipR1 [chip fastq replicate1, first read] [chip fastq replicate2, first read] … 
            --chipR2 [chip fastq replicate1, second read] [chip fastq replicate2, second read] … 
            --ctrlR1 [control fastq replicate1, first read] [control fastq replicate2, first read] … 
            --ctrlR2 [control fastq replicate1, second read] [control fastq replicate2, second read] … 
            --genome [path to output save folder] 
            --output [path to output save folder] 
            --setname [dataset name]

### Example : To process single/paired-end aligned reads with default settings:
    chipap	--mode single / paired 
            --chipR1 [chip bam replicate1] [chip bam replicate2] … 
            --ctrlR1 [control bam replicate1] [control bam replicate2] … 
            --genome [path to genome folder] 
            --output [path to output save folder]
            --setname [dataset name]

<br>

### Pipeline Run Command
The command line that was used to call the pipeline will be located in a text file: _[setname]_command_line.txt_ in the output save folder. This is useful for documentation of the run, and for re-running of the pipeline after a run failure or some tweaking if need be.  For example:

    chipap_xx.py --mode paired --ref hg38 --genome [path_to_genome_folder] --output [full_path_to_output_save_folder] --setname [dataset name] --sample_table [path_to_sample_table_file] --custom_setting_table [path_to_setting_table_file].tsv --motif [path_to_known_motif_file] --fcmerge --goann --pathann --deltemp --thread [#_of_threads_to_use] --run

<br>


## Usage notes and Command Line Flags / Parameters
### Required Arguments

| <nobr><b>Argument/Flag</b></nobr> | <nobr><b>Possible Values</b></nobr> | <nobr><b>Detailed Explanation</b></nobr> |
|-|-|-|
| <nobr>--mode</nobr> | <nobr>single / paired</nobr> | Single-end or paired-end sequencing analysis. If a paired-end run, files will typically be labelled ending in *_R1 and *_R2 before the file extension. If these labels aren’t present then likely, you have single-ended sequencing data and select the “single” option. |
| <nobr>--genome</nobr> | <nobr>[directory]</nobr> | Your genome folder directory. Requires full path, not relative path. This is the folder where the pre-computed genome alignment and processing files are saved. These genomes can be downloaded from [dropbox](https://www.dropbox.com/s/ioqd3hwdahh9xon/genomes.zip) or instead you can compute your own (a guide for this is coming soon actually… Keep an eye on the wiki…) |
| <nobr>--output</nobr> | <nobr>[directory]</nobr> | Your desired output folder. Requires full path, not relative path. |
| <nobr>--setname</nobr> | <nobr>[text]</nobr> | The prefix to label output and intermediate files (no space allowed). ChIP-AP will rename all processed data files to have this “setname” prefix. |

### Optional Arguments

| <nobr><b>Argument/Flag</b></nobr> | <nobr><b>Possible Values</b></nobr> | <nobr><b>Detailed Explanation</b></nobr> |
|-|-|-|
| <nobr>--peak</nobr> | <nobr>narrow / broad</nobr> | Narrow peaks for transcription factors (default). Broad peaks for histone modifiers. <br>If unsure what will work best for your dataset, you may need to run ChIP-AP once in each mode and inspect the output to ensure you get what you want from the results. Pay close attention to the width of the peaks called and regions of enrichment. |
| <nobr>--chipR1</nobr> | <nobr>[repl1 repl2 …]</nobr> | Your ChIP datasets: ordered by replicate, separated by space. Best to include full path rather than relative paths. |
| <nobr>--chipR2</nobr> | [repl1 repl2 …] | [ Paired-end Only ] Your ChIP datasets second read: ordered by replicate, separated by space. Best to include full path rather than relative paths. |
| <nobr>--ctrlR1</nobr> | <nobr>[repl1 repl2 …]</nobr> | Your control datasets: ordered by replicate, separated by space. Best to include full path rather than relative paths. |
| <nobr>--ctrlR2</nobr> | <nobr>[repl1 repl2 …]</nobr> | [ Paired-end Only ] Your control datasets second read: ordered by replicate, separated by space. Best to include full path rather than relative paths. |
| <nobr>--sample</nobr> | <nobr>[file]</nobr> | Rather than including the input sample file names and paths on the commandline, one can make a sample table containing the same information and input this instead. The sample-table is a 4*n-sized table (with n = number of replicates) containing the absolute paths to each of your ChIP and control replicates (See below for more information regarding this file and its layout). <br>When this option is used, this table will disregard any assigned values to --chipR1, --chipR2, --ctrlR1, and --ctrlR2. |
| <nobr>--setting</nobr> | <nobr>[file]</nobr> | [ For Advanced Users ONLY ] The settings-table allows you to fine-tune the behaviour of every program as part of ChIP-AP. Hence very disasterous if you get wrong! If you are unsure even a little about what you’re doing, then stick to default settings please – this goes even for bioinformaticians. <br>This txt file is a 2*n-sized table (with n = number of replicates) containing custom arguments for every program as part of ChIP-AP. The default settings table is provided in the genome folder. You can COPY this file and make changes as necessary. To use your customs settings table, provide full path to updated txt file. |
| <nobr>--ref</nobr> | hg19 / hg38 /<br> mm9 / mm10 /<br> dm6 / sacCer3 | Your sample genome reference build. Default is hg38 (human). The genomes listed to the left are provided, pre-calculated by us and are the only genomes used and tested for now. We will provide added functionality soon to add your own custom genome to ChIP-AP. <br>Watch this space (actually the github…)! |
| <nobr>--motif</nobr> | <nobr>[file]</nobr> | Your predicted/known motif file, in HOMER matrix format (.motif). If provided, once peaks are called, HOMER motif discovery will be run on the total called peak set for this motif. Many users prefer MEME-ChIP instead… we know… WIP :) |
| <nobr>--fcmerge</nobr> | | This flag will force fold change analysis to be computed based on merged replicates instead of on each replicate seperately.
| <nobr>--goann</nobr> | | This flag will instruct to annotate peaks with all relevant GO terms as provided by HOMER. |
| <nobr>--pathann</nobr> | | This flag will instruct to annotate peaks with all relevant pathway and interaction enrichment terms as provided by HOMER. |
| <nobr>--deltemp</nobr> | | This flag will instruct to delete large intermediary files right after they are not going to be used for further processes (eg intermediary fq files. This option will save a significant amount of space by the end of your analysis, so recommended. |
| <nobr>--thread</nobr> | <nobr>[integer]</nobr> | Maximum number of processes to use. Default is half the maximum available on your system so as to not choke it during the run. If running on a laptop or low-thread count cpu, best to push this up to maximum number threads available -1 – but this will significantly slow your laptop if attempting other tasks while ChIP-AP is running. |
| <nobr>--run</nobr> | | Use to immediately run the suite by running the master script. This is the big red button ok... Use at your own risk! <br> When not used, the generated master script (MASTER_script.sh) in the output folder can be run manually by user. We have made sure that when you tell ChIP-AP “do not press the big red button” that it will behave and do as you say. |
| <nobr>--homer_motif</nobr> | consensus /<br> union / both | This flag will instruct HOMER (findMotifsGenome.pl) to perform motif enrichment analysis at the end of the pipeline run. |
| <nobr>--meme_motif</nobr> | consensus /<br> union / both | This flag will instruct MEME (meme-chip) to perform motif enrichment analysis at the end of the pipeline run. |

<br>


## Running ChIP-AP with Graphical Interfaces
ChIP-AP offers 2 graphical interfaces for running – depending on a users proficiency with ChIP-AP. They are the _chipap_dashboard_ and the _chipap_wizard_.

We recommend new users to ChIP-AP first use the wizard. Proficient users of ChIP-AP should use the dashboard as it enables inputting the data more quickly.

### Running the Wizard
Once ChIP-AP is setup, at the command line type: 

	chipap_wizard.py

to use the wizard graphical user interface (GUI). The multiple windows that appear ask users questions sequentially regarding information required for a run. This ensure data is provided in the correct order without overwhelming the user with all questions and options simultaneously. This option is therefore recommended for inexperienced ChIP-AP users. Tooltip boxes appear over entry fields providing brief explanations of input requirements.
<p align="center">
<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/chipap_wizard_guide.png>
</p>

### Running the Dashboard
Once ChIP-AP is setup, at the command line type: 

	chipap_dashboard.py

The GUI below will appear enabling users to input everything required for a run. Tooltip boxes will appear over each field providing brief explanations of whats required in said field.
<p align="center">
<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/chipap_dashboard_guide.png>
</p>

### For advanced users
A customized pipeline run is possible by providing the required flag arguments for each individual program in the pipeline at the * stage above. By selecting the right option, a secondary window will appear containing multiple entry fields for such flag arguments. This option is for advanced users only. For valid flags, please refer to the documentation of each program specifically as each program has different flags and syntax for its parameters.
<p align="center">
<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/setting_table_guide.png width="500">
</p>
In each of these entry fields, user may key in the flag arguments as they are written for each corresponding program (e.g., “-q 20” for samtools view which will filter reads with a MAPQ > 20). This requires user to read dedicated official manuals to understand what and how these arguments can be given to the corresponding program in the pipeline. This, accompanied with the fact that bad/wrong flag arguments may break the pipeline run, shows that this manual customization of pipeline programs settings is recommended for experienced users only, and that inexperienced users should stick with the default settings or get a bioinformaticians input.
<br>



## ChIP-AP Graphical Overview
<p align="center">
<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/chipap_workflow_guide.png>
</p>
<b>For detailed explanations of all the steps and methodologies used throughout ChIP-AP refer to our documentation (https://github.com/JSuryatenggara/ChIP-AP/wiki/ChIP-AP-Guide)</b>
<br>


## Main Pipeline Output
### Final Analysis Table (including supplementary annotations)
The table below shows the contents of [setname]_all_peaks_go_pathway_annotated.tsv. Smaller sized and less verbose variants of this table are saved in the output folder with suffixes: concatenated, annotated, and calculated (see Source in the table below)

<table>
    <thead>
        <tr>
            <th>Column #</th>
            <th>Peak Attribute</th>
            <th>Source</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td><nobr><b>Column 1 (A)</td>
            <td><nobr>Peak ID (unique peak ID)</td>
            <td rowspan=6>
                <b>Pipeline script:</b><br>22_peaks_ processing _script.sh
                <br><b>Called program:</b><br>cat (Bash)
                <br><b>Output file:</b><br>[setname]_all_peaks_concatenated.tsv
                <br><b>Output folder:</b><br>22_peaks_processing
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 2 (B)</td>
            <td><nobr>Chr (chromosome)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 3 (C)</td>
            <td><nobr>Start (peak start coordinate)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 4 (D)</td>
            <td><nobr>End (peak end coordinate)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 5 (E)</td>
            <td><nobr>Strand (on which peak is found)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 6 (F)</td>
            <td><nobr>Peak Caller Combination</td>
        </tr>
        <tr>
            <td><nobr><b>Column 7 (G)</td>
            <td><nobr>Peak Caller Overlaps</td>
            <td rowspan=5>
                <b>Pipeline script:</b><br>22_peaks_processing_script.sh
                <br><b>Called script:</b><br>fold_change_calculator_suite_xx.py
                <br><b>Output file:</b><br>[setname]_all_peaks_calculated.tsv
                <br><b>Output folder:</b><br>22_peaks_processing
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 8 (H)</td>
            <td><nobr>ChIP Tag Count</td>
        </tr>
        <tr>
            <td><nobr><b>Column 9 (I)</td>
            <td><nobr>Control Tag Count</td>
        </tr>
        <tr>
            <td><nobr><b>Column 10 (J)</td>
            <td><nobr>Fold Change</td>
        </tr>
        <tr>
            <td><nobr><b>Column 11 (K)</td>
            <td><nobr>Number of Motifs</td>
        </tr>
        <tr>
            <td><nobr><b>Column 12 (L)</td>
            <td><nobr>Annotation</td>
            <td rowspan=14>
                <b>Pipeline script:</b><br>22_peaks_processing_script.sh
                <br><b>Called program:</b><br>HOMER annotatePeaks
                <br><b>Output file:</b><br>[setname]_all_peaks_annotated.tsv
                <br><b>Output folder:</b><br>22_peaks_processing
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 13 (M)</td>
            <td><nobr>Detailed Annotation</td>
        </tr>
        <tr>
            <td><nobr><b>Column 14 (N)</td>
            <td><nobr>Distance to TSS</td>
        </tr>
        <tr>
            <td><nobr><b>Column 15 (O)</td>
            <td><nobr>Nearest PromoterID</td>
        </tr>
        <tr>
            <td><nobr><b>Column 16 (P)</td>
            <td><nobr>Entrez ID</td>
        </tr>
        <tr>
            <td><nobr><b>Column 17 (Q)</td>
            <td><nobr>Nearest Unigene</td>
        </tr>
        <tr>
            <td><nobr><b>Column 18 (R)</td>
            <td><nobr>Nearest Refseq</td>
        </tr>
        <tr>
            <td><nobr><b>Column 19 (S)</td>
            <td><nobr>Nearest Ensembl</td>
        </tr>
        <tr>
            <td><nobr><b>Column 20 (T)</td>
            <td><nobr>Gene Name</td>
        </tr>
        <tr>
            <td><nobr><b>Column 21 (U)</td>
            <td><nobr>Gene Alias</td>
        </tr>
        <tr>
            <td><nobr><b>Column 22 (V)</td>
            <td><nobr>Gene Description</td>
        </tr>
        <tr>
            <td><nobr><b>Column 23 (W)</td>
            <td><nobr>Gene Type</td>
        </tr>
        <tr>
            <td><nobr><b>Column 24 (X)</td>
            <td><nobr>CpG%</td>
        </tr>
        <tr>
            <td><nobr><b>Column 25 (Y)</td>
            <td><nobr>GC%</td>
        </tr>
        <tr>
            <td><nobr><b>Column 26 (Z)</td>
            <td><nobr>Biological Process</td>
            <td rowspan=3>
                <b>Pipeline script:</b><br>23_go_annotation_script.sh
                <br><b>Called script:</b><br>GO_annotate_suite_xx.py
                <br><b>Output file:</b><br>[setname]_all_peaks_go_annotated.tsv
                <br><b>Output folder:</b><br>23_supplementary_annotations
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 27 (AA)</td>
            <td><nobr>Molecular Function</td>
        </tr>
        <tr>
            <td><nobr><b>Column 28 (AB)</td>
            <td><nobr>Cellular Component</td>
        </tr>
        <tr>
            <td><nobr><b>Column 29 (AC)</td>
            <td><nobr>Interaction with Common Protein</td>
            <td rowspan=8>
                <b>Pipeline script:</b><br>23_pathway_annotation_script.sh
                <br><b>Called script:</b><br>pathway_annotate_suite_xx.py
                <br><b>Output file:</b><br>[setname]_all_peaks_pathway_annotated.tsv
                <br><b>Output folder:</b><br>23_supplementary_annotations
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 30 (AD)</td>
            <td><nobr>Somatic Mutations (COSMIC)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 31 (AE)</td>
            <td><nobr>Pathway (KEGG)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 32 (AF)</td>
            <td><nobr>Pathway (BIOCYC)</td>
        </tr>
        <tr>
           <td><nobr><b>Column 33 (AG)</td>
            <td><nobr>Pathway (pathwayInteractionDB)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 34 (AH)</td>
            <td><nobr>Pathway (REACTOME)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 35 (AI)</td>
            <td><nobr>Pathway (SMPDB)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 36 (AJ)</td>
            <td><nobr>Pathway (Wikipathways)</td>
        </tr>
    </tbody>
</table>

<br>


### Pipeline Run Info
This file summarizes the assignment of the files (IP sample or control, read 1 or 2; replicate number) and the file name conversion for every unaligned or aligned sequencing reads to be processed. Each line tells the user what the original files have been renamed into. Check this file if you suspect the order of samples were incorrectly entered (i.e., swapped chip with control)

    Chromatin IP dataset replicate 1, 1st read : Original filename = a.fastq --> New filename = setname_chip_rep1_R1.fq.gz

    Chromatin IP dataset replicate 2, 1st read : Original filename = b.fastq --> New filename = setname _chip_rep2_R1.fq.gz  

    Chromatin IP dataset replicate 1, 2nd read : Original filename = c.fastq --> New filename = setname _chip_rep1_R2.fq.gz  

    Chromatin IP dataset replicate 2, 2nd read : Original filename = d.fastq --> New filename = setname _chip_rep2_R2.fq.gz  

    Control dataset replicate 1, 1st read : Original filename = e.fastq --> New filename = setname _ctrl_rep1_R1.fq.gz  

    Control dataset replicate 2, 1st read : Original filename = f.fastq --> New filename = setname _ctrl_rep2_R1.fq.gz  

    Control dataset replicate 1, 2nd read : Original filename = g.fastq --> New filename = setname _ctrl_rep1_R2.fq.gz  

    Control dataset replicate 2, 2nd read : Original filename = h.fastq --> New filename = setname _ctrl_rep2_R2.fq.gz

<br>


### Sample Table
Contains the full path of each input ChIP and control sample in the pipeline run in a tab-separated value file: [setname]_sample_table.tsv in the output save folder in ChIP-AP sample table format. This is useful for documentation of the run, and for re-running of the pipeline after a run failure or some tweaking if need be. Below is an example of sample table file content (header included), given paired-end samples with two ChIP replicates and two control replicates:

| chip_read_1 | chip_read_2 | ctrl_read_1 | ctrl_read_2 |
|-|-|-|-| 
| ... /a.fastq | ... /c.fastq | ... /e.fastq | ... /g.fastq |
| ... /b.fastq | ... /d.fastq | ... /f.fastq | ... /h.fastq |

If your sample is single-ended, then the '..._read_2' columns can be simple be left blank, or the sample table can be simplified to as follows:

| chip_read_1 | ctrl_read_1 |
|-|-|
| ... /a.fastq | ... /e.fastq |
| ... /b.fastq | ... /f.fastq |

<br>

### Setting Table & Default Parameters
A cornerstone of ChIP-AP’s functionality is the settings table. ChIP-AP, with the raw fq files and the settings table, is able to reproduce (near) identically any analysis that was performed (provided the same program version numbers are used). The ‘near identically’ statements is owing to the fact that reported alignments of multi-mappers may, in some cases, give every so slightly different results. This ambiguity can be alleviated however by filtering out alignments with low MAPQ scores in the corresponding alignment filter step, post-alignment to ensure consistent results from every analysis run. The provision of the settings table therefore ensures reproducibility of any analysis with minimal effort and bypasses the usually sparse and significantly under-detailed methods sections of publications. Science is supposed to be reproducible, yet bioinformatics analysis are typically black-boxes which are irreproducible. This 1 file, changes that!  

The structure of the settings table is simple. It is a 2 column tab-separated value file with the names of the programs on the 1st column, and the necessary flags required or changed in the 2nd column. If making your own custom table, then the 1st column below must be copied as-is and not changed. These 2 columns together, list the flags and argument values for each program used in the pipeline.  

When ChIP-AP is run, a copy of the used settings table is saved as a tab-separated value file: [setname]_ setting_table.tsv in the output save folder. 

If you have a custom settings table made and provided it as input, then ChIP-AP will make a 2nd copy of this table in the same output save folder. This decision is made as it is useful documentation of the run performed. This file is also useful for re-running of the pipeline after run failure or some tweaking if necessary. If submitting an issue request on Github, you must provide us your settings table used as well as all other requested information. See Github for details regarding this.

**We consider the dissemination of the information of this file as vital and essential along with results obtained. The table can be included as a supplemental table in a manuscript or can be included as a processed data file when submitting data to GEO – either way, the information of this file must be presented when publishing data.**
		    
Below is an example of setting table file in its default-setting state:

**Note:** When running GEM, to ensure absolute reproducibility of all results every single run, you **MUST** add the flag "--t 1" to the DST below. When GEM is run in multi-threaded mode, its results are not reproducible over multiple runs even with the same parameters (although they are very close to each other). This is an issue we raised with the developers and is out of our hands. We highly doubt they will fix this aspect of GEM's behaviour. Therefore if 100% reproducibility of results is key for you, run in single-threaded mode only - even if it will take longer to run.
		    

| program | argument |
|-|-|
| <b>fastqc1 | -q |
| <b>clumpify | dedupe spany addcount qout=33 fixjunk |
| <b>bbduk | ktrim=l hdist=2 |
| <b>trimmomatic | LEADING:20 SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:20 |
| <b>fastqc2 | -q |
| <b>bwa_mem |  |
| <b>samtools_view | -q 20 |
| <b>plotfingerprint |  |
| <b>fastqc3 | -q  |
| <b>macs2_callpeak |  |
| <b>gem | -Xmx10G --k_min 8 --k_max 12 |
| <b>sicer2 |  |
| <b>homer_findPeaks |  |
| <b>genrich | --adjustp -v |
| <b>homer_mergePeaks |  |
| <b>homer_annotatePeaks |  |
| <b>fold_change_calculator | --normfactor uniquely_mapped |
| <b>homer_findMotifsGenome | -size given -mask | 
| <b>meme_chip | -meme-nmotifs 25 |

As can be seen, certain flags and values for some programs have been preset as per our testing and opinions. A point to note however, some flags for programs, such as -BAMPE in MACS2, are not listed since they are “hard-coded” into the pipeline and cannot be modified. For this example of -BAMPE in MACS2, this is “hard-coded” because this flag is essential for running peak calling in paired-end datsets. Parameters and flags like this that must be set are “hard-coded” and hidden and cannot be changed unless by choosing the appropriate narrow/broad run modes. A listing of all these “hard-coded” parameters can be found in our [documentation](https://github.com/JSuryatenggara/ChIP-AP/wiki/ChIP-AP-Guide).
<br>


## ChIP vs Control fold change calculation
For transcription factor samples (narrow peaks), ChIP weighted peak center coordinate is determined by the median coordinate of all reads in the peak region. Fold change was then calculated as read depth in ChIP sample divided by non-zero read depth of control sample at the weighted peak center coordinate.

For histone modifier samples (broad peak type), fold change was simply calculated based on average read depth in ChIP sample, divided by non-zero average read depth of control sample, along the same peak region.

To correct for read depth discrepancies due to imbalanced number of reads between ChIP and control samples, ChIP-AP uses one of three available normalization factors depending on the commandline flag specified:
- (Default) Based on only uniquely mapped reads in ChIP vs control bam files (filtered using samtools view -F256 and counted using samtools view -c)
- Based on only successfully mapped reads in ChIP vs control bam files (filtered using samtools view -F4 and counted using samtools view -c) 
- Based on user-determined value.

Advanced users may choose to change this by changing the argument field for fold_change_calculator in the settings table as follows:
- --normfactor uniquely_mapped 
    - _Default_ - based on uniquely mapped reads; default setting
- --normfactor mapped 
    - change to all mapped reads
- --normfactor user_value --chip_norm [x] --ctrl_norm [y] 
    - change to user-determined normalization factor, where x / y is the user-determined ratio of the number of ChIP reads to the number of control reads
<br>

## Genrich p-value threshold adjustment formula
Based on our testing, Genrich tends to misbehave when processing datasets with low read depth. It starts to call low enrichment regions as peaks, leaving users with an impossibly large number of called peaks. We contacted the developer of Genrich and confirmed such behaviour with them. Through our testing, we derived equations that allow us to curtail Genrich’s aberrant behaviour in such scenarios. This adjustment is performed and derived by us and is not attributed to Genrich and/or its developer(s).

To avoid the aforementioned peak calling depth issue, we curtail such behaviour by setting up an auto-adjusting p value threshold that responsively raises the default limit of peak’s minimum p value as the read depth gets lower using the following equations Q1, Q2, Q3.
<p align="center">
<img src="https://github.com/JSuryatenggara/ChIP-AP/blob/storage/images/genrich_equation_guide.PNG" width="500">
</p>
<br>

## Interpreting ChIP-AP Output
Ok so ChIP-AP does report a fair amount of stuff. If you ran it locally you have a swath of folders and you have nooooo clue what to look for and its all confusing. We get that. The reality though its very simple to know what to look for to know your experimental run worked and in this section were going to walk you through that!

### Did my analysis work?
There are a couple of things to look for to answer this question. 1, the fingerprint plot and 2, the venn diagram of the merged peaks. Let's begin…

1. <b>The Fingerprint Plot</b>

    The fingerprint plot tells us how well the enrichment of your samples worked. It is generated by the function from the deeptools package and is generated after the alignment files. As such, the plots are found in the “08_results” folder and are labelled “fingerprint_xxxxx.png/svg.” The PNG files allow you to view them in any image viewer, the SVG files are for opening in Adobe Illustrator or Inkscape to make HQ publication figures later if you need.
<p align="center">
    <img src="https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/fingerprint_plot_guide.png" width="500"> 
</p>
    To interpret the fingerprint plot, (more information can be found on the deeptools documentation site), but  put simply, the input control should be a diagonal line as close as possible toward the 1:1 diagonal. Your ChIP sample should have a bend/kink towards the bottom right corner. The greater the separation between the input and the chip sample, the greater the enrichment you will see in the final result (i.e., lots of peaks). If the lines are overlapping, then you will see little enrichment and your experiment didn’t work that well. If you’re sample lines are switched – then you probably switched the sample names and we recommend doing the right thing and repeating the experiment and not simply switch sample names for the sake of a publication.

   In this example, there is reasonable enrichment in our chip samples. And so we are confident we can see enrichment.


2. <b>The Venn Diagram (well Venn Text)</b>

    In the folder “21_peaks_merging” folder, you will find the “venn.txt” file. This will show you a textual venn diagram of the overlap between the called peaks across all peak callers. To know your experiment worked well, then you should see a full list with combinations of all peak callers and relatively large numbers for the consensus peak sets (ie peaks called by multiple peak callers) – this is the ideal case. However, from our experiences, there will almost always be 1 maybe 2 peak callers that don’t like a dataset for some reason and so you may find a peak caller performed poorly but the others performed admirably. This is still a good and valid result. 
<p align="center">
    <img src="https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/venn_diagram_guide.png" width="500">
</p>
    If you look at this file and only see a small number of peaks and little overlap, and only 1 peak caller seems to have dominated peak calling, then likely your experiment didn’t work that great. Just because only 1 peak caller performed well though, doesn’t mean the experiment is a write-off and a failure. It can still be valid and so doings some manual validations on the top FC differential peaks by chip-PCR might give you an indication whether there is salvageable data or not. Also if you have other confirmatory experimental evidence then even 1 peak calling getting results is fine. This is why we implemented multiple peak callers, because there are many instances where the signal:noise just creates a mess for most peak callers but generally 1 will be the super-hero of the day in such a situation.

3. <b>What Results Files Do I Look At Exactly?</b>

    Valid question. In the folder “22_peak_processing,” open the “xxxx_all_peaks_calculated.tsv” file in excel and you’re good to go. Now to open it there is a little step to do…

    Open a new blank workbook in excel
<p align="center">
    <img src="https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/excel_screenshot_1_guide.png" width="500">
</p>
    In the ribbon at the top, go to “Data”, then select “From Text/CSV”
<p align="center">
    <img src="https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/excel_screenshot_2_guide.png" width="500">
</p>
    In the dialog box that opens up, find and open the peaks files “xxxx_all_peaks_calculated.tsv.” Follow all the prompts and keep pressing “Next” / “Proceed” till the end and the file opens. Opening the peak file this way circumvents an issue that Excel constantly makes which is it will interpret some gene names such as OCT1 as a date, when its not. So by following the aforementioned steps, excel will not do this stupid conversion and instead, when you save the file as an xlsx, it will ensure that this issue doesn’t happen (seen it in sooooo many publications its not funny – just import data this way please people?)

   From this file, you can view all the results and data for you analysis. Refer to Interpreting ChIP-AP Output for the definition of what each column means.

4. <b>How Do I View My Alignments And Data?</b>
    
    People typically want to view their results on UCSC or other genome browsers. As we don’t have a web-server to host such coverage files (and making an accessible ucsc hub is a real pain and we don’t want to implement that), the onus is on you to view them locally on your machine. All laptops, whether then can run ChIP-AP or not can run [IGV](https://software.broadinstitute.org/software/igv/download) and view the coverage and bam files. The coverage and bam failes can be located in the “08_results” folder. 

    Download [IGV](https://software.broadinstitute.org/software/igv/download), install it (super easy) and then load the coverage and bam files needed. Make sure you load the right genome build however! That’s critical. From the ChIP-AP main output file (see the _Main Pipeline Output_ section), you can copy columns B,C,D straight into IGV and it will take you to the peak region. 

5. <b>In Short, What's Relevant?</b>  
    Easy answers
    1. Check fingerprint plot and make sure it looks good
    2. Check venn.txt file and make sure you get good spread of peaks
    Together points 1 and 2 tell you your experiment worked!
    3. Your final peak file is in “22_peak_processing” open the “xxxx_all_peaks_calculated.tsv” – This is the file you need to upload to GEO as your processed data file for your analysis and the only file you need to work with when looking through your data.
    4. Also as part of your submission to GEO or as a supplemental table in your manuscript, you MUST include the settings table named “default_settings_table.txt” located in the root analysis directory. This provided with the raw fq files, which must be uploaded to GEO, will ensure complete reproducibility of the analysis performed.
    5. Manuscript details for M&M. A statement such as the following should suffice:
    
        For processing our ChIP-Seq analysis, we utilized ChIP-AP (https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbab537/6489109). Raw fq files are uploaded to GEO with accession number GSE172355, and the custom settings table utilized for analysis can be found on GEO as a processed settings file and also in Table 1 in our manuscript. Full details of ChIP-AP and its function can be found in its corresponding manuscript (https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbab537/6489109).

<br>
		    
<b>We have a Q&A section located as a wiki page with frequently asked questions and answers. This list will be updated as we receive more questions to advise the wider community.</b>
