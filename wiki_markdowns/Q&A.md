<b>Below is a listing of frequently asked questions and their answers.  We will keep updating this section as we get feedback and common questions from the community.</b>

## Q&A
1. <b>Q – What does ChIP-AP stand for?</b>  
	<b>A –</b> ChIP-Seq Analysis Pipeline

2. <b>Q – What should I use as my control? Input or IgG?</b>  
	<b>A –</b> The consensually accepted control to be used for ChIP-Seq experiments is input control and not IgG. You’re free to do so if you wish but that’s not whats commonly done in the literature. Take a look at consortia such as ENCODE and you won’t find IgG controls, only input.

3. <b>Q – How are peaks called? In what order?</b>  
	<b>A –</b> Peaks are called as your ChIP-IP over Input control.

4. <b>Q – Can ChIP-AP handle an un-balanced number of replicates? (ie 3 ChIP and 1 input?)</b>  
	<b>A –</b> Yes it can. The initial QC and cleanup steps are done on a per replicate basis, but after the alignment step, samples are merged and the remaining peak-calling, annotation, FC calculations etc… are done on the merged sample sets.
5. <b>Q – I’ve only input 1 sample ChIP and 1 Input control, why are there “_merged” files in the output directories?</b>  
	<b>A –</b> Yeah… ideally ChIP-AP shouldn’t do this and we’re aware of this behaviour. The contents of the single replicate file and the merged are the same though so theres not really an issue which file you refer to as long as you’re consistent. We will address this small issue in future updates.

6. <b>Q – Why are you using BWA instead of other aligners such as Bowtie2?</b>  
	<b>A –</b> BWA is used in preference to other aligners such as Bowtie2 as in benchmarking papers such as (Thankaswamy-Kosalai et al., 2017) and others, we were more statisfied with the results and performance of BWA over others, hence its inclusion in this pipeline. Personal preference ppl.

7. <b>Q – Can ChIP-AP handle bam files from other aligners as input, or only BWA?</b>  
	<b>A –</b> Technically it should be able to handle alignment bam files from any aligner provided the files are formatted correctly. We have not tested it thoroughly though, but should still be fine. If there are issues let us know through the githib and we can work together to resolve such matters.
8. <b>Q – I have developed this incredibly awesome new peak caller called Peak2DaMax, can it be incorporated into ChIP-AP?</b>  
	<b>A –</b> Owing to the modular nature of ChIP-AP, it is realtively easy for additional programs to be added/removed if need be. You can attempt to modify the source files to get this working and reach out if you have any issues. As for officially including it into ChIP-AP… Reach out to us and we can zoom.
9. <b>Q – Can ChIP-AP be setup on a shared computing cluster?</b>  
	<b>A –</b> Technically yes it can. Have we done it? Yes (and no, we have access to 2 clusters, we’ve set it up on 1 and not the other yet). Was the installation more involved? Yes. Every shared computing cluster has its own configuration in terms of available software and personal profile configuration requirements. Its too difficult for us to do every possible combination in advance. Refer to the References and Citations section to see what programs need to be installed/available to for ChIP-AP to run and setup your environment accordingly to access those programs. ChIP-AP will also run a pre-flight check before each run to ensure everything is available and accessible. So check the logs to make sure ChIP-AP is happy and then you will be ok. If unfamiliar with working in said environments, contact the required IT support for the shared cluster to help you.
10. <b>Q – The installer script downloads genomes for which I’m never going to use (like fly), what gives? Why cant I chose what to install?</b>  
	<b>A –</b> We have tried to streamline the installation as much as possible so we just install everything by default. Maybe in the future we will go through and make a more thorough installer with more fine-grained options but for now this is how it is. You can modify the installer script if you are savy enough to remove those commands – that’s up to you. This is why we provide all the code used for everything, so you can tinker with the pipeline how you want. If you break it though… The other option is wait till everything installs and then delete what you don’t want – the HOMER install directories are quite large for sure.

    For now its an all or nothing installer. Future updates should address being able to add/remove aspects of the installation without having to wipe the slate clean and start over. For a 1st release though, we are happy with this behaviour for now.
11. <b>Q – If I get ChIP-AP working on my institutes computing cluster SuperUltraMegaVoltron, can you include the guide on the wiki?</b>  
	<b>A –</b> If you’re computing cluster is called SuperUltraMegaVoltron then most certainly YES we can work together to get your installation guide up on the wiki simply to say we have it working on SuperUltraMegaVoltron (but does it have a V-MAX option??). Even if It doesn’t have an exciting name as that and is simply called BigPuddle, we will still work with you to get the installation guide up on the wiki. We want ChIP-AP to become the go-to tool for analyses. So we will work towards that goal with anyone who is willing to help and collaborate.
12. <b>Q – Is there a difference in the output between the GUI and command line versions of ChIP-AP?</b>  
	<b>A –</b> Nope!  We have designed ChIP-AP to be just as functional from the GUI or command line.  Any feature accessible on the command line is also accessible from the GUI. However, with great power comes great responsibility! If you modify the settings table without knowing what you’re doing you will break the pipeline run and get wrong results. 
    
    This option, while available for modification in the GUI, really is intended for the most advanced of users. Even most “bioinformaticians” wont know how to take full advantage of this functionality and will break the pipeline more times than improving its efficiency.  So, the ability is there, but don’t use it unless you REALLLLLYYYYY know what you’re doing. This really goes for everyone, biologist or bioinformatician.

<br>