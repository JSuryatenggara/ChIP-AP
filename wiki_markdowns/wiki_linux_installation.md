# ChIP-AP
## Linux Installation Guide v1
This installation guide was performed on the Ubuntu VM provided for your use (i.e., if you use the provided VM, this is how it was set up). Other Ubuntu based variants may require some slight modifications to certain steps (as some base dependencies might differ between variants), so keep that in mind. Unfortunately, it is impossible for us to make a guide for every variant out there. Also, some installations steps might have issues depending on your current system setup and configuration. This  installation guide was done on a fresh install of Ubuntu 20.04 LTS but we have tested it on the 16.04 and 18.0x LTS versions as well. To resolve any issues that may arise, Google (not Bing, sorry Microsoft) and Stack Exchange are your best friends really, so look online and you should be able to find a solution that works for you.

A point of note: Ideally, for you to install ChIP-AP you need to have administrator access to install the required files (so this is if it’s your own personal laptop, this is the ideal situation). It is possible to install in shared computing
clusters but there are usually restrictions as to what can be installed, so the installation will be trickier unfortunately. You may need to get assistance from the maintaining IT team to help you set it up if in this situation. We have done
our best to ensure that all can be installed without outside help, but depending how computer systems are setup, we cannot guarantee this. As a point of note, this tutorial doesn’t require administrator access to install anything, but your computer setup might vary. So, attempt to install in your base directory first without administrator help, but if you get stuck with administrator issues, then you will need to contact your institutes IT team for them to help you install the required dependencies.

A reminder again about the installation requirements. In terms of storage space, the installation of ChIP-AP and all its required generated genome files etc... will take up ~50-60Gb alone – this is NOT including any of your chip or analysis files. So, make sure you set up this on a drive/computer with enough storage space to accommodate these heavy requirements. Off-loading some of the installation directories to an external drive for example is possible (using hard/symbolic links) but is an advanced topic and is beyond the scope of this tutorial. Again, google is your friend here.

<br>

### Prerequisite Installation
Before installing ChIP-AP, you need to make sure that Anaconda 3 is installed. There are thorough guides online but we will cover the basic steps here.

1. Firstly, in your web-browser of choice, search for “anaconda3 individual install”

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_1.PNG>
    
<br>

2. Download the right installer for your system (should be a 64-bit linux installation) and save it in your Downloads folder (a) → (b) → (c)

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_2.PNG>
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_3.PNG>
    
<br>

3. Navigate to the folder where the installer “.....sh” file is downloaded, in this example, it was downloaded to the “Downloads” folder. 
    
    So, open up the file explorer (a) by clicking the right icon. Navigate to “Downloads” (b), and right mouse click on the white space of the folder. Next select “Open in Terminal” (c).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_4.PNG>
    
<br>
    
4. The following terminal window should now appear.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_5.PNG>
    
<br>

5. Next, type the following command

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_5.5.PNG>

    ...where the filename entered after “bash” is the exact same filename as the 1.3 downloaded from anaconda. For ease of writing, you can write “bash Ana” and then press the TAB key on your keyboard and the terminal should auto-complete the rest of the filename for you – so you don’t have to do it manually. Once done, press ENTER
    
<br>

6. Press ENTER to view the license agreement which you must agree to for installation to proceed. To progress through each page of the agreement press the SPACE bar - read it, if you want...

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_6.PNG>

    Once you get to the end... you must type “yes” in full and press ENTER

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_7.PNG>
    
<br>

7. Next, you will be asked where to install anaconda3. Unless you are an advanced user and know how to configure anaconda properly, accept defaults and press ENTER again.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_8.PNG>

    Anaconda will then go ahead and install itself. This will take a few minutes. Remember though, no eating in the lab while this is working!
    
<br>

8. When complete, you will be asked to initialize Anaconda3, type “yes” and press ENTER

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_9.PNG>
    
<br>

9. Great! Anaconda 3 is now installed. Now you must close and re-start the terminal before the next steps.

<br>

### Setting Up ChIP-AP
Now down to the “meat” of the installation!

1. Download the latest version of ChIP-AP from our github page: (https://github.com/JSuryatenggara/ChIP-AP) and place it in its own folder in the location of your choosing.

    Our recommendation is as follows - in your home directory, create a new folder named “tools” and create a sub-folder named “chip_ap” (notice the underscore here rather than dash, this is to simplify a number of steps later unbeknownst to you).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_10.PNG>

    Also, while you’re at it, download the pre-configured genome folders required for ChIP-AP from (https://www.dropbox.com/s/ioqd3hwdahh9xon/genomes.zip). Unzip the file and folders and put them in a “genomes” folder. They are required for running ChIP-AP successfully.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_11.PNG>
    
<br>

2. Make sure you unzip everything and your folder structure looks similar to the following

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_12.PNG>
    
<br>

3. Next, right mouse-click on the white-space of the folder, and pull up the context menu, and click on “Open in Terminal” to bring up a new terminal window.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_13.PNG>
    
<br>

4. Now to install ChIP-AP and its dependencies, from within the extracted folder

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_14.PNG>

    ..type the following command and press ENTER

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_14.5.PNG>
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_15.PNG>
    
<br>

5. As part of the installation process, ChIP-AP will ask a couple of questions.
    
    The first question asked is whether you want to install ChIP-AP in its own environment. This is so ChIP-AP will sit in its own contained “capsule” and won’t be affected by other installed programs and so this will ensure that all the dependencies of ChIP-AP will be correct and met all the time. So, we highly recommend this option

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_16.PNG>

    If you answer “n”, then ChIP-AP will be installed in the base environment along with its dependencies.

    If answer “y”, you will be asked to name the environment and then press ENTER. A confirmation statement will be printed for you to confirm with an additional ENTER press. The installation will take a while to run as a lot of packages will need to be installed, so this is highly dependent on your internet and cpu  peeds as to how long it will take. Take note of the environment name you set.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_17.PNG>
    
<br>

6. Once the installation completes, open up the terminal

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_18.PNG>
    
<br>

7. Depending on your proficiency, there are 3 ways to use ChIP-AP. If you setup ChIP_AP in an environment in step 5, once you open the terminal, you will need to type the following command (below) where xxxx is the name of the environment as you defined in step 5 before running ChIP-AP EVERYTIME. Otherwise, this command isn’t necessary.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_19.PNG>
    
    <br>
    <br>

    - To use ChIP-AP using the command line, refer to the documentation on our github (xxx) for full details on how to setup a run with appropriate flags/parameters (For advanced/proficient users).

    <br>

    - For most non-bioinformatician users of ChIP-AP, we recommend using the wizard. To use this, at the opened terminal command line type:

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_20.PNG>

        ... and the first window of the wizard will appear (below). The wizard will ask you questions sequentially to acquire all the data it needs for a successful run. Below also is a snapshot of the windows you will see throughout the wizard process.

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_21.PNG>

    <br>
    <br>

    - For users who are more experienced with ChIP-AP, we recommend using the dashboard. At the command line, simply type:

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_22.PNG>

        ... and the dashboard will appear (below) wherein you can fill in all required information.  NOTE:  For you to use the dashboard you need to have a minimum screen resolution of at least 1920x1080.  If it is less than this, then all the dashboard elements will not appear on the screen and hence will be unusable.  So either use an external monitor with a resolution of 1920x1080 or larger, or the wizard.

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_23.PNG>
    
<br>

8. Once the run is complete, please refer to our github wiki for the tutorial on navigating and interpreting the results.

<br>

### Subsequent Running of ChIP-AP

Ok so now that everything is set up to run ChIP-AP each time is relatively simple.  Open a new terminal window…

If you setup ChIP-AP in its own environment then you need to write the following command first

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_24.PNG>
    
where xxxxx is the name of the environment you created in the installation.

The next command you need to write at the command line will start chipap depending on what mode you want

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Linux/Linux_25.PNG>
    
NOTE:  For you to use the dashboard you need to have a minimum screen resolution of at least 1920x1080.  If it is less than this, then all the dashboard elements will not appear on the screen and hence will be unusable.  So either use an external monitor with a resolution of 1920x1080 or larger, or the wizard.