## ChIP-AP Setup Guide for MacOS v1

This installation guide was performed on MacOS Catalina (10.15.7) but should be compatible with all versions 10.15 (Catalina) and later. Please note you will need administrator installation access for this installation, there is no way around that for macOS installations unfortunately. If this is your personal computer/laptop, this is just the password you use for any installation/login as normal. In the following tutorial, we will ignore any steps where you are asked for your password (as it may/will appear multiple times), when it does appear just enter the password and continue with the remainder of the instructions.

Additionally, you need to make sure your Mac is running an Intel processor. ChIP-AP will not currently run on the new Apple Silicon (AS) chips (such as the M1 or later). This is because the underlying software running ChIP-AP hasn’t been updated yet for AS chips as per our understanding. Apple’s Rosetta 2 might be able to run ChIP-AP but we cannot guarantee that it will do so flawlessly or produce consistent results, so we use the blanket statement of “it won’t run” to avoid ambiguity till we can be confident in its running (also we don’t have an M1 Mac to test it
thoroughly – any volunteers to test this for us?).

To make sure you’re Mac has the right processor, go to the Apple menu in the top left (a), click “About This Mac” (b) and then check that you have an “Intel” processor (c).

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_1.PNG>

Newer versions of MacOS have updated security to restrict you installing software from anywhere other than the AppStore – because Apple is the gatekeeper of course. So, you may need to turn this restriction off before proceeding in the installation. To do this, go to the Apple menu in the top left (a), click “System Preferences” (b), “Security and Privacy” (c).

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_2.PNG>

Next, go to “General” (a), then make sure that you can make changes by ensuring the lock is open (b), then finally click Allow apps downloaded from “App Store and identified developers” (c).

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_3.PNG>

Now you can proceed with the installation and should be ok.

A reminder about the installation requirements. In terms of storage space, the installation of ChIP-AP and all its required generated genome files etc... will take up ~50-60Gb alone – this is NOT including any of your chip or analysis files. So, make sure you set up this on a drive/computer with enough storage space to accommodate these heavy requirements. Off-loading some of the installation directories to an external drive for example is possible (using hard/symbolic links) but is an advanced topic and is beyond the scope of this tutorial. Again, google is your friend here.

<br>

### Prerequisite Installation
Before installing ChIP-AP, you need to make sure that Anaconda 3 is installed. There are thorough guides online but we will cover the basic steps here.

1. Firstly, in your web-browser of choice, search for “anaconda3 individual install”
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_4.PNG>

<br>

2. Download (a) the right installer for your system (should be a 64-bit Graphical Installer [go command line installer if you’re feeling gung-ho!]) (b) and save it in your Downloads folder (c). When the download completes, you can press the magnifying glass icon to be taken to where the installer was downloaded (d)
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_5.PNG>
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_6.PNG>

<br>

3. Next, double click on the installer package to start it
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_7.PNG>

<br>

4. At the first window, click “Continue” on the first sheet
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_8.PNG>

<br>

5. Click “Continue” again
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_9.PNG>

<br>

6. Accept the license and End User agreements and click “Continue”
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_10.PNG>

<br>

7. Agree to the agreement
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_11.PNG>

<br>

8. Next click “Install.” If you are more familiar with setting up Anaconda you can change the installation location, or if you want to install it to an external drive if your built-in storage isn’t large enough, then do so as required for your system.
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_12.PNG>

<br>

9. Once complete, you will be asked about the PyCharm IDE environment as well – click “Continue”
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_13.PNG>

<br>

10. On the next window, click “Close”
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_14.PNG>

<br>

11. Great! Anaconda 3 is now installed. If you are asked to move the installer to the bin, do so since you will no longer need it, this is only done on some versions of Safari so if you’re not asked about it don’t worry. 

    If you installed Anaconda from the command line you must choose “yes” to initialize anaconda and then close the Terminal entirely and reopen it before proceeding. If you don’t you the remaining installation steps will fail to complete properly.

<br>

12. The next thing we need to make sure to install is Java since some newer versions of macOS don’t come with it pre-installed – thanks Apple!!! Ok so you need to go to (Java SE Development Kit 16 - Downloads (oracle.com)) and download the macOS  Installer (a), accept the license terms (b), download it (c), and allow downloads (e)
    
    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_15.PNG>

    Next, start the downloaded package by double-clicking on it.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_16.PNG width=200>

    Run through the installer process accepting all defaults unless you know how to configure the installer yourself. Enter the password when prompted to. You can delete the installer when complete as its no longer needed.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_17.PNG>

<br>

### Setting up ChIP-AP
Now down to the “meat” of the installation! Actually, when we first started writing we thought ”pffttt it works fine in Linux, macOS is Unix it’ll be a breeze to adapt....” Oh boy were we wrong!!! OK let's start...

1. To start, we recommend you create a sub-folder in your “Documents” directory named “tools”

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_18.PNG>

<br>

2. Within the “tools” folder, create a sub-folder named “chip_ap.” Notice the
underscore here rather than dash, this is to simplify a number of steps
unbeknownst to you later 

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_19.PNG width=300>

<br>

3. Download the latest version of ChIP-AP from our github page (https://github.com/JSuryatenggara/ChIP-AP) and place it in the “chip_ap” folder. 

    Also, download the pre-configured genome folders required for ChIP-AP from (https://www.dropbox.com/s/ioqd3hwdahh9xon/genomes.zip) and place in the same folder. Once everything is unzipped, make sure the “chip_ap” folder structure looks as follows:

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_20.PNG>

<br>

4. Next, you need to go up 1 folder

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_21.PNG>

    Next, right mouse click on the “chip_ap” folder and the last item in the pop-up menu will be “New Terminal at Folder”

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_22.PNG>

    This will bring up the terminal window we now need to complete the installation

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_23.PNG>

<br>

5. At the terminal, type the following command and press ENTER

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_24.PNG>

<br>

6. As part of the installation process, ChIP-AP will ask a couple of questions

    The first question asked is whether you want to install ChIP-AP in its own environment. This is so ChIP-AP will sit in its own contained “capsule” and won’t be affected by other installed programs and so this will ensure that all the dependencies of ChIP-AP will be correct and met all the time. So, we highly recommend this option

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_25.PNG>

    If you answer “n”, then ChIP-AP will be installed in the base environment along with its dependencies.

    If you answer “y”, you will be asked to name the environment and then press ENTER. A confirmation statement will be printed for you to confirm with an additional ENTER press. If doing this, make sure to remember the name of the environment you entered as its very important!

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_26.PNG>

<br>

7. Once you press ENTER, you will be asked to enter your password for the next installation step. When you start typing your password NOTHING will appear. The text entry is still working but it wont show you any character entered for security reasons. So, make sure you type in your password correctly since you can’t see what you’re doing. When entered press ENTER to continue.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_27.PNG>

<br>

8. You will be then asked to confirm all details with an ENTER press on the keyboard.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_28.PNG>

    The installation will now start and will take a while to run as a lot of packages will need to be installed, so this is highly dependent on your internet and cpu speeds as to how long it will take. Take note of the environment name you set. For now, go take a break but remember, no eating in the lab!!

<br>

9. Once the installation completes close the terminal entirely and open up a new terminal window by searching for it in spotlight (a), and typing in “terminal” (b), and selecting it in the results window (c)

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_29.PNG>

<br>

10. Depending on your proficiency, there are 3 ways to use ChIP-AP. If you setup ChIP-AP in an environment in step 6, once you open the terminal, you will need to type the following command (below) where xxxx is the name of the environment as you defined in step 6 before running ChIP-AP EVERYTIME. Otherwise, this command isn’t necessary.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_30.PNG>

    <br>
    <br>

    - To use ChIP-AP using the command line, refer to the documentation on our github (xxx) for full details on how to setup a run with appropriate flags/parameters (For advanced/proficient users).

        <br>

    - For most non-bioinformatician users of ChIP-AP, we recommend using the wizard. To use this, at the opened terminal command line type:

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_31.PNG>

        ... and the first window of the wizard will appear (below). If running ChIP-AP for the first time, it may take 1-2 minutes for the first window to appear as the pipeline checks for any missing dependencies and installs them in the background, so be patient till everything appears. There will be a message printed to the command line if this is the case if you notice it, otherwise just wait for the initial check to complete. This delay will only be for the first time running either the dashboard or the wizard though.
        
        Once the wizard appears, the wizard will ask you questions sequentially to acquire all the data it needs for a successful run. Below also is a snapshot of the windows you will see throughout the wizard process.

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_32.PNG>

        <br>
        <br>

    - For users who are more experienced with ChIP-AP, we recommend using the dashboard. At the command line, simply type:

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_33.PNG>

        ... and the dashboard will appear (below) wherein you can fill in all required information. 

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_34.PNG>

        NOTE: For you to use the dashboard you need to have a minimum screen resolution of at least 1920x1080. If it is less than this, then all the dashboard elements will not appear on the screen and hence will be unusable. So either use an external monitor with a resolution of 1920x1080 or larger, or the wizard.

<br>

11. Once the run is complete, please refer to our github wiki for the tutorial on navigating and interpreting the results.

<br>

### Subsequent Running of ChIP-AP
OK so now that everything is set up to run ChIP-AP each time is relatively simple. Open a new terminal window... 

If you setup ChIP-AP in its own environment then you need to write the following command first:

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_35.PNG>

where xxxxx is the name of the environment you created in the installation.

The next command you need to write at the command line will start chipap depending on what mode you want:

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/MacOS/MacOS_36.PNG>

NOTE: For you to use the dashboard you need to have a minimum screen resolution of at least 1920x1080. If it is less than this, then all the dashboard elements will not appear on the screen and hence will be unusable. So either use an external monitor with a resolution of 1920x1080 or larger, or the wizard.

<br> 

### Installation Troubleshooting
So MacOS seems to be picky with a number of things and this is getting worse with each new version of macOS released.  Despite out best efforts to streamline the installation process, depending on your systems configuration and previous installations there may be issues.  We will keep updating this list as we come across issues users are coming across.

1. ChIP-AP fails to load with a “pandas not found” error or similar message.
    
    We have specifically told macOS to install this but for some reason it sometimes doesn’t register.  So try re-running again and it should work.  In our experiences if it complains, an immediate re-run will work fine.