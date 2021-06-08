## ChIP-AP Setup Guide for Windows v1

Of all the installation guides we’ve prepared, the Windows guide is by far the longest and most involved with most command line dabbling till set up.  This is because the software used for ChIP-AP was never destined to run on Windows but there is a way around it if you’re tech savy enough.  Fortunately for you, we think we are 😉 

This installation guide was performed on a Windows 10 installation running v20H2, although will work with v1903 or later (Homer or Professional). Basically, as long as you have access to install the Linux Subsystem you will be fine.  Older versions of Windows, like Windows 7/8, are not compatible. For this installation, you will need administrator access to your machine to install and make modifications as well as an active Microsoft Store Account. If you are setting it up on your own personal laptop, this should be fine as it’s your usual login password.  If this is a computer owned/maintained by your institute’s IT department, you will either need them to set it up for you or give you temporary administrator access to set it up (which they can if they are nice).

<br>

### Prerequisite Installation

1. To check which version of Windows you are running, go to the search icon (a), and type “About” in search segment (b), and then click “About your PC” in (c).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_1.PNG>

<br>

2. In the window that appears, you will see which version of Windows you are running (a).  Your version needs to be 1903 or greater.  In our example we are running 20H2. The 19/20 refers to the year of release (2019/2020).  So, if your “year” is 2019 onwards you should be set to proceed.  If not, then you need to update to the latest version of Windows 10 to continue.  

    You will also see how much RAM is installed, as a minimum you need 8Gb, but 16Gb or more is recommended.  You will also need at least 60Gb of storage space on your SSD/HDD just to setup ChIP-AP, not including any additional space you need for your actual analysis, which we recommend as another 100Gb or so.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_2.PNG>

<br>

3. Now that were sure you have the right build of Windows installed; we need to setup the Linux Subsystem component of windows to allow us to install everything.  To do this, you need to go back to the search icon (a) and search for “Windows Features” (b), and select it (c).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_3.PNG>

<br>

4. In the window that appears, scroll down until you find the “Windows Subsystem for Linux” option and tick the box next to it (a), then press OK (b).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_4.PNG>

<br>

5. When you press OK, the following window will appear (a) and Windows will go through to setup what it needs to.  This can take a while depending on the speed of your machine.  When complete you MUST restart your computer before you can continue the rest of the installation process (b). 


    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_5.PNG>

<br>

6. Alright now that we have the Linux subsystem component installed, we now need to actually install Linux as a subsystem (ok that sounds confusing but it actually is right).  To do this you need to go to the search icon again (a) and search for the “Store” (b), and select it (c). 

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_6.PNG>

<br>

7. Once loaded, you need to search (a) for “Ubuntu 20.04 LTS” (b) and press ENTER or select the search result.


    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_7.PNG>

<br>

8. In the following window, you will need to “Get” this Linux distriubtion (a).  You will then be prompted to enter your login and purchase details.  This is FREE to download and install so there is no charge but you still need to go through the steps to get it added to your account be it add password/PIN or other means. Ubuntu 16.04/18.04 are also compatible and will work just fine.


    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_8.PNG>

<br>

9. Once Ubuntu 20.04 LTS is added to your account you then need to click the “Install” icon (a).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_9.PNG>

<br>

10. Finally once its installed, you need to “Launch” (a) the installation. 

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_10.PNG>

<br>

11. In the termnal window that appears, you’re first message will be “Installing, this may take a few minutes”… which it will….

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_11.PNG>

<br>

12. After a few moments you will be asked to enter a new username (which can be different than your windows login name) and then press ENTER.  If after a few minutes nothing appears at all, press ENTER once, this should then bring up the following username question if it didn’t appear (for whatever reason that may happen).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_12.PNG>

<br>

13. Next, you will be asked for a new password (a).  This can be different than your windows password. When you start typing, NOTHING will appear for security reasons.  So make sure you enter the password exactly as you want it.  Press ENTER when done.  You will be asked to re-enter the password to confirm (b).  Press ENTER again when done.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_13.PNG>

<br>

14. This now completes the Linux subsystem installation and setup.  You will see information comparable to the following to validate everything worked fine.  Once you see the green final prompt (c), everything is good to go and you can close the Linux terminal (d).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_14.PNG>

<br>

15. So, up until now, we have installed the Linux subsystem to allow us to install all the necessary programs under the hood.  There are now 2 more things to install before we can get to installing ChIP-AP.  The first thing to install is MobaXterm.  To do this, open up your browser of choice and search for “MobaXterm”

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_15.PNG>

<br>

16. Next, navigate to the download section and download the “Home Edition” which is free.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_16.PNG>

<br>

17. On the next web-page, since you have administrator rights, pick the installer edition. 

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_17.PNG>

<br>

18. When you choose to download it, you will be downloading a *.zip folder.  Unzip this in a location of your choosing.  You will then need to run the *.msi installer.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_18.PNG>

<br>

19. Next, run through the installer accepting all default settings (unless you know what you’re doing and want to configure it yourself). 

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_19.PNG>

<br>

20. Next, run MobaXterm by going to the search icon (a) and search for “MobaXterm” (b) and select it (c).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_20.PNG>

<br>

21. The first time you run MobaXterm, the initial splash screen will take a little longer to set everything up.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_21.PNG width=400>

<br>

22. You will then be asked to give MobaXterm access through the firewall, which you must accept to (a).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_22.PNG>

<br>

23. Alright, the next couple steps are a little more complicated so take special care. We will need to switch between your web-browser and MobaXterm.  So first, your web-browser of choice, search for “anaconda3 individual install”

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_23.PNG>

<br>

24. Next, go to the bottom of the page by clicking, “Download” (a) and get to the stage where you can select the right installer (b).

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_24.PNG>

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_25.PNG>

<br>

25. Next, right click on the “64-Bit (x86) Installer” under Linux and select “copy link.”

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_26.PNG>

<br>

26. Now, go back to MobaXterm and double mouse click “WSL-Ubuntu-20.04” in the top left.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_27.PNG>

<br>

27. This will open the panel shown below

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_28.PNG>

<br>

28. Next, type the following commands, 1 line at a time, pressing ENTER after the end of each line.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_29.PNG>

    What you should now see is something as follows…

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_30.PNG>

<br>

29. So, what we did is made a new directory (folder) named “tools”, changed directory (cd) into it.  Then we made another new directory named “chipap” and then changed directory (cd) into that.  The next step is to download the anaconda installation file and run it in Linux.  Now, at the command prompt, type “wget” and then paste the link path from the web-browser above (right-click and go to paste).  It will look something like what’s below.  Now press ENTER.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_31.PNG>

<br>

30. The installer file will then download and give you back the command line prompt.  You then need to type “bash Ana” and then press TAB on the keyboard and the terminal will auto-complete the rest of the command required so you don’t have to complete it manually.  It will look as follows (below). Then press ENTER to start the installation process.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_32.PNG>

<br>

31. Press ENTER to view the license agreement which you must agree to for installation to proceed.  To progress through each page of the agreement press the SPACE bar - read it, it if you want…

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_33.PNG>

    Once you get to the end… you must type “yes” in full and press ENTER

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_34.PNG>

<br>

32. Next, you will be asked where to install anaconda3.  Unless you are an advanced user and know how to configure anaconda properly, accept defaults and press ENTER again.  

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_35.PNG>

    Anaconda will then go ahead and install itself.  This will take a few minutes.  Remember though, no eating in the lab while this is working!

<br>

33. When complete, you will be asked to initialize Anaconda3, type “yes” and press ENTER

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_36.PNG>

<br>

Great! Anaconda3 is now installed.  Now you must close and re-start MobaXterm before the next steps of actually installing ChIP-AP!!!

Please read carefully the final few lines if it says something along the lines of “you have chosen to not initialize anaconda….” You need to re-run the anaconda installation using the command below making sure to add the “-u” at the end.

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_37.PNG>

<br>

If you don’t make sure to initialize anaconda post-installation, nothing down-stream of this will work.



### Setting up ChIP-AP
Now down to the “meat” of the installation!
1. First thing, you need to open a new terminal of MobaXterm as before.  Double click on Ubuntu 20.04.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_38.PNG>

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_39.PNG>

<br>

2. In the terminal window, you need to type the following commands 1 per line followed by ENTER.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_40.PNG>

<br>

3. You now need to download the latest version of ChIP-AP from our github page by typing “wget” followed by a space, then pasting the following url and press ENTER after to initiate the download: https://github.com/JSuryatenggara/ChIP-AP


    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_41.PNG>

    This shouldn’t take long to download.  
    
    Next, we need to download the chipap genomes required for the run.  As before, typing “wget” followed by a space, then pasting the following url and press ENTER after to initiate the download: https://www.dropbox.com/s/ioqd3hwdahh9xon/genomes.zip

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_42.PNG>

    This WILL take a long time to download

<br>

4. Once the previous download has finished, enter the following command followed by ENTER.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_43.PNG>

    This will ask you to enter the password you set in the beginning when we installed Linux.

<br>

5. Next, enter the following commands 1 at a time followed by ENTER.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_44.PNG>

    Next…type the following command and press ENTER, to start the installation process

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_45.PNG>

<br>

6. As part of the installation process, ChIP-AP will ask a couple of questions

    The first question asked is whether you want to install ChIP-AP in its own environment.  This is so ChIP-AP will sit in its own contained “capsule” and won’t be affected by other installed programs and so this will ensure that all the dependencies of ChIP-AP will be correct and met all the time.  So, we highly recommend this option

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_46.PNG>

    If you answer “n”, then ChIP-AP will be installed in the base environment along with its dependencies.

    If answer “y”, you will be asked to name the environment and then press ENTER.  A confirmation statement will be printed for you to confirm with an additional ENTER press.  The installation will take a while to run as a lot of packages will need to be installed, so this is highly dependent on your internet and cpu speeds as to how long it will take. Take note of the environment name.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_47.PNG>

<br>

7. Once the installation completes, close MobaXterm and reopen it.  Open up a new Ubuntu 20.04 terminal as before.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_48.PNG>

<br>

8. Depending on your proficiency, there are 3 ways to use ChIP-AP. If you setup ChIP_AP in an environment in step 6, once you open the terminal, you will need to type the following command (below) where xxxx is the name of the environment as you defined in step 6 before running ChIP-AP EVERYTIME. Otherwise, this command isn’t necessary.

    <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_49.PNG>

    <br>
    <br>

    - To use ChIP-AP using the command line, refer to the documentation on our github (xxx) for full details on how to setup a run with appropriate flags/parameters (For advanced/proficient users).

    <br>

    - For most non-bioinformatician users of ChIP-AP, we recommend using the wizard.  To use this, at the opened terminal command line type:

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_50.PNG>

        ... and the first window of the wizard will appear (below).  The wizard will ask you questions sequentially to acquire all the data it needs for a successful run.  Below also is a snapshot of the windows you will see throughout the wizard process. 

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_51.PNG>

        <br>
        <br>

    - For users who are more experienced with ChIP-AP, we recommend using the dashboard.  At the command line, simply type 
    
        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_52.PNG>

        ... and the dashboard will appear (below) wherein you can fill in all required information.  NOTE:  For you to use the dashboard you need to have a minimum screen resolution of at least 1920x1080.  If it is less than this, then all the dashboard elements will not appear on the screen and hence will be unusable.  So either use an external monitor with a resolution of 1920x1080 or larger, or the wizard. 

        <img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_53.PNG>

    <br>

8. Once the run is complete, please refer to our github wiki for the tutorial on navigating and interpreting the results.

<br>
### Subsequent Running of ChIP-AP 
Ok so now that everything is set up to run ChIP-AP each time is relatively simple.  
Open a new MobaxTerm terminal window…

If you setup ChIP-AP in its own environment then you need to write the following command first

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_54.PNG>

where xxxxx is the name of the environment you created in the installation

The next command you need to write at the command line will start chipap depending on what mode you want

<img src=https://raw.githubusercontent.com/JSuryatenggara/ChIP-AP/storage/images/Windows/Windows_55.PNG>

NOTE:  For you to use the dashboard you need to have a minimum screen resolution of at least 1920x1080.  If it is less than this, then all the dashboard elements will not appear on the screen and hence will be unusable.  So either use an external monitor with a resolution of 1920x1080 or larger, or the wizard.