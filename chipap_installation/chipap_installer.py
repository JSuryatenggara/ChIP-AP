#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


# script_version = '2.1'    - When installing ChIP-AP in its own environment, the environment name will now be the name
#                               user typed in when prompted in the terminal, instead of "python3_chipap"
#                           - Replaces the "name:" and "prefix:" lines in the .yml file(s) if they already exist
#                           - Adds the "name:" and "prefix:" lines into the .yml file(s) if they did not yet exist
#                           - Motif enrichment analysis suite meme (version 5.0.5) has been added to both .yml files


import shutil
import subprocess
import sys
import os


# Function to search the input_file for a line containing old_string, then replace that whole line with the new_string
def replace_path(input_file, old_string, new_string):
    old_file = open(input_file, "r")

    new_file = ""

    for line in old_file.readlines():
        if old_string in line:
            new_line = new_string + "\n"
        else:
            new_line = line
        
        new_file += new_line
        old_file.close()

    writing_file = open(input_file, "w")
    writing_file.write(new_file)
    writing_file.close()


conda_file = shutil.which('conda') # Get the full path to anaconda 3
conda_dir = '{}/anaconda3'.format(conda_file.split('/anaconda3/')[0]) # Get the parent directory of anaconda 3

while True: # Infinite loop until broken by user action
    user_answer = input('Do you want to install ChIP-AP in a specific environment? (Y/N) ') # Proceed only after the user provided Y or N answer to this question

    if user_answer.lower() == 'y': # If user is installing ChIP-AP in a newly created environment
        user_environment = input('Please type in the name of your ChIP-AP environment: ') # Ask what the user wants the newly created conda environment to be named with
        chipap_env_dir = '{}/envs/{}'.format(conda_dir, user_environment)
        prefix_string = 'prefix: {}'.format(chipap_env_dir) # Path to the directory in which the newly created environment components are installed
        name_string = 'name: {}'.format(user_environment) # The name of the newly created environment
        input('ChIP-AP will be installed in {} environment. Press ENTER to continue'.format(user_environment)) # Notify the user
        break # Dialogue finished. Proceed with the installation.

    if user_answer.lower() == 'n': # If user is installing ChIP-AP in the base environment (where Anaconda 3 is installed)
        user_environment = 'base' # The default name of the environment where Anaconda 3 is installed
        chipap_env_dir = '{}'.format(conda_dir)
        prefix_string = 'prefix: {}'.format(chipap_env_dir) # Path to the directory in which the newly created environment components are installed
        name_string = 'name: {}'.format(user_environment) # The name of the newly created environment
        input('ChIP-AP will be installed in base environment. Press ENTER to continue') # Notify the user
        break # Dialogue finished. Proceed with the installation.

    else: # If user typed in anything else other than Y or N (case-insensitive)
        print('Invalid input. Please try again.') # Prompt them to try again


if sys.platform == "linux" or sys.platform == "linux2": # Execute this if the installation platform is Linux
    shutil.copy('{}/chipap_env_linux.yml'.format(sys.path[0]), '{}/chipap_{}.yml'.format(sys.path[0], user_environment)) # Copy paste the provided chipap_env_linux.yml file
    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'r') as original_yml: yml_contents = original_yml.readlines() # Read the copy pasted file contents as list of lines

    new_yml_contents = []    
    for yml_contents_row in yml_contents:
        if 'name: ' not in yml_contents_row and 'prefix: ' not in yml_contents_row:
            new_yml_contents.append(yml_contents_row) # Append all lines into the list except the line that defines the name and the directory of the environment

    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'w') as modified_yml: modified_yml.write('{}\n'.format(name_string) + ''.join(new_yml_contents) + '{}\n'.format(prefix_string)) # Insert the environment name as the first item in the list and the environment directory (prefix) as the last item in the list


elif sys.platform == "darwin": # Execute this instead if the installation platform is MacOS
    subprocess.run('/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"', shell = True)
    subprocess.run('brew install wget', shell = True)
    shutil.copy('{}/chipap_env_macos.yml'.format(sys.path[0]), '{}/chipap_{}.yml'.format(sys.path[0], user_environment)) # Copy paste the provided chipap_env_macos.yml file
    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'r') as original_yml: yml_contents = original_yml.readlines() # Read the copy pasted file contents as list of lines

    new_yml_contents = []    
    for yml_contents_row in yml_contents:
        if 'name: ' not in yml_contents_row and 'prefix: ' not in yml_contents_row:
            new_yml_contents.append(yml_contents_row) # Append all lines into the list except the line that defines the name and the directory of the environment

    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'w') as modified_yml: modified_yml.write('{}\n'.format(name_string) + ''.join(new_yml_contents) + '{}\n'.format(prefix_string)) # Insert the environment name as the first item in the list and the environment directory (prefix) as the last item in the list
    subprocess.run('pip install pandas', shell = True)


if user_answer.lower() == 'y': # Execute this if user is installing ChIP-AP in a newly created environment
    subprocess.run('conda env create -f {}/chipap_{}.yml'.format(sys.path[0], user_environment), shell = True)

if user_answer.lower() == 'n': # Execute this if user is installing ChIP-AP in the base environment (where Anaconda 3 is installed)
    subprocess.run('conda env update -f {}/chipap_{}.yml'.format(sys.path[0], user_environment), shell = True)

subprocess.run('chmod +x {}/chipap_scripts/*'.format(sys.path[0]), shell = True) # Mark all scripts in folder "chipap_scripts" as executable

subprocess.run('chmod +x {}/homer_genome_update.sh'.format(sys.path[0]), shell = True) # Mark the script "homer_genome_update.sh" as executable


# Find every instance of idr.py script in the ChIP-AP installation conda environment and put them all into a list
popen_idr = subprocess.Popen('find {} -name "idr.py"'.format(chipap_env_dir), shell = True, stdout = subprocess.PIPE)
idr_full_path_list = popen_idr.communicate()[0].decode("utf-8").split()

for idr_full_path in idr_full_path_list: # For every instance of idr.py script found above
    print("Replacing installed stock IDR script ({}) with ChIP-AP's patched IDR script ({}/idr.py)".format(idr_full_path, sys.path[0]))
    shutil.copy('{}/idr.py'.format(sys.path[0]), idr_full_path) # Replace the "broken" stock IDR script with the patched one provided by ChIP-AP installation package


genome_folder_full_path = "genome_folder_full_path = '{}/genomes'".format(sys.path[0]) # Define the new path to genome folder based on installation directory

dashboard_full_path = '{}/chipap_scripts/chipap_dashboard.py'.format(sys.path[0]) # Define the script file in which the old path to the genome folder is to be replaced
wizard_full_path = '{}/chipap_scripts/chipap_wizard.py'.format(sys.path[0]) # Define the script file in which the old path to the genome folder is to be replaced

replace_path(dashboard_full_path, 'genome_folder_full_path =', genome_folder_full_path) # Replace the old path to the genome folder with the new one
replace_path(wizard_full_path, 'genome_folder_full_path =', genome_folder_full_path) # Replace the old path to the genome folder with the new one


homer_genome_update_script_full_path = '{}/homer_genome_update.sh'.format(sys.path[0]) # Define the script file in which the old path to the genome folder is to be replaced

hg38_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install hg38'.format(chipap_env_dir) # Define the new path to homer based on installation directory
hg19_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install hg19'.format(chipap_env_dir) # Define the new path to homer based on installation directory
mm9_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install mm9'.format(chipap_env_dir) # Define the new path to homer based on installation directory
mm10_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install mm10'.format(chipap_env_dir) # Define the new path to homer based on installation directory
sacCer3_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install sacCer3'.format(chipap_env_dir) # Define the new path to homer based on installation directory
dm6_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install dm6'.format(chipap_env_dir) # Define the new path to homer based on installation directory

replace_path(homer_genome_update_script_full_path, '-install hg38', hg38_homer_genome_command_line) # Replace the old path to homer with the new one
replace_path(homer_genome_update_script_full_path, '-install hg19', hg19_homer_genome_command_line) # Replace the old path to homer with the new one
replace_path(homer_genome_update_script_full_path, '-install mm9', mm9_homer_genome_command_line) # Replace the old path to homer with the new one
replace_path(homer_genome_update_script_full_path, '-install mm10', mm10_homer_genome_command_line) # Replace the old path to homer with the new one
replace_path(homer_genome_update_script_full_path, '-install sacCer3', sacCer3_homer_genome_command_line) # Replace the old path to homer with the new one
replace_path(homer_genome_update_script_full_path, '-install dm6', dm6_homer_genome_command_line) # Replace the old path to homer with the new one

subprocess.run('bash {}/homer_genome_update.sh'.format(sys.path[0]), shell = True) # Run the provided homer_genome_update.sh to update required HOMER genomic databases


path_to_chipap_scripts = 'PATH=$PATH:{}/chipap_scripts'.format(sys.path[0]) # Define the path to the directory "chipap_scripts"

if sys.platform == "linux" or sys.platform == "linux2":  # Execute this instead if the installation platform is Linux
    path_exist = 0
    with open(os.path.expanduser('~/.bashrc'), "r+") as file:
        for line in file:
            if '{}\n'.format(path_to_chipap_scripts) == line: # Look for the path to the directory "chipap_scripts"
               path_exist = 1 # If it is already there
               break # Stop looking and proceed with the next process

        if path_exist == 0: # If it is not found until the end of the lines
            file.write('{}\n'.format(path_to_chipap_scripts)) # Write the path to the directory "chipap_scripts" so the scripts within can be accessed from any directory
    subprocess.run('exec bash', shell = True) # reload bashrc


elif sys.platform == "darwin":  # Execute this instead if the installation platform is MacOS
    subprocess.run('cp ./gem_macOS.zip {}/'.format(chipap_env_dir), shell = True)
    # shutil.copy('./gem_macOS.zip'.format("./",chipap_env_dir))
    subprocess.run('unzip -o {}/gem_macOS.zip'.format(chipap_env_dir), shell = True)

    if os.path.isfile('~/.bash_profile'): # If the computer already has the .bash_profile file
        print (".bash_profile exist") # Notify the user
        path_exist = 0
        with open(os.path.expanduser('~/.bash_profile'), "r+") as file: # Open the .bash_profile file
            for line in file:
                if '{}\n'.format(path_to_chipap_scripts) == line: # Look for the path to the directory "chipap_scripts"
                    path_exist = 1 # If it is already there
                    break # Stop looking and proceed with the next process

            if path_exist == 0: # If it is not found until the end of the lines     
                cwd = os.getcwd()
                file.write('export PATH="{}/chipap_scripts:$PATH"'.format(cwd)) # Write the path to the directory "chipap_scripts" so the scripts within can be accessed from any directory
    else: # If the computer does not have the .bash_profile file
        print (".bash_profile doesnt exist, making it") # Notify the user
        with open(os.path.expanduser('~/.bash_profile'), "w+") as file: # Make a new  .bash_profile file
            cwd = os.getcwd()
            file.write('export PATH="{}/chipap_scripts:$PATH"'.format(cwd)) # Write the path to the directory "chipap_scripts" so the scripts within can be accessed from any directory
    subprocess.run('source ~/.bash_profile', shell = True) # reload bash_profile
