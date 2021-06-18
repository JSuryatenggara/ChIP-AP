#!/usr/bin/env python3
#pyright: reportUnboundVariable=false

import shutil
import subprocess
import sys
import os
import time
# try:
#     import pandas as pd
# except ImportError:
#     print("Dependency missing, downloading and installing now")
#     import pip
#     pip.main(['install', '--user', 'pandas'])
#     time.sleep(5) # Sleep for 3 seconds

import pandas as pd

def replace_path(input_file, old_string, new_string):
    old_file = open(input_file, "r")

    new_file = ""

    for line in old_file.readlines():
        # stripped_line = line.strip()
        if old_string in line:
            new_line = new_string + "\n"
        else:
            new_line = line
        
        new_file += new_line
        old_file.close()

    writing_file = open(input_file, "w")
    writing_file.write(new_file)
    writing_file.close()



conda_file = shutil.which('conda')
conda_dir = '{}/anaconda3'.format(conda_file.split('/anaconda3/')[0])

while True:
    user_answer = input('Do you want to install ChIP-AP in a specific environment? (Y/N) ')

    if user_answer.lower() == 'y':
        user_environment = input('Please type in the name of your ChIP-AP environment: ')
        chipap_env_dir = '{}/envs/{}'.format(conda_dir, user_environment)
        prefix_string = 'prefix: {}'.format(chipap_env_dir)
        name_string = 'name: {}'.format(user_environment)
        input('ChIP-AP will be installed in {} environment. Press ENTER to continue'.format(user_environment))
        break

    if user_answer.lower() == 'n':
        user_environment = 'base'
        chipap_env_dir = '{}'.format(conda_dir)
        prefix_string = 'prefix: {}'.format(chipap_env_dir)
        name_string = 'name: {}'.format(user_environment)
        input('ChIP-AP will be installed in base environment. Press ENTER to continue')
        break

    else:
        print('Invalid input. Please try again.')


if sys.platform == "linux" or sys.platform == "linux2":
    shutil.copy('{}/chipap_env_linux.yml'.format(sys.path[0]), '{}/chipap_{}.yml'.format(sys.path[0], user_environment))
    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'r') as original_yml: yml_contents = original_yml.read()
    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'w') as named_yml: named_yml.write('{}\n'.format(name_string) + yml_contents + '\n{}'.format(prefix_string))
    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'a') as prefixed_yml: prefixed_yml.write('\n{}'.format(prefix_string))

elif sys.platform == "darwin":
    subprocess.run('/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"', shell = True)
    subprocess.run('brew install wget', shell = True)
    shutil.copy('{}/chipap_env_macos.yml'.format(sys.path[0]), '{}/chipap_{}.yml'.format(sys.path[0], user_environment))
    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'r') as original_yml: yml_contents = original_yml.read()
    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'w') as named_yml: named_yml.write('{}\n'.format(name_string) + yml_contents + '\n{}'.format(prefix_string))
    with open('{}/chipap_{}.yml'.format(sys.path[0], user_environment), 'a') as prefixed_yml: prefixed_yml.write('\n{}'.format(prefix_string))
    subprocess.run('pip install pandas', shell = True)


subprocess.run('conda env create -f {}/chipap_{}.yml'.format(sys.path[0], user_environment), shell = True)

subprocess.run('chmod +x {}/chipap_scripts/*'.format(sys.path[0]), shell = True)

subprocess.run('chmod +x {}/homer_genome_update.sh'.format(sys.path[0]), shell = True)


genome_folder_full_path = "genome_folder_full_path = '{}/genomes'".format(sys.path[0])

dashboard_full_path = '{}/chipap_scripts/chipap_dashboard.py'.format(sys.path[0])
wizard_full_path = '{}/chipap_scripts/chipap_wizard.py'.format(sys.path[0])

replace_path(dashboard_full_path, 'genome_folder_full_path =', genome_folder_full_path)
replace_path(wizard_full_path, 'genome_folder_full_path =', genome_folder_full_path)



homer_genome_update_script_full_path = '{}/homer_genome_update.sh'.format(sys.path[0])

hg38_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install hg38'.format(chipap_env_dir)
hg19_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install hg19'.format(chipap_env_dir)
mm9_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install mm9'.format(chipap_env_dir)
mm10_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install mm10'.format(chipap_env_dir)
sacCer3_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install sacCer3'.format(chipap_env_dir)
dm6_homer_genome_command_line = 'perl {}/share/homer*/.//configureHomer.pl -install dm6'.format(chipap_env_dir)

replace_path(homer_genome_update_script_full_path, '-install hg38', hg38_homer_genome_command_line)
replace_path(homer_genome_update_script_full_path, '-install hg19', hg19_homer_genome_command_line)
replace_path(homer_genome_update_script_full_path, '-install mm9', mm9_homer_genome_command_line)
replace_path(homer_genome_update_script_full_path, '-install mm10', mm10_homer_genome_command_line)
replace_path(homer_genome_update_script_full_path, '-install sacCer3', sacCer3_homer_genome_command_line)
replace_path(homer_genome_update_script_full_path, '-install dm6', dm6_homer_genome_command_line)

subprocess.run('bash {}/homer_genome_update.sh'.format(sys.path[0]), shell = True)


path_to_chipap_scripts = 'PATH=$PATH:{}/chipap_scripts'.format(sys.path[0])

if sys.platform == "linux" or sys.platform == "linux2":
    # linux
    with open(os.path.expanduser('~/.bashrc'), "r+") as file:
        for line in file:
            if '{}\n'.format(path_to_chipap_scripts) == line:
               break
            else: # not found, we are at the eof
                file.write('{}\n'.format(path_to_chipap_scripts)) # append the new path to chipap_scripts
    subprocess.run('exec bash', shell = True) # reload bashrc

elif sys.platform == "darwin":
    # OS X
    subprocess.run('cp ./gem_macOS.zip {}/'.format(chipap_env_dir), shell = True)
    # shutil.copy('./gem_macOS.zip'.format("./",chipap_env_dir))
    subprocess.run('unzip -o {}/gem_macOS.zip'.format(chipap_env_dir), shell = True)

    if os.path.isfile('~/.bash_profile'):
        print (".bash_profile exist")
        with open(os.path.expanduser('~/.bash_profile'), "r+") as file:
            for line in file:
                if '{}\n'.format(path_to_chipap_scripts) == line:
                    break
                else: # not found, we are at the eof
                    cwd = os.getcwd()
                    file.write('export PATH="{}/chipap_scripts:$PATH"'.format(cwd))
    else:
        print (".bash_profile doesnt exist, making it")
        with open(os.path.expanduser('~/.bash_profile'), "w+") as file:
            cwd = os.getcwd()
            file.write('export PATH="{}/chipap_scripts:$PATH"'.format(cwd))
    subprocess.run('source ~/.bash_profile', shell = True) # reload bash_profile
