#!/usr/bin/env python
 
import os
import sys
import subprocess

import configure

# Thanks to http://patorjk.com/software/taag/  
# and http://www.chris.com/ascii/index.php?art=creatures/dragons
# for ASCII art inspiriation

pyne_logo = """\

                                  /   \       
 _                        )      ((   ))     (                          
(@)                      /|\      ))_((     /|\                          
|-|                     / | \    (/\|/\)   / | \                      (@) 
| | -------------------/--|-voV---\`|'/--Vov-|--\---------------------|-|
|-|                         '^`   (o o)  '^`                          | |
| |                               `\Y/'                               |-|
|-|                                                                   | |
| |        /\             ___           __  __             /\         |-|
|-|       /^~\           / _ \_   _  /\ \ \/__\           /^~\        | |  
| |       /^~\          / /_)/ | | |/  \/ /_\             /^~\        |-|
|-|       /^~\         / ___/| |_| / /\  //__             /^~\        | | 
| |       ^||`         \/     \__, \_\ \/\__/             ^||`        |-|  
|-|        ||                |____/                        ||         | | 
| |       ====                                            ====        |-|
|-|                                                                   | |
| |                                                                   |-|
|-|___________________________________________________________________| |
(@)              l   /\ /         ( (       \ /\   l                `\|-|
                 l /   V           \ \       V   \ l                  (@)
                 l/                _) )_          \I                   
                                   `\ /'
                                     `  
"""

def main_body():
    if not os.path.exists('build'):
        os.mkdir('build')
    makefile = os.path.join('build', 'Makefile')
    if not os.path.exists(makefile):
        rtn = subprocess.check_call(['cmake', '..'], cwd='build')
    rtn = subprocess.check_call(['make'], cwd='build')
    cwd = os.getcwd()
    os.chdir('build')
    configure.setup()
    os.chdir(cwd)

def main():
    success = False
    try:
        main_body()
        success = True
    finally:
        configure.final_message(success)

if __name__ == "__main__":
    main()
