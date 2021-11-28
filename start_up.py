# -*- coding: utf-8 -*-
"""
@author: Popova A.S. 
"""

import subprocess
import numpy as np 


program_list_1 = ["Matrix.py", "Exact", "Minors" ,"Moments.py", "Approximation.py", "Result_exact.py"]

program_list_2 = ["Matrix.py", "Minors" ,"Moments.py", "Approximation.py", "Result.py"]

def run(program_list):
    for program in program_list:
        if program.endswith('.py'):
            subprocess.call(["python3", program])
        else:
            subprocess.call("./"+ program)        
        if program ==  "Result_exact.py":
                print("\nFinished: Result.py")
        else:
            print("\nFinished: " + program)

str_ = input("\nDo you need the exact calculation? (yes/no): ")

if str_.lower().startswith('y'): 
    run(program_list_1)
else:
    run(program_list_2)                           
