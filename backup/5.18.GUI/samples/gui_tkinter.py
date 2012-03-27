from sys import exit
from Tkinter import *
import subprocess
root = Tk() 
Button(root, text='Click to Exit !', command=subprocess.call('./a.out', shell=True)).pack() 
# Simple command

root.mainloop()
