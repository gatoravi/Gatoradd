# ======== Select a file for opening:
import Tkinter,tkFileDialog
import subprocess
from Tkinter import *

#filename1 = tkFileDialog.askopenfilename()
#file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose a file')
#if filename1 != None:
    #data = file.read()
    #file.close()
    #print "The name of the file is "+filename1

#filename2 = tkFileDialog.askopenfilename()
#file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose a file')
#if filename2 != None:
    #data = file.read()
    #file.close()
    #print "The name of the file is "+filename2



class App:
  def __init__(self, master):

    frame = Frame(master)
    frame.pack()

    self.w = Label(root, text="PHYADD beta ver!!")
    self.w.pack()

    self.button = Button(frame, text="QUIT", fg="red", command=frame.quit)
    self.button.pack(side=LEFT)
  
    self.button = Button(frame, text="Taxa File", fg="red", command=self.select_taxa_file)
    self.button.pack(side=RIGHT) 
 
    self.button = Button(frame, text="Tree File", fg="red", command=self.select_tree_file)
    self.button.pack(side=RIGHT)
    
    self.hi_there = Button(frame, text="Run <----!!", command=self.run)
    self.hi_there.pack(side=LEFT)
    
    self.entry = Entry(root, width=10).pack(side=TOP,padx=10,pady=10)

  def say_hi(self):
    print "hi there, everyone!"
  
  def select_tree_file(self):
    self.tree_filename = tkFileDialog.askopenfilename()
    print self.tree_filename
    
  def select_taxa_file(self):
    self.taxa_filename = tkFileDialog.askopenfilename()
    print self.taxa_filename
    
  def run(self):
    subprocess.call("./RandomAdd "+self.tree_filename+" "+self.taxa_filename+" 10 rop_pgui", shell=True)

root = Tk()

app = App(root)

root.mainloop()


#subprocess.call("./RandomAdd "+filename1+" "+filename2+" 10 rop_p", shell=True)    
#subprocess.call("./RandomAdd %s %S", filename1, filename2) 

