# ======== Select a file for opening:
import Tkinter,tkFileDialog
import subprocess

root = Tkinter.Tk()
filename1 = tkFileDialog.askopenfilename()
#file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose a file')
if filename1 != None:
    #data = file.read()
    #file.close()
    print "The name of the file is "+filename1
    
    
    
    
filename2 = tkFileDialog.askopenfilename()
#file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose a file')
if filename2 != None:
    #data = file.read()
    #file.close()
    print "The name of the file is "+filename2
    



subprocess.call("./RandomAdd "+filename1+" "+filename2+" 10 rop_p", shell=True)    
#subprocess.call("./RandomAdd %s %S", filename1, filename2) 

