from Tkinter import *
from tkMessageBox import *
from tkColorChooser import askcolor              
from tkFileDialog   import askopenfilename      

def callback():
    askopenfilename() 
    
    
errmsg = 'Error!'
Button(text='Open', command=callback).pack(fill=X)
Button(text='Spam', command=(lambda: showerror('Spam', errmsg))).pack(fill=X)
mainloop()
