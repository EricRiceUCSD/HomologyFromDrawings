from tkinter import *
from PIL import Image
import HomologyFunctions as funcs

# GUI Main Window
master = Tk()
master.title("Homology Calculator")

# Partitioning the GUI into left and right sides
leftFrame = Frame(master,
                  width=500,
                  height=500,
                  cursor="pencil")
leftFrame.pack(side=LEFT)
rightFrame = Frame(master,
                   width=250,
                   height=400)
rightFrame.pack(side=RIGHT)
rightFrame.pack_propagate(0)

# Initializing a canvas
drawing = Canvas(leftFrame,
                 width=500,
                 height=500)
drawing.pack(side=BOTTOM)

# Variable to track the left mouse button
left_but = "up"

# Variables to track the x & y when drawing
x_pos = None
y_pos = None

def left_but_down(event):
    global left_but
    left_but = "down"

def left_but_up(event):

    # Acquiring global variables
    global x_pos
    global y_pos
    global left_but
    
    left_but = "up"

    # Clearing position
    x_pos = None
    y_pos = None

# Drawing function which loops constantly while the mouse is moving
def draw(event):

    # Acquiring global variables
    global x_pos
    global y_pos
    global left_but
    
    if left_but == "down":
        
        # Checking for valid x & y, then drawing
        if x_pos is not None and y_pos is not None:
            drawing.create_line(x_pos, y_pos,\
                    event.x, event.y,\
                    smooth=TRUE, width=3.0)

        # Update position
        x_pos = event.x
        y_pos = event.y

# Binding our functions to the mouse buttons
drawing.bind("<Motion>", draw)
drawing.bind("<ButtonPress-1>", left_but_down)
drawing.bind("<ButtonRelease-1>", left_but_up)

# Adding labels
prompt = Label(leftFrame, text = "Press and drag the"+\
               " mouse to draw,"+'\n'+" then click Submit to"+\
               " calculate the homology of your drawing.")
prompt.pack(side=TOP)

outputLabel = Label(rightFrame, text="")
outputLabel.pack(side=BOTTOM)

bettiLabel = Label(rightFrame, text="")
bettiLabel.pack()
bettiNumbers = Canvas(rightFrame,
                      width=200,
                      height=200)
bettiNumbers.pack()

componentHoles = Label(leftFrame, text="")
componentHoles.pack()


# Defining functions that we will assign to buttons
def submit(event=None):
    global drawing
    global rightFrame
    global outputLabel
    global bettiLabel
    global bettiNumbers
    global componentHoles
    global blockSize
    global radius

    print("Calculating...")
    drawing.postscript(file="CanvasPostscript.eps", colormode="monochrome")
    img = Image.open("CanvasPostscript.eps")
    img.save("CanvasPNG.png", "png")

    K = funcs.CechComplex(
            funcs.getReducedVertexSet(
                funcs.getValsByCoord("CanvasPNG.png"), int(blockSize.get())), float(radius.get()))

    if K.betti()[1] == 1:
        outputLabel.config(text="Your drawing has "+str(K.betti()[1])+" hole.")
    else:
        outputLabel.config(text="Your drawing has "+str(K.betti()[1])+" holes.")

    bettiString = ""
    for p in range(0, min(6, len(K.betti()))):
        bettiString += "β_"+str(p)+" = "+str(K.betti()[p])+"\n"
    bettiLabel.config(font=("TkDefaultFont", 12),
                      text="The Betti numbers are:")
    bettiNumbers.delete("all")
    bettiNumbers.create_text(100, 0, anchor=N,
                             font=("Times New Roman", 20, "bold"),
                             text=bettiString)

    if K.betti()[0] > 1:
        C = funcs.splitByComponent(K)
        componentsString = "From left to right, the components of your\n"+\
                "drawing have "
        if len(C) == 2:
            componentsString += str(C[0].betti()[1])+" "
        else:
            for j in range(0, len(C) - 1):
                componentsString += str(C[j].betti()[1])+", "
        componentsString += "and "+str(C[-1].betti()[1])+" holes, respectively."
        componentHoles.config(text=componentsString)

    print("Done.")

def clear(event=None):
    global drawing
    global outputLabel
    global closestMatchLabel
    global closestMatchDisplay
    global componentHoles
    
    drawing.delete("all")
    outputLabel.config(text="")
    bettiLabel.config(text="")
    bettiNumbers.delete("all")
    componentHoles.config(text="")

# Creating buttons and variable entries
submit_but = Button(rightFrame, text="Submit", command=submit)
submit_but.pack(ipadx=50)
master.bind("<Return>", submit)
clear_but = Button(rightFrame, text="Clear", command=clear)
clear_but.pack(ipadx=55)

blockSizeLabel = Label(rightFrame, text="Block Size")
blockSizeLabel.pack()
blockSize = Entry(rightFrame)
blockSize.pack()
blockSize.insert(0, "20")

radiusLabel = Label(rightFrame, text="Radius")
radiusLabel.pack()
radius = Entry(rightFrame)
radius.pack()
radius.insert(0, "0.8")

def __init__():   
    mainloop()
