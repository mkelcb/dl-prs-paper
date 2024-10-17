# -*- coding: utf-8 -*-

#MIT License

#Copyright (c) 2017 Marton Kelemen

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import matplotlib

import pandas as pd
import numpy as np

import platform

#data = results
# location = "../../../0cluster/results/local/test"
#acc_measure = "r^2"
def exportNNPlot(data,location, acc_measure = "r^2") :

    content = None
    cols = list()
    
    if platform.system().find('Windows') == -1 :
        print("Matlab uses Agg as we are on a *nix")
        matplotlib.use('Agg')
        
    import matplotlib.pyplot as plt
    
    if "train_accuracy" in data:
        #print("traing exists")
        cols.append('Traning')
        if(content is None) : content = data["train_accuracy"]
        else : content = np.column_stack( (content, data["train_accuracy"] ) )
        
    if "test_accuracy" in data:
        #print("test exists")
        cols.append('Validation')
        if(content is None) : content = data["test_accuracy"]
        else : content = np.column_stack( (content, data["test_accuracy"] ) ) 
        
        
    df = pd.DataFrame(content, index=data["epochs"], columns=cols )
    
  
    #df = df.cumsum()
    
    plt.figure()
    
    ax = df.plot(title = "Knet Accuracy") 
    ax.set_xlabel("epochs")
    ax.set_ylabel("Accuracy ("+acc_measure+")")
    

    fig = ax.get_figure()
    fig.savefig(location + '.eps', format='eps', dpi=1000)
    fig.savefig(location + '.png', dpi=300)
    
    
    
 
#content = np.random.randn(1000, 2)
#ind = range(1000)
#df = pd.DataFrame(content, index=ind, columns=['Traning', 'Validation'])
#df = df.cumsum()

#plt.figure()
#ax = df.plot() 

#fig = ax.get_figure()
#fig.savefig('../../../0cluster/results/local/asdf.eps', format='eps', dpi=1000)
#fig.savefig('../../../0cluster/results/local/asdf.png', dpi=300)
 
