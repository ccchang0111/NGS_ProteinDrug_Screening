# v3, 2016/10/31, take the clustered_stat file and submit to plotly
# v2: allow user to choose 1. percentage of peptide count or 2. ranking score for the heatmap
# Andrew Chang on 9/21/2016

# Purpose: visualize the clustered file using plotly
# note: need to register plotly account and install its package

import Tkinter
import tkFileDialog
import os
import pandas as pd
# import seaborn as sns
import plotly.plotly as py
import plotly.graph_objs as go

root3=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root3.withdraw()
root3.update()
file_path = tkFileDialog.askopenfilename(title='Please select the Clustered_Stats file for visualization')
root3.destroy()

file_name = os.path.basename(file_path).split('.')[0]

## import the test file that need to be compared ##
test = pd.read_excel(file_path, header=0)

##############################################
## Sum over the Rank Score for each cluster ##
##############################################
cluster_no = test['Cluster_ID'].unique()
idx_column = test.columns.str.contains('No.')

df = test.loc[:,idx_column].T

## Plotly for interactive heatmap
py.sign_in('ccchang0111', 'lqp8iy3ibx')
data = [go.Heatmap(
    x=df.columns.values.tolist(),
    y=df.index.values.tolist(),
    z=df.values.tolist(), 
    colorscale='spectral')]
py.iplot(data, filename='Cluster-Score_heatmap')

# https://plot.ly/61/~ccchang0111/#