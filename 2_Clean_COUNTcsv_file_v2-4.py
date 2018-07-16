# v2-4 Andrew Chang on 7/8/2018: allow user to assign specific length of the sequence
# v2-3 Andrew Chang on 12/5/2016: allow user to remove sequences with COUNT frequency less than X. 
# v2-2 Andrew Chang on 10/18/2016: for each column name (percentage and no. peptide..), added the filename at the end of the name. This will help clarify ambiguity for the next step
# v2-1 Andrew Chang on 4/15/2016: added new column that calculates the percentage of each seq. relative to the non-cleaned data. 
# v2 Andrew Chang on 4/11/2016

# Purpose: clean the .csv files in the COUNT folder; 
# 1. remove asterisk
# 2. remove frame shift 
# 3. add template info (frameworks of Reflexion or CKP) 
# 4. calculate the percentage of each sequence in this Extract file (note: it does't consider those short sequences removed in 1_Extract_Demultiplex_data)
# 5. remove those sequences with copy number less or equal to set number (default is 0, meaning no removal)

# Procedure:
# 1. User select all .csv files that need to be cleaned
# 2. User input whether to remove sequences with stop codon (but not the amber stop)
# 3. User input whether to remove sequences with frame shift. If yes, then provide the correct last two  amino acids sequences
# 4. User input whether to map each sequence to its template (only Reflexion and CKP are considered here, but the list can be expanded) 

import Tkinter
import tkFileDialog
import os
import pandas as pd

root6=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root6.withdraw()
root6.update()
file_paths = tkFileDialog.askopenfilenames(title='Please select the COUNT files that need to be cleaned (multiple files ok)')
root6.destroy()
rmv_stop = raw_input("Remove sequences contain stop codon (Y or N): ")
rmv_fshift = raw_input("Remove frame-shift (Y or N): ")
if (rmv_fshift == 'Y') or (rmv_fshift == 'y'):
    fshift = raw_input("Last two amino acids when there is no frame-shift (ex: SG for Relfx): ")
frm_search = raw_input("Map the framework (Y or N): ")
if (frm_search == 'Y') or (frm_search =='y'):
    lib = raw_input("Please specify the library: Relexion = 1, CKP = 2 ")
rmv_count  = raw_input("Remove sequences with copy number less or equal to X (integer, defaul X = 0): ")  

peptide_len_ext = raw_input("Extract specific peptide length (Y or N): ")
if (peptide_len_ext == 'Y') or (peptide_len_ext =='y'):
    peptide_len = raw_input("Please specify the length(s) of peptide (if multiple, separated each length by ','. ex: 10,11,12): ")
## clean the peptide_len
peptide_len = peptide_len.replace(" ","").split(",")
if "" in peptide_len:
    peptide_len.remove("")
peptide_len = map(int, peptide_len)

for file_path in file_paths:  
    file_name = os.path.basename(file_path).split('.')[0]
    test = pd.read_csv(file_path, header=0)
    test['frmw']='NA' # initialize the new column for storing framework info
    test_totcnt = test['No. peptide'].sum() # total number of sequences in the file before any cleaning
    test['Percentage'] = (test['No. peptide']/test_totcnt)*100 # add a new column of 'normalized counts' or 'percentage' based on the total seq. read in that sample (note: this percentate does not consider those short reads that were removed in 1_Extract)
    
    if (rmv_stop == 'Y') or (rmv_stop =='y'):
        idx=test['Peptide'].str.contains('\*')
        test = test.drop(test[idx].index)
    
    if (rmv_fshift == 'Y') or (rmv_fshift =='y'):
        test = test[test['Peptide'].str[-2:] == fshift]
    
    if (frm_search == 'Y') or (frm_search =='y'):
        Reflexion_dict={'2':'PPTE','3':'PQEL','4':'SEAL','7':'PDKL','8':'KVTV','23':'TVSE','24':'MMGR','27':'ISDT','28':'EACR','29':'NDDH','32':'TIDQ','37':'SSVK','42':'RPKT','47':'YKWE','53':'EGGE','56':'CSQN','75':'GAME','40':'GVEV','44':'QILF','64':'VQAM','65':'NCRS','66':'GHPL','70':'HDNY'}
        CKP_dict = {'EETI-II':'GCVC','Min23':'GSVC','AVR9':'CGRC','Circulin-A':'CSCK','Conotoxin-MVIIA':'CCTG','Huwentoxin':'CCPN','Charybdotoxin':'HNTS','CBD':'CASG','CBD_amylose':'CSGT','Mch1':'CDAG','Gurmarin':'CCEP','Asteropsin-A':'YPCC','AMP-1':'CCSQ','CPI':'CSGA','MCoTI-II':'SGSD','Kalata':'VCGE','Kalata_Thrombin':'MCGE'}
        
        if int(lib) == 1:
            for key, value in Reflexion_dict.iteritems():
                index = test['Peptide'].str.contains(value)
                test.loc[index,'frmw']=key
        
        if int(lib) == 2:
            for key, value in CKP_dict.iteritems():
                index = test['Peptide'].str.contains(value)
                test.loc[index,'frmw']=key
    
    minPercentage = test.Percentage.min() # remember the minimum percentage value, this will be used in the S/N calculation

    # double check if user has entered rmv_count or not. If not, do nothing and print out the following message.
    if rmv_count == '':
        print "No Copy Number Threshold is assigned, so all clean sequences are kepted."
    else: # if user does input something
        rmv_count = int(rmv_count) # convert string to int
        idx1 = test['No. peptide'] > rmv_count
        test = test[idx1]
    
    ## select specific peptide length
    if (peptide_len_ext == 'Y') or (peptide_len_ext =='y'):
        p_lengths = test.Peptide.str.len()
        idx2 = p_lengths.isin(peptide_len)
        test = test[idx2]

    test['L2'] = 'NA' # initialize new column for
    test['L2'] = test['Peptide'].str[-2:] # extract the last two amino acid (check if there is any frame-shift) into a new column
    test['min%'] = minPercentage
    test = test[['No. peptide','Percentage','Peptide','frmw','L2','min%']]
    test.columns = ['No. peptide_'+file_name,'Percentage_'+file_name, 'Peptide', 'frmw', 'L2', 'min%']
    test.to_csv(os.path.dirname(file_path) + "/Cleaned_" + file_name + ".csv", index=False) 