import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# Reads csv table that contains gene expression data
filename = input("Enter the name of the CSV file that contains regions, C9orf72, and the SOD1 data columns: ")

if os.path.exists(filename) and filename.endswith(".csv"):
    data = pd.read_csv(filename)

    if(data.empty):
        print("Does not have the needed data")
        sys.exit(1)

else:
    print("File does not exist in folder")
    sys.exit(1)


# Initializes counters
indexLength = 0
sodLength = 0
c9orLength = 0

expressionData = pd.DataFrame(data)

# gets the length of the different columns in the csv table
for i, j in expressionData.iterrows():
    if(pd.notna(j["Region"]) and j["Region"] != ""):
        if(pd.notna(j["SOD1"]) and j["SOD1"] != ""):
            if(pd.notna(j["C9orf72"]) and j["C9orf72"] != ""):
                indexLength += 3
                sodLength += 1
                c9orLength += 1


# initializes numpy array for each column in the csv file
loadDataValsGene = np.empty(indexLength, dtype=object)
sod1ValsAct = np.zeros(sodLength)
c9orf72ValsAct = np.zeros(c9orLength)

index = 0
sodIndex = 0
c9orIndex = 0

# translates data values in the csv file into numpy arrays
for k, m in expressionData.iterrows():
    if(pd.notna(m["Region"]) and m["Region"] != ""):
        if(pd.notna(m["SOD1"]) and m["SOD1"] != ""):
            if(pd.notna(m["C9orf72"]) and m["C9orf72"] != ""):
                loadDataValsGene[index] = m["Region"]
                index += 1
                loadDataValsGene[index] = m["SOD1"]
                index += 1
                loadDataValsGene[index] = m["C9orf72"]
                index += 1
                sod1ValsAct[sodIndex] = m["SOD1"]
                sodIndex += 1
                c9orf72ValsAct[c9orIndex] = m["C9orf72"]
                c9orIndex += 1

# turns the 1d array it initally is into 2d arrays
rowAmount = 0
if(len(loadDataValsGene) % 3 != 0):
    rowAmount = int(np.ceil(len(loadDataValsGene) / 3))
else:
    rowAmount = int(np.ceil(len(loadDataValsGene) / 3))
index = 0

array_2d = np.empty((rowAmount, 3), dtype= object)
for i in range(array_2d.shape[0]):
    for j in range(array_2d.shape[1]):
        if (index < len(loadDataValsGene)):
            array_2d[i][j] = loadDataValsGene[index]
            index += 1
        
print("done")

# creates array for the connectivity strength and the names of each region
valsEachRow = np.zeros(len(sod1ValsAct))
nameEachRow = np.zeros(len(sod1ValsAct))


#accesses the connectome csv file
filename = input("Enter the name of the CSV file that contains the regions(Connectome)): ")

if os.path.exists(filename) and filename.endswith(".csv"):
    data = pd.read_csv(filename) 

    if(data.empty):
        print("Does not have the needed data")
        sys.exit(1)

else:
    print("File does not exist")
    sys.exit(1)

vals = 0
sumVals = 0
index = 0

# loop that adds values to the connectivity strength of each column
# at the end get the total connectivity strength of all the data in the csv file

nameEachRow = data["Region"].values
for i, j in data.iterrows():
    vals = 0
    for k in data.columns:
        if(k != "Region"):
            vals += j[k]
        if(k == data.columns[-1]):
            valsEachRow[index] = vals
            sumVals += vals
            index += 1


actAtrophyVals = np.empty(len(sod1ValsAct), dtype= object)
actRegionVals = np.empty(len(sod1ValsAct), dtype= object)

#loads a csv of atrophy rates

filename = input("Enter the name of the CSV file that contains Atrophy Rates and Regions: ")

if os.path.exists(filename) and filename.endswith(".csv"):
    dataAtrophy = pd.read_csv(filename)

    if(dataAtrophy.empty):
        print("Does not have the needed data")
        sys.exit(1)
else:
    print("File does not exist")
    sys.exit(1)

regionIndex = 0
atrophyIndex = 0
for i, j in dataAtrophy.iterrows():
    if((pd.notna(j["Region"]) or j["Region"] != "") and pd.notna(j["AtrophyRate"]) or j["AtrophyRate"] != ""):
        actRegionVals[regionIndex] = j["Region"]
        regionIndex += 1
        actAtrophyVals[atrophyIndex] = j["AtrophyRate"]
        atrophyIndex += 1


# initializes a 2d array
rowAmount = 0
rowAmount = int(np.ceil(len(loadDataValsGene) / 2))

index = 0

array_2d = np.empty((rowAmount, 2), dtype= object)

# sorts it into a 2d array
for i in range(array_2d.shape[0]):
    for j in range(array_2d.shape[1]):
        if (index < len(loadDataValsGene)):
            array_2d[i][j] = loadDataValsGene[index]
            index += 1


# puts all the values into a oandas dataframe
allDataVals = ({'Region': nameEachRow, 'SOD1': sod1ValsAct, 'C9orf72': c9orf72ValsAct, 'Connectivity Strength': valsEachRow, 'Atrophy Rates': actAtrophyVals})

print(allDataVals)

valsDataFrame = pd.DataFrame(allDataVals)

# makes variables that get the values of correlation between each value 
corrConnectSOD1 = valsDataFrame["Connectivity Strength"].corr(valsDataFrame["SOD1"])
corrConnectC9or = valsDataFrame["Connectivity Strength"].corr(valsDataFrame["C9orf72"])

corrAtrophySOD1 = valsDataFrame["Atrophy Rates"].corr(valsDataFrame["SOD1"])
corrAtrophyC9or = valsDataFrame["Atrophy Rates"].corr(valsDataFrame["C9orf72"])


print(corrConnectC9or)
print(corrConnectSOD1)
print(corrAtrophyC9or)
print(corrAtrophySOD1)


# puts the correlation values into a dataframe
corrDataVals = ({"corrConnSOD1": corrConnectSOD1, "corrConnC9orf72": corrConnectC9or, "corrAtroC9orf72": corrAtrophyC9or, "corrAtroSOD1": corrAtrophySOD1})

corrDataFrame = pd.DataFrame([corrDataVals])



valsDataFrame.to_csv('scores.csv', index=False)

sod1ValsAct = sod1ValsAct.astype(float)
c9orf72ValsAct = c9orf72ValsAct.astype(float)
valsEachRow = valsEachRow.astype(float)
actAtrophyVals = actAtrophyVals.astype(float)

#creates a scatter plot for the correlation values and regression line

slopeSOD1Conn, intercept1 = np.polyfit(sod1ValsAct, valsEachRow, 1)
regressLine1 = slopeSOD1Conn * sod1ValsAct + intercept1

plt.scatter(sod1ValsAct, valsEachRow)
plt.plot(sod1ValsAct, regressLine1, linestyle='--', linewidth=1)
plt.title("Basic Scatter Plot")
plt.xlabel("SOD1 Values")
plt.ylabel("Connectivity Strength Values")
plt.show()



slopec9orfConn, intercept2 = np.polyfit(c9orf72ValsAct, valsEachRow, 1)
regressLine2 = slopec9orfConn * c9orf72ValsAct + intercept2

plt.scatter(c9orf72ValsAct, valsEachRow)
plt.plot(c9orf72ValsAct, regressLine2, linestyle='--', linewidth=1)
plt.title("Basic Scatter Plot")
plt.xlabel("C9orf72 Values")
plt.ylabel("Connectivity Strength Values")
plt.show()


slopeSOD1Atro, intercept3 = np.polyfit(sod1ValsAct, actAtrophyVals, 1)
regressLine3 = slopeSOD1Atro * sod1ValsAct + intercept3

plt.scatter(sod1ValsAct, actAtrophyVals)
plt.plot(sod1ValsAct, regressLine3, linestyle='--', linewidth=1)
plt.title("Basic Scatter Plot")
plt.xlabel("SOD1 Values")
plt.ylabel("Atrophy Values")
plt.show()


slopec9orfAtro, intercept4 = np.polyfit(c9orf72ValsAct, actAtrophyVals, 1)
regressLine4 = slopec9orfAtro * c9orf72ValsAct + intercept4

plt.scatter(c9orf72ValsAct, actAtrophyVals)
plt.plot(c9orf72ValsAct, regressLine4, linestyle='--', linewidth=1)
plt.title("Basic Scatter Plot")
plt.xlabel("C9orf72 Values")
plt.ylabel("Atrophy Values")
plt.show()





