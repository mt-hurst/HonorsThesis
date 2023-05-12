# This python code was run from a Google Colaboratory Notebook with access to my 
# csv files for training. If you want to run it locally, you'll need
# sklearn >= 1.2.1
# 

import numpy as np
import matplotlib.pyplot as plt
#Colaboratory was loading the wrong sklearn earlier, so I was manually setting it to version 1.2.1.
#sklearn version 1.2.2 is now out, and I don't need the !pip install command anymore
#!pip install scikit-learn==1.2.1
import sklearn
from sklearn import inspection

from sklearn.ensemble import RandomForestClassifier
from sklearn import tree

from sklearn import metrics

#Importing the csv files from my Google Drive
#To get this to work, you need to have the python
from google.colab import files
f = open('/content/drive/MyDrive/Spring 23/Thesis ML/CSVs/pionrate.csv')
g = open('/content/drive/MyDrive/Spring 23/Thesis ML/CSVs/mollrate.csv')

pionlines = f.readlines()
molllines = g.readlines()
combined = pionlines + molllines

#Function for separating data into training and testing set
def train_dat_create(datalines):
  '''Takes in CVS data in a list and returns two lists of Training Data and Testing Data.
      Training Data is 80% of the input data randomly selected, and Testing Data is 20% of the input data randomly selected'''
  #Initializing lists
  test = []
  train = []
  for datline in datalines:
    rand = np.random.rand()#Create a random float between 0 and 1
    if int(datline.rsplit(',')[0]) == 0:#If the event is a Moll Generated Event
      if rand < .8:
        train.append((int(datline.rsplit(',')[0]),int(datline.rsplit(',')[1]),int(datline.rsplit(',')[2]),int(datline.rsplit(',')[3]),float(datline.rsplit(',')[4])))
      else: test.append((int(datline.rsplit(',')[0]),int(datline.rsplit(',')[1]),int(datline.rsplit(',')[2]),int(datline.rsplit(',')[3]),float(datline.rsplit(',')[4])))
    if int(datline.rsplit(',')[0]) == 1:#If the event is a Pion Generated event
      if rand <.8:#Put 80% of the pion events into the training set
        train.append((int(datline.rsplit(',')[0]),int(datline.rsplit(',')[1]),int(datline.rsplit(',')[2]),int(datline.rsplit(',')[3]),float(datline.rsplit(',')[4])))
      else: test.append((int(datline.rsplit(',')[0]),int(datline.rsplit(',')[1]),int(datline.rsplit(',')[2]),int(datline.rsplit(',')[3]),float(datline.rsplit(',')[4])))#Put the other 20% of pion events in the teaching set
  return test,train



testset, trainset = train_dat_create(combined)
molltestset, molltrainset = train_dat_create(molllines)
piontestset, piontrainset = train_dat_create(pionlines)

pionarr = np.array(piontrainset)
mollarr = np.array(molltrainset)

testarr = np.array(testset) # This is the array that I use for testing
trainarr = np.array(trainset) # This is the array that I use for Training

traininput = trainarr[0:,1:3]
trainassignments = trainarr[0:,0]
trainweights = trainarr[0:,4]


# testing impact of number of estimators on results
estimators = 1
est_x = [] #x coord list
y = [] # no pion exit classifier performace
yy = [] # with pion exit classifier performance
for i in range(1,21):

  #Creates Binary Classifier where n_estimators is the number of Decision Trees
  clf_simple = RandomForestClassifier(n_estimators = estimators, random_state = 1111);
  #Fits the classifier with only SHMX Detector Signal and the Pion Detector signal
  clf_simple = clf_simple.fit(traininput,np.ravel(trainassignments));

  #Creates Binary Classifier where n_estimators is the number of Decision Trees
  clf_trig = RandomForestClassifier(n_estimators = estimators, random_state = 1111);
  #Fits the classifier with the SHMX Detector Signal, Pion Detector Signal, 
  # AND the boolean of whether the Pion Exit Scintillator saw a charged particle
  clf_trig = clf_trig.fit(trainarr[0:,1:4], np.ravel(trainassignments));
  est_x.append(estimators)
  estimators = 3*i
  y.append(clf_simple.score(testarr[0:,1:3],testarr[0:,0]))
  yy.append(clf_trig.score(testarr[0:,1:4],testarr[0:,0]))


plt.plot(est_x,y,'--ob',label = 'No Exit Scint')
plt.plot(est_x,yy,'-or',label = 'With Exit Scint')
plt.ylabel("Percent of events Correctly Classified", size = 15)
plt.xlabel("Number of Estimator Decision Trees", size = 15)
plt.legend()


#Create Classifiers

#Classifier with NO rate weighting

clf_old = RandomForestClassifier(n_estimators = 20, random_state = 1111)
clf_old = clf_old.fit(traininput, np.ravel(trainassignments))

#Creates Binary Classifier where n_estimators is the number of Decision Trees
clf = RandomForestClassifier(n_estimators = 20, random_state = 1111) 
#Fits the classifier with only SHMX Detector Signal and the Pion Detector signal
clf = clf.fit(traininput,np.ravel(trainassignments),sample_weight = trainweights)

print('The accuracy of the Random Forest Classifier without the Pion Exit data is WITHOUT RATING : '+ str(clf_old.score(testarr[0:,1:3],testarr[0:,0])))
print('The accuracy of the Random Forest Classifier without the Pion Exit data is : '+ str(clf.score(testarr[0:,1:3],testarr[0:,0])))

#Plotting the boundaries of the two Classifiers

d = inspection.DecisionBoundaryDisplay.from_estimator(clf_old,X =(testarr[0:,1:3]),grid_resolution = 100,cmap = plt.cm.RdBu,alpha = 0.8)
#display.plot()
plt.ylim(0,1000)
plt.xlim(0,500)
pion_1 = []
moll_1 = []
for ent in trainarr[0:]:
  if ent[0] == 0:
    moll_1.append((ent[1],ent[2]))
  if ent[0] == 1:
    pion_1.append((ent[1],ent[2]))
  
pion_1_arr = np.array(pion_1)
moll_1_arr = np.array(moll_1)

plt.scatter(moll_1_arr[0:,0],moll_1_arr[0:,1],s=2,c='r',marker = 'x', label = 'Moller Event')
plt.scatter(pion_1_arr[0:,0],pion_1_arr[0:,1],s=2,c='b',marker = '^', label = 'Pion Event')
plt.xlabel("SHMX Data",size = 14)
plt.ylabel("PionDetector Data",size = 14)
plt.title("Boundary of Classifier Regions")
plt.legend()
plt.colorbar(d.ax_.collections[1])
d.ax_.collections[1].set_cmap(plt.cm.RdBu)

h = inspection.DecisionBoundaryDisplay.from_estimator(clf,X =(testarr[0:,1:3]),grid_resolution = 100,cmap = plt.cm.RdBu,alpha = 0.8)
#display.plot()
plt.ylim(0,1000)
plt.xlim(0,500)
pion_2 = []
moll_2 = []
for ent in trainarr[0:]:
  if ent[0] == 0:
    moll_2.append((ent[1],ent[2]))
  if ent[0] == 1:
    pion_2.append((ent[1],ent[2]))
pion_2_arr = np.array(pion_2)
moll_2_arr = np.array(moll_2)

plt.scatter(moll_2_arr[0:,0],moll_2_arr[0:,1],s=2,c='r',marker = 'x', label = 'Moller Event')
plt.scatter(pion_2_arr[0:,0],pion_2_arr[0:,1],s=2,c='b',marker = '^', label = 'Pion Event')
plt.xlabel("SHMX Data",size = 14)
plt.ylabel("PionDetector Data",size = 14)
plt.title("Boundary of Classifier Regions with Rate Weighting")
plt.legend()
plt.colorbar(h.ax_.collections[1])
h.ax_.collections[1].set_cmap(plt.cm.RdBu)

#Replace the clf_old with whatever classifier you want to see the confusion matrix of
cm = metrics.confusion_matrix(y_labels, clf_old.predict(testarr[0:,1:3]))
cm_display = metrics.ConfusionMatrixDisplay(cm).plot()
plt.title("Confusion Matrix", size = 16)
