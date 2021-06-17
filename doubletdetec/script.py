import os
import pandas as pd
import numpy as np
import doubletdetection as dd
os.chdir("U:\\GitHub\\My_Code_Collection\\doubletdetec")

counts=pd.read_csv("Xsmall.csv").values

#clf = dd.BoostClassifier()
#labels = clf.fit(raw_counts).predict()
#scores = clf.doublet_score()


#counts = np.random.poisson(size=(500, 100))

    # no phenograph
clf = dd.BoostClassifier(n_iters=2, use_phenograph=False, standard_scaling=True)
labels=clf.fit(counts.T).predict(p_thresh=1e-16, voter_thresh=0.5)
pd.DataFrame(clf.doublet_score().mask).to_csv('output.csv',index=False,header=False)


# pd.DataFrame(clf.doublet_score().mask).to_csv('output.csv',index=False,header=['isdoublet'])

# with phenograph
#    clf = doubletdetection.BoostClassifier(n_iters=2, use_phenograph=True, standard_scaling=True)
#    clf.fit(counts).predict(p_thresh=1e-16, voter_thresh=0.5)
#    clf.doublet_score()
#    doubletdetection.plot.convergence(clf, show=False, p_thresh=1e-16, voter_thresh=0.5)
#    doubletdetection.plot.threshold(clf, show=False, p_step=6)