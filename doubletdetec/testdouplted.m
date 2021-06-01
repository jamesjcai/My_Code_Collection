% system('"d:\Miniconda3\envs\harmonypy\python.exe" script.py')

%%
pd = py.importlib.import_module('pandas');
np = py.importlib.import_module('numpy');
counts=pd.read_csv("Xsmall.csv").values;
dd = py.importlib.import_module('doubletdetection');

%%
return;

clf = dd.BoostClassifier(pyargs('n_iters',2,'use_phenograph',false,'standard_scaling',true));
% n_iters=2, use_phenograph=False, standard_scaling=True);
labels=clf.fit(counts.T).predict(pyargs('p_thresh',1e-16,'voter_thresh',0.5));
clf.doublet_score().mask

