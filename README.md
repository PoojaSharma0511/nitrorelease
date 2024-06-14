# nitrorelease

Step 1: Curation of Data
Training data files are created for nitrobenzene, o-nitroanisole, o-nitrophenol and m-nitrophenol with two coordinates C-N and C-O as input feature (in angstrom) and potential (in hartree) as output. (files are saved with name mol_train_data.out)


Step 2: Training the Machine learning model
Hyperparameter optimization is done while training to get optimized length value for the molecules using program 'training.py'.

Step 3: generation of kernel and alpha files for prediction
Optimized hyperparametrs are used in program kernel.py to get kernel_mol.inp and alpha_mol.inp files.

Step 4:
