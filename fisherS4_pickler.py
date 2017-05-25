import pickle
import numpy as np

fisherS4=np.array(np.loadtxt("fisherS4.txt"))

pickle.dump(fisherS4, open('fisherS4_matrix.p', 'wb'))
