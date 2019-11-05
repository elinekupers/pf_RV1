import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
from model import EMBurak
from analyzer import DataAnalyzer
from scipy.io import loadmat
from utils.image_gen import ImageGenerator

# Tests the algorithm using a '0' from mnist and a sparse coding dictionary
#  Also tests using the motion prior
try:
    data = loadmat('data/mnist_dictionary.mat')
    D = data['D']
except IOError:
    print 'Need to have a dictionary file'
    raise IOError

_, N_pix = D.shape
L_I = int(np.sqrt(N_pix))  # Linear dimension of image

ig = ImageGenerator(L_I)
# ig.make_digit(mode = 'random')
ig.make_digit()
ig.normalize()

S_gen = ig.img
S_gen_name = ig.img_name

emb = EMBurak(S_gen, D, n_t=30, save_mode=True,
              s_gen_name=S_gen_name, n_itr=3,
              path_mode='Experiment',
              motion_prior='VelocityDiffusion')

emb.gen_data()
emb.run_EM()

emb.save()

da = DataAnalyzer(emb.data)

# Plot the Estimated Image and Path after the algorithm ran
da.plot_EM_estimate(da.N_itr - 1)
plt.show()
