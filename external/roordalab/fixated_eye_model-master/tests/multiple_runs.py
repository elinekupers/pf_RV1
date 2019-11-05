from model import *
from analyzer import *

# Tests the algorithm's ability to run multiple times without
# reinitializing the objects

N_replicates = 2  # Run each image <> times

try:
    data = loadmat('data/mnist_dictionary.mat')
    D = data['D']
except IOError:
    print 'Need to have a dictionary file'
    raise IOError

_, N_pix = D.shape
L_I = int(np.sqrt(N_pix))  # Linear dimension of image

ig = ImageGenerator(L_I)

ig.make_digit(mode='random')
ig.normalize()
S_gen = ig.img
emb = EMBurak(S_gen, D, N_T=20, save_mode=True)
filenames = []
for _ in range(N_replicates):
    emb.reset()
    emb.gen_data()
    emb.run_EM(N_itr=5)
    filename = emb.save
    filenames.append(filename)

da_s = []
for fn in filenames:
    da = DataAnalyzer.fromfilename(fn)
    da_s.append(da)
