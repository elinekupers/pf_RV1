from model import *
from analyzer import *

# Tests the algorithm using a variety of mnist digits and a sparse coding dictionary

N_diff_images = 10 # Run for <> different images
N_replicates = 10 # Run each image <> times
N_itr = 20

try:
    data = loadmat('data/mnist_dictionary.mat')
    D = data['D']
except IOError:
    print 'Need to have a dictionary file'
    raise IOError

_, N_pix = D.shape
L_I = int(np.sqrt(N_pix)) # Linear dimension of image

ig = ImageGenerator(L_I)

for _ in range(N_diff_images):
    ig.make_digit(mode = 'random')
    ig.normalize()
    S_gen = ig.img
    S_gen_name = ig.img_name
    emb = EMBurak(S_gen, D, N_T = 100, save_mode = True, S_gen_name = S_gen_name, N_itr = 20)
    for _ in range(N_replicates):
        emb.reset()
        emb.gen_data()
        emb.run_EM()
        filename = emb.save
        
