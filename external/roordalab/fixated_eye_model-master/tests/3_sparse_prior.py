"""
Code to compare different prior strengths
"""

import numpy as np
from scipy.io import loadmat

from src.model import EMBurak
from utils.image_gen import ImageGenerator

output_dir = 'sparsity'
DC = 50.
n_t = 100
n_repeats = 5

# Sparse coding dictionary prior
data = loadmat('sparse_coder/output/mnist_dictionary.mat')
D = data['D']

_, N_pix = D.shape
L_I = int(np.sqrt(N_pix))  # Linear dimension of image

ig = ImageGenerator(L_I)
ig.make_digit()
ig.normalize()


emb = EMBurak(ig.img, D, n_t=100, save_mode=True,
              s_gen_name=ig.img_name, dc_gen=DC, dc_infer=DC,
              output_dir=output_dir, LAMBDA=0.)

for _ in range(n_repeats):
    emb.reset()
    emb.gen_data()
    emb.run_EM()
    emb.save()

# Independent Pixel Prior
D = np.eye((L_I ** 2))
emb = EMBurak(ig.img, D, n_t=n_t, LAMBDA=0.0, save_mode=True,
              s_gen_name=ig.img_name, dc_gen=DC, dc_infer=DC,
              output_dir=output_dir)

for _ in range(n_repeats):
    emb.reset()
    emb.gen_data()
    emb.run_EM()
    emb.save()




# convert -set delay 30 -colorspace GRAY -colors 256 -dispose 1 -loop 0 -scale 50% *.png alg_performance.gif

# convert -set delay 30 -colors 256 -dispose 1 -loop 0 *.jpg alg_performance.gif
