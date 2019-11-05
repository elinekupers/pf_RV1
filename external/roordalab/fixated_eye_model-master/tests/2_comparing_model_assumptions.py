"""
The goal of this script is to run the model using different assumptions
and to then compare the results
"""

import numpy as np
from scipy.io import loadmat

import sys
sys.path.append('..')

from src.model import EMBurak
from utils.image_gen import ImageGenerator

output_dir = 'comparing_models'

L_I = 14

ig = ImageGenerator(L_I)
ig.make_big_E()
ig.normalize()
DC = 100.
n_t = 100
D = np.eye((L_I ** 2))
LAMBDA = 0.0

# Motion Prior: no motion, Image Prior: Independent Pixels

emb = EMBurak(ig.img, D, n_t=n_t, LAMBDA=LAMBDA, save_mode=True,
              s_gen_name=ig.img_name, dc_gen=DC, dc_infer=0.01,
              output_dir=output_dir)
emb.gen_data()
emb.run_EM()
emb.save()

# Motion Prior: diffusion, Image Prior: Independent Pixels

emb = EMBurak(ig.img, D, n_t=n_t, LAMBDA=LAMBDA, save_mode=True,
              s_gen_name=ig.img_name, dc_gen=DC, dc_infer=DC,
              output_dir=output_dir)
emb.gen_data()
emb.run_EM()
emb.save()

# Motion Prior: diffusion, Image Prior: MNIST

data = loadmat('data/mnist_dictionary.mat')
D = data['D']

_, N_pix = D.shape
L_I = int(np.sqrt(N_pix))  # Linear dimension of image

ig = ImageGenerator(L_I)
ig.make_digit()
ig.normalize()


emb = EMBurak(ig.img, D, n_t=n_t, save_mode=True,
              s_gen_name=ig.img_name, dc_gen=DC, dc_infer=DC,
              output_dir=output_dir)
emb.gen_data()
emb.run_EM()
emb.save()
