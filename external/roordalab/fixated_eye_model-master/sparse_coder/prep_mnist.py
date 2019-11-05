# Script takes in data from mnist mat file and creates a pkl file
# with all of the training set in one matrix
# labels has the labels of the digits

import numpy as np
import cPickle as pkl
from scipy.io import loadmat
from skimage.measure import block_reduce
from scipy.io import savemat


data_dir = "data/final/"

data = loadmat('data/mnist_all.mat')

labels = []
labels.append(0)
for i in range(10):
    labels.append(data['train' + str(i)].shape[0])
labels = np.array(labels)

idx = np.cumsum(labels)
IMAGES = np.zeros((idx[-1], 784), dtype = 'int')
LABELS = np.zeros((idx[-1],), dtype = 'int')

for i in range(10):
    IMAGES[idx[i]:idx[i + 1], :] = data['train' + str(i)]
    LABELS[idx[i]:idx[i + 1]] = i
 
IMAGES = IMAGES.reshape(60000, 28, 28)

IMAGES1 = block_reduce(IMAGES.astype('float32'), 
                       block_size = (1, 2, 2), 
                       func = np.mean)   

savemat(data_dir + 'mnist.mat', {'IMAGES': IMAGES1.astype('int16'),
                      'LABELS': LABELS.astype('int')})

#pkl.dump(LABELS, open(data_dir + 'mnist_labels.pkl', 'wb'))

#pkl.dump(IMAGES1.astype('int16'), open(data_dir + 'mnist.pkl', 'wb'))
