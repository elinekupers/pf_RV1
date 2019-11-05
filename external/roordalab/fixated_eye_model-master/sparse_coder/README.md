Code to run sparse coding on MNIST

REQUIRES:

Create a .mat file with the MNIST digits (see first few lines of code for
how the data is loaded. Assumes that each image is a row of a matrix)

----

Using this code for other data:

Open sparse_coding_FISTA_theano.ipynb in ipython notebook. 

Choose the correct parameters in the first box after the package loading 
and run the rest of the script.

Notes: 

-- Put your dataset in a .mat file, to be loaded. See prep_mnist.py for
how to save in a matfile using scipy.io.savemat. The program expects an array
of dimension (N_examples, Image_width, Image_width)

-- Change the pos_only flag to False (usual application)

-- Set parameters: number of dictionary elements, sparse penalty strength (eg. lambda in (I-DA) + lambda |A|), batch size, etc

-- To print costs, set the flag show_costs=True when running func:`train`
