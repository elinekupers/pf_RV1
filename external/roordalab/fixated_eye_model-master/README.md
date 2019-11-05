Decoding (and simulating) foveal retinal ganglion cells. 

To use the code, go to the directory of this file (README.md) and run the following code:

ipython test/example_mnist.py -i

In interactive python mode, you can play around with the output of the algorithm by 
da.plot_EM_estimate(<iteration_number>)
where <iteration_number> is by default 0,1,...,9. 

Installation:

In addition to standard scientific computing packages, the code requires theano.

To run on GPU, make sure to add the .theanorc file in your home directory given here:

http://deeplearning.net/software/theano/library/config.html

I hope that you like the code!

-- Alex