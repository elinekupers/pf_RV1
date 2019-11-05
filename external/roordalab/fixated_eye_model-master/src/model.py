# Python script containing a class that does expectation maximization
#   to estimate an image from simulated LGN cell responses
# See the end of the script for a sample usage

import numpy as np
import theano
import theano.tensor as T
import os
from utils.path_generator import (DiffusionPathGenerator,
                                  ExperimentalPathGenerator)
import utils.particle_filter_new as pf
# from utils.theano_gradient_routines import ada_delta
from utils.fista import fista_updates
from utils.SNR import SNR
import cPickle as pkl
from utils.time_filename import time_string
from utils.BurakPoissonLP import PoissonLP


class TheanoBackend(object):
    """
    Theano backend for executing the computations
    """
    def __init__(self, XS, YS, XE, YE, IE, Var, A, d, pos_only=True):
        """
        Initializes all theano variables and functions
        """
        self.XS = XS
        self.YS = YS
        self.XE = XE
        self.YE = YE
        self.IE = IE
        self.Var = Var
        self.A = A
        self.d = d
        self.l_i = XS.shape[0]

        # Define Theano Variables Common to Generation and Inference
        self.t_XS = theano.shared(self.XS, 'XS')
        self.t_YS = theano.shared(self.YS, 'YS')
        self.t_XE = theano.shared(self.XE, 'XE')
        self.t_YE = theano.shared(self.YE, 'YE')
        self.t_IE = theano.shared(self.IE, 'IE')
        self.t_Var = theano.shared(self.Var, 'Var')

        self.t_XR = T.matrix('XR')
        self.t_YR = T.matrix('YR')
        self.t_R = T.matrix('R')

        #  Parameters
        self.t_L0 = T.scalar('L0')
        self.t_L1 = T.scalar('L1')
        self.t_DT = T.scalar('DT')
        self.t_DC = T.scalar('DC')
        self.t_G = T.scalar('G')

        def inner_products(t_S, t_Var, t_XS, t_YS, t_XE, t_YE, t_XR, t_YR):
            """
            Take dot product of image with an array of gaussians
            t_S - image variable shape - (i2, i1)
            t_Var - variances of receptive fields
            t_XS - X coordinate for image pixels for dimension i1
            t_YS - Y coordinate for image pixels for dimension i2
            t_XE - X coordinate for receptive fields j
            t_YE - Y coordinate for receptive fields j
            t_XR - X coordinate for retina in form batch, timestep, b,t
            t_YR - Y coordinate ''
            """

            # Note in this computation, we do the indices in this form:
            #  b, i, j, t
            #  batch, pixel, neuron, time step

            # indices: b, i1, j, t
            t_dX = (t_XS.dimshuffle('x', 0, 'x', 'x') -
                    t_XE.dimshuffle('x', 'x', 0, 'x') -
                    t_XR.dimshuffle(0, 'x', 'x', 1))
            t_dX.name = 'dX'
            # indices: b, i2, j, t
            t_dY = (t_YS.dimshuffle('x', 0, 'x', 'x') -
                    t_YE.dimshuffle('x', 'x', 0, 'x') -
                    t_YR.dimshuffle(0, 'x', 'x', 1))
            t_dY.name = 'dY'

            # Use outer product trick to dot product image with point filters
            t_PixRFCouplingX = T.exp(-0.5 * t_dX ** 2 /
                                     t_Var.dimshuffle('x', 0, 'x', 'x'))
            t_PixRFCouplingY = T.exp(-0.5 * t_dY ** 2 /
                                     t_Var.dimshuffle('x', 0, 'x', 'x'))
            t_PixRFCouplingX.name = 'PixRFCouplingX'
            t_PixRFCouplingY.name = 'PixRFCouplingY'

            # Matrix of inner products between the images and the retinal RFs
            # indices: b, j, t
            # Sum_i2 T(i2, i1) * T(b, i2, j, t) = T(b, i1, j, t)
            t_IpsY = T.sum(t_S.dimshuffle('x', 0, 1, 'x', 'x') *
                           t_PixRFCouplingY.dimshuffle(0, 1, 'x', 2, 3),
                           axis=1)
            # Sum_i1 T(b, i1, j, t) * T(b, i2, j, t) = T(b, j, t)
            t_Ips = T.sum(t_IpsY * t_PixRFCouplingX, axis=1)
            t_Ips.name = 'Ips'

            # For the gradient, we also prepare d Ips / dS
            # This is in the form b, i2, i1, j, t
            t_PixRFCoupling = (t_PixRFCouplingX.dimshuffle(0, 'x', 1, 2, 3) *
                               t_PixRFCouplingY.dimshuffle(0, 1, 'x', 2, 3))

            return t_Ips, t_PixRFCoupling

        # self.inner_products = inner_products

        def firing_prob(t_Ips, t_G, t_IE, t_L0, t_L1, t_DT):
            # Firing probabilities indexed by b, j, t
            # t_Ips - Image-RF inner products indexed as b, j, t
            # t_G - gain constant
            # t_IE - identity of retinal ganglion cells
            # t_L0, t_L1 - min, max firing rate
            # t_DT - time step size

            t_IEr = t_IE.dimshuffle('x', 0, 'x')
            t_Gen = t_IEr + (1 - 2 * t_IEr) * t_G * t_Ips  # Generator signal

            t_FP_0 = t_DT * T.exp(T.log(t_L0) + T.log(t_L1 / t_L0) * t_Gen)

            t_FP = T.switch(t_FP_0 > 0.9, 0.9, t_FP_0)
            return t_FP

        self.firing_prob = firing_prob

        # Simulated Spike Generation

        self.t_S_gen = T.matrix('S_gen')  # Image dims are i2, i1
        self.t_Ips_gen, _ = inner_products(self.t_S_gen, self.t_Var,
                                           self.t_XS, self.t_YS,
                                           self.t_XE, self.t_YE,
                                           self.t_XR, self.t_YR)
        self.t_FP_gen = firing_prob(self.t_Ips_gen, self.t_G, self.t_IE,
                                    self.t_L0, self.t_L1, self.t_DT)

        # Computes image-RF inner products and the resulting firing
        # probabilities
        self.RFS = theano.function(inputs=[self.t_S_gen, self.t_XR, self.t_YR,
                                           self.t_L0, self.t_L1,
                                           self.t_DT, self.t_G],
                                   outputs=[self.t_Ips_gen, self.t_FP_gen])
        self.rng = T.shared_randomstreams.RandomStreams(seed=10)
        self.t_R_gen = (self.rng.uniform(size=self.t_FP_gen.shape) <
                        self.t_FP_gen).astype('float32')

        self.spikes = theano.function(inputs=[self.t_S_gen,
                                              self.t_XR, self.t_YR,
                                              self.t_L0, self.t_L1,
                                              self.t_DT, self.t_G],
                                      outputs=self.t_R_gen)

        # Latent Variable Estimation

        def spiking_cost(t_R, t_FP):
            """
            Returns the negative log likelihood of the spikes given
            the inner products
            t_R - spikes in form j, t
            t_FP - Firing probabilities in form b, j, t
            Returns -log p(R|X,S) with indices in the form b, j, t
            """
            #            t_E_R_f = -(t_R.dimshuffle('x', 0, 1) * T.log(t_FP)
            #         + (1 - t_R.dimshuffle('x', 0, 1)) * T.log(1 - t_FP))
            #         Try using poisson loss instead of bernoulli loss
            t_E_R_f = -t_R.dimshuffle('x', 0, 1) * T.log(t_FP) + t_FP

            t_E_R_f.name = 'E_R_f'

            return t_E_R_f

        self.spiking_cost = spiking_cost

        self.t_A = theano.shared(self.A, 'A')
        self.t_D = theano.shared(self.d, 'D')
        self.t_S = T.dot(self.t_A, self.t_D).reshape((self.l_i, self.l_i))
        self.image_est = theano.function(inputs=[], outputs=self.t_S)
        # FIXME: shouldn't access L_I, Image dims are i2, i1

        self.t_GAMMA = T.scalar('GAMMA')
        self.t_LAMBDA = T.scalar('LAMBDA')

        self.t_Ips, _ = inner_products(self.t_S, self.t_Var,
                                       self.t_XS, self.t_YS,
                                       self.t_XE, self.t_YE,
                                       self.t_XR, self.t_YR)

        self.t_FP = firing_prob(self.t_Ips, self.t_G, self.t_IE,
                                self.t_L0, self.t_L1, self.t_DT)

        # Compute Energy Functions (negative log-likelihood) to minimize
        # Weights (batch, timestep) from particle filter
        self.t_Wbt = T.matrix('Wbt')
        self.t_E_R_f = spiking_cost(self.t_R, self.t_FP)

        self.t_E_R = T.sum(T.sum(self.t_E_R_f, axis=1) * self.t_Wbt)
        self.t_E_R.name = 'E_R'

        self.t_E_bound = self.t_GAMMA * (
            T.sum(T.switch(self.t_S < 0., -self.t_S, 0)) +
            T.sum(T.switch(self.t_S > 1., self.t_S - 1, 0)))
        self.t_E_bound.name = 'E_bound'

        self.t_E_sp = self.t_LAMBDA * T.sum(T.abs_(self.t_A))
        self.t_E_sp.name = 'E_sp'

        self.t_E_rec = self.t_E_R + self.t_E_bound
        self.t_E_rec.name = 'E_rec'

        self.t_E = self.t_E_rec + self.t_E_sp
        self.t_E.name = 'E'

        # Cost from poisson terms separated by batches for particle filter log
        # probability
        self.t_E_R_b = T.sum(self.t_E_R_f, axis=(1, 2))
        self.spike_energy = theano.function(inputs=[self.t_XR, self.t_YR,
                                                    self.t_R,
                                                    self.t_L0, self.t_L1,
                                                    self.t_DT, self.t_G],
                                            outputs=self.t_E_R_b)

        # Generate costs given a path, spikes, and time-batch weights
        energy_outputs = [self.t_E, self.t_E_bound, self.t_E_R, self.t_E_sp]
        self.costs = theano.function(
            inputs=[
                self.t_XR, self.t_YR, self.t_R, self.t_Wbt,
                self.t_L0, self.t_L1, self.t_DT, self.t_G,
                self.t_GAMMA, self.t_LAMBDA],
            outputs=energy_outputs)

        # Define theano variables for gradient descent
        # self.t_Rho = T.scalar('Rho')
        # self.t_Eps = T.scalar('Eps')
        # self.t_ada_params = (self.t_Rho, self.t_Eps)

        # self.grad_updates = ada_delta(self.t_E, self.t_A, *self.t_ada_params)
        # self.t_A_Eg2, self.t_A_EdS2, _ = self.grad_updates.keys()

        # Define variables for FISTA minimization
        self.t_L = T.scalar('L')

        self.grad_updates = fista_updates(
            self.t_A, self.t_E_rec, self.t_LAMBDA,
            self.t_L, pos_only=pos_only)

        _, self.t_fista_X, self.t_T = self.grad_updates.keys()

        # Initialize t_A, and extra variables

        inputs = [self.t_XR, self.t_YR, self.t_R, self.t_Wbt,
                  self.t_L0, self.t_L1, self.t_DT, self.t_G,
                  # self.t_GAMMA, self.t_LAMBDA, self.t_Rho, self.t_Eps]
                  self.t_GAMMA, self.t_LAMBDA, self.t_L]
        self.img_grad = theano.function(inputs=inputs,
                                        outputs=energy_outputs,
                                        updates=self.grad_updates)

    def get_A(self):
        return self.t_A.get_value()

    def reset_image_estimate(self):
        """
        Resets the value of the image as stored on the GPU
        """
        self.t_A.set_value(np.zeros_like(self.A).astype('float32'))

    def calculate_L(self, n_t, n_n, l0, l1, dt, d_scl, ctant):
        """
        Return the value of the Lipschitz constant of the smooth part
            of the cost function

        Parameters
        ----------
        n_t : int
            number of timesteps
        n_n : int
            number of neurons
        l0 : float
            baseline firing rate
        l1 : float
            maximum firing rate
        dt : float
            timestep
        d_scl : float
            scale of dictionary elements (mean sum of squares)
        ctant : float
            constant to loosen our bound
        """
        return (n_t * n_n * l1 * dt * d_scl ** 2 *
                np.log(l1/l0) ** 2 * ctant).astype('float32')

    def reset_m_aux(self):
        """
        Resets auxillary gradient descent variables for the M step
            eg. fista we reset the copy of A and the step size
        """
        # self.t_A_Eg2.set_value(np.zeros_like(self.A).astype('float32'))
        # self.t_A_EdS2.set_value(np.zeros_like(self.A).astype('float32'))
        self.t_T.set_value(np.array([1.]).astype(theano.config.floatX))
        # self.t_fista_X.set_value(np.zeros_like(self.A).astype('float32'))
        self.t_fista_X.set_value(self.t_A.get_value())


class EMBurak(object):

    def __init__(self, s_gen, d, dt=0.001,
                 n_t=50,
                 l_n=14, a=1., LAMBDA=0.,
                 save_mode=False,
                 n_itr=20, s_gen_name=' ',
                 motion_gen_mode='Diffusion', dc_gen=100.,
                 motion_prior='PositionDiffusion', dc_infer=100.,
                 output_dir=''):
        """
        Initializes the parts of the EM algorithm
            -- Sets all parameters
            -- Initializes Dictionary
            -- Compiles the theano backend
            -- Sets the gain factor for the spikes
            -- Initializes the object that generates the paths
            -- Initializes the Particle Filter object
            -- Checks that the output directory exists

        Parameters
        ----------
        s_gen : array, float32, shape (l_i, l_i)
            Image that generates the spikes -
        d : array, float32, shape (n_l, n_pix)
            Dictionary used to infer latent factors
        s_gen_name : str
            Name of image (eg. label)
        dt : float
            Timestep for Simulation
        n_t : int
            Number of timesteps of Simulation
        l_n : int
            Linear dimension of neuron array
        a : float
            Image pixel spacing / receptive field spacing
        LAMBDA: float
            Strength of sparse prior
        save_mode : bool
            True if you want to save the data
        n_itr : int
            Number of iterations to break the EM into
        s_gen_name : str
            Name for the generating image
        motion_gen_mode : str
            Method to generate motion. Either Diffusion or Experiment
        dc_gen : float
            Diffusion constant for generating motion
        motion_prior : str
            Prior to use for infering motion
        dc_infer : float
            Diffusion constant for inference
        output_dir : str
            Files saved to 'output/output_dir' If none, uses a time string

        Checks for consistency L_I ** 2 = N_pix

        Note that changing certain parameters without reinitializing the class
        may have unexpected effects (because the changes won't necessarily
        propagate to subclasses.
        """

        self.data = {}

        if output_dir is None:
            output_dir = time_string()
        self.output_dir = os.path.join('output/', output_dir)

        self.save_mode = save_mode  # If true, save results after EM completes
        print 'The save mode is ' + str(save_mode)
        self.d = d.astype(theano.config.floatX)  # Dictionary
        self.d_scl = np.sqrt((d ** 2).sum(1).mean())
        self.n_l, self.n_pix = d.shape
        # n_l - number of latent factors
        # n_pix - number of pixels in the image

        self.s_gen = s_gen.astype(theano.config.floatX)
        self.s_gen_name = s_gen_name
        self.l_i = s_gen.shape[0]  # Linear dimension of the image

        if not self.l_i ** 2 == self.n_pix:
            raise ValueError('Mismatch between dictionary and image size')

        # Simulation Parameters
        self.dt = dt  # Simulation timestep
        self.dc_gen = dc_gen  # Diffusion Constant for Generating motion
        self.dc_infer = dc_infer  # Diffusion Constant for Infering motion
        print 'The diffusion constant is ' + str(self.dc_gen)
        self.l0 = 10.
        self.l1 = 100.

        # Problem Dimensions
        self.n_t = n_t  # Number of time steps
        self.l_n = l_n  # Linear dimension of neuron receptive field grid

        # Image Prior Parameters
        self.gamma = 100.  # Pixel out of bounds cost parameter
        self.LAMBDA = LAMBDA  # the sparse prior is delta (S-DA) + LAMBDA * |A|

        # EM Parameters
        # M - Parameters (ADADELTA)
        # self.rho = 0.4
        # self.eps = 0.001

        # M - parameters (FISTA)
        self.fista_c = 0.5  # Constant to multiply fista L
        self.n_g_itr = 5
        self.n_itr = n_itr

        # E Parameters (Particle Filter)
        self.n_p = 15  # Number of particles for the EM

        # Initialize pixel and LGN positions
        self.a = a  # pixel spacing
        # Position of pixels
        self.XS = np.arange(- self.l_i / 2, self.l_i / 2)
        self.YS = np.arange(- self.l_i / 2, self.l_i / 2)
        self.XS = self.XS.astype('float32')
        self.YS = self.YS.astype('float32')
        self.XS *= self.a
        self.YS *= self.a

        # Position of LGN receptive fields
        self.init_rf_centers()

        # Variances of Gaussians for each pixel
        self.Var = 0.25 * np.ones((self.l_i,)).astype('float32')

        # Gain factor (to be set later)
        self.G = 1.

        # X-Position of retina (batches, timesteps), batches used for inference
        # only
        self.XR = np.zeros((1, self.n_t)).astype('float32')
        # Y-Position of retina
        self.YR = np.zeros((1, self.n_t)).astype('float32')

        # Spikes (1 or 0)
        self.R = np.zeros((self.n_n, self.n_t)).astype('float32')

        # Sparse Coefficients
        self.A = np.zeros((self.n_l,)).astype('float32')

        # Shapes of other variables used elsewhere

        # Pixel values for generating image (same shape as estimated image)
        #        self.s_gen = np.zeros((self.l_i, self.l_i)).astype('float32')
        # Assumes that the first dimension is 'Y'
        #    and the second dimension is 'X'

        # Weighting for batches and time points
        # self.Wbt = np.ones((self.n_b, self.n_t)).astype('float32')

        # Dictionary going from latent factors to image
        #       self.d = np.zeros((self.n_l, self.n_pix)).astype('float32')

        self.tc = TheanoBackend(self.XS, self.YS,
                                self.XE, self.YE,
                                self.IE, self.Var,
                                self.A, self.d)
        self.set_gain_factor()

        if motion_gen_mode == 'Diffusion':
            self.pg = DiffusionPathGenerator(
                self.n_t, self.l_i, self.dc_gen, self.dt)
        elif motion_gen_mode == 'Experiment':
            self.pg = ExperimentalPathGenerator(
                self.n_t, 'data/resampled_paths.mat', self.dt)
        else:
            raise ValueError('motion_gen_mode must' +
                             'be Diffusion of Experiment')

        self.motion_prior = motion_prior
        self.init_particle_filter()

        if (self.save_mode):
            self.init_output_dir()

    def gen_data(self):
        """
        Generates a path and spikes
        Builds a dictionary saving these data
        """
        # Generate Path
        path = self.pg.gen_path()
        self.XR[0, :] = path[0]
        self.YR[0, :] = path[1]

        self.calculate_inner_products()

        self.gen_spikes()
        self.pf.Y = self.R.transpose()  # Update reference to spikes for PF
        # TODO: EWW

        if self.save_mode:
            self.build_param_and_data_dict()

    def init_rf_centers(self):
        """
        Initialize the centers of the receptive fields of the neurons
        """
        self.n_n = 2 * self.l_n ** 2
        self.XE, self.YE = np.meshgrid(
            np.arange(- self.l_n / 2, self.l_n / 2),
            np.arange(- self.l_n / 2, self.l_n / 2)
            )

        self.XE = self.XE.ravel().astype('float32')
        self.YE = self.YE.ravel().astype('float32')
#        self.XE *= self.a
#        self.YE *= self.a

        def double_array(m):
            """
            m - 1d array to be doubled
            :rtype : array that is two copies of m concatenated
            """
            l = m.shape[0]
            res = np.zeros((2 * l,))
            res[0:l] = m
            res[l: 2 * l] = m
            return res

        self.XE = double_array(self.XE)
        self.YE = double_array(self.YE)

        # Identity of LGN cells (ON = 0, OFF = 1)
        self.IE = np.zeros((self.n_n,)).astype('float32')
        self.IE[0: self.n_n / 2] = 1

    def set_gain_factor(self):
        """
        Sets the gain factor so that an image with pixels of intensity 1
            results in spikes at the maximum firing rate
        """
        self.G = 1.
        Ips, FP = self.tc.RFS(
            np.ones_like(self.s_gen), self.XR, self.YR,
            self.l0, self.l1, self.dt, self.G)
        self.G = (1. / Ips.max()).astype('float32')

    def init_particle_filter(self):
        """
        Initializes the particle filter class
        Requires spikes to already be generated
        """
        # Define necessary components for the particle filter
        if self.motion_prior == 'PositionDiffusion':
            # Diffusion
            D_H = 2  # Dimension of hidden state (i.e. x,y = 2 dims)
            sdev = np.sqrt(self.dc_infer * self.dt / 2) * np.ones((D_H,))
            ipd = pf.GaussIPD(D_H, self.n_n, sdev * 0.001)
            tpd = pf.GaussTPD(D_H, self.n_n, sdev)
            ip = pf.GaussIP(D_H, sdev * 0.001)
            tp = pf.GaussTP(D_H, sdev)
            lp = PoissonLP(self.n_n, D_H, self.l0, self.l1,
                           self.dt, self.G, self.tc.spike_energy)

        elif self.motion_prior == 'VelocityDiffusion':
            # FIXME: save these params
            D_H = 4   # Hidden state dim, x,y,vx,vy
            v0 = 30.  # Initial Estimate for velocity
            dcv = 6.  # Velocity Diffusion Constant
            st = np.sqrt(dcv * self.dt)

            eps = 0.00001  # Small number since cannot have exact zero
            sigma0 = np.array([eps, eps, v0, v0])  # Initial sigmas
            sigma_t = np.array([eps, eps, st, st])  # Transition sigmas

            # Transition matrix
            A = np.array([[1, 0, self.dt, 0],
                          [0, 1, 0, self.dt],
                          [0, 0, 1, 0],
                          [0, 0, 0, 1]])

            ipd = pf.GaussIPD(D_H, self.n_n, sigma0)
            tpd = pf.GaussTPD(D_H, self.n_n, sigma_t, A=A)
            ip = pf.GaussIP(D_H, sigma0)
            tp = pf.GaussTP(D_H, sigma_t, A=A)
            lp = PoissonLP(self.n_n, D_H, self.l0, self.l1,
                           self.dt, self.G, self.spike_energy)

        else:
            raise ValueError(
                'Unrecognized Motion Prior ' + str(self.motion_prior))

        self.pf = pf.ParticleFilter(ipd, tpd, ip, tp, lp,
                                    self.R.transpose(), self.n_p)

    def gen_spikes(self):
        """
        Generate LGN responses given the path and the image
        """
        self.R[:, :] = self.tc.spikes(
            self.s_gen, self.XR, self.YR, self.l0,
            self.l1, self.dt, self.G)[0]
        print 'Mean firing rate ' + str(self.R.mean() / self.dt)

    def true_costs(self):
        """
        Prints out the negative log-likelihood of the observed spikes given the
            image and path that generated them
        Note that the energy is normalized by the number of timesteps
        """
        print 'Pre-EM testing'
        args = (self.XR, self.YR, self.R, self.Wbt,
                self.l0, self.l1, self.dt, self.G,
                self.gamma, self.LAMBDA)

        out = self.tc.costs(*args)

        E, E_bound, E_R, E_rec, E_sp = out
        print ('Costs of underlying data ' +
               str((E / self.n_t, E_rec / self.n_t)))

    def reset(self):
        """
        Resets the class between EM runs
        """
        self.pf.reset()
        # self.c.reset()
        self.data = {}
        self.tc.reset_image_estimate()
        self.tc.reset_m_aux()

    def run_E(self, t):
        """
        Runs the the particle filter until it has run a total of t time steps
        t - number of timesteps
        The result is saved in self.pf.XS,WS,means
        """
        if t > self.n_t:
            raise IndexError('Maximum simulated timesteps exceeded in E step')
        if self.pf.t >= t:
            raise IndexError(
                'Particle filter already run past given time point')
        while self.pf.t < t:
            self.pf.advance()
        self.pf.calculate_means_sdevs()

        print 'Path SNR ' + str(SNR(self.XR[0][0:t], self.pf.means[0:t, 0]))

    def run_M(self, t, N_g_itr=5):
        """
        Runs the maximization step for the first t time steps
        resets the values of auxillary gradient descent variables at start
        t - number of time steps
        result is saved in t_A.get_value()
        """
        self.tc.reset_m_aux()
        print ('Total Energy / t | Bound. Energy / t ' +
               '| Spike Energy / t | + Sparse E / t |  SNR')
        fista_l = self.tc.calculate_L(
            t, self.n_n, self.l0, self.l1, self.dt, self.d_scl, self.fista_c)

        for v in range(N_g_itr):
            args = (self.pf.XS[0:t, :, 0].transpose(),
                    self.pf.XS[0:t, :, 1].transpose(),
                    self.R[:, 0:t], self.pf.WS[0:t].transpose(),
                    self.l0, self.l1, self.dt,
                    self.G, self.gamma, self.LAMBDA,
                    # self.rho, self.eps)
                    fista_l)

            out = self.tc.img_grad(*args)
            self.img_SNR = SNR(self.s_gen, self.tc.image_est())

            print self.print_costs(out, t) + str(self.img_SNR)

    @staticmethod
    def print_costs(out, t):
        """
        Prints costs given output of img_grad
        Cost divided by the number of timesteps
        out - tuple containing the differet costs
        t - number to timesteps
        """
        strg = ''
        for item in out:
            strg += str(item / t) + ' '
        return strg

    def run_EM(self):
        """
        Runs full expectation maximization algorithm
        self.N_itr - number of iterations of EM
        self.N_g_itr - number of gradient steps in M step
        Saves summary of run info in self.data
        Note running twice will overwrite this run info
        """
        self.tc.reset_image_estimate()

        EM_data = {}

        print '\n' + 'Running full EM'

        for u in range(self.n_itr):
            t = self.n_t * (u + 1) / self.n_itr
            print ('\n' + 'Iteration number ' + str(u) +
                   ' Running up time = ' + str(t))

            # Run E step
            self.run_E(t)

            # Run M step
            if u <= 2:
                c = 4
            else:
                c = 1
            self.run_M(t, N_g_itr=self.n_g_itr * c)

            iteration_data = {'time_steps': t, 'path_means': self.pf.means,
                              'path_sdevs': self.pf.sdevs,
                              'image_est': self.tc.image_est(),
                              'coeff_est': self.tc.get_A()}

            EM_data[u] = iteration_data

        if self.save_mode:
            self.data['EM_data'] = EM_data

    def init_output_dir(self):
        """
        Create an output directory
        """
        if not os.path.exists('output'):
            os.mkdir('output')

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

    def build_param_and_data_dict(self):
        """
        Creates a dictionary, self.data, that has all of the parameters
            of the model
        """
        # Note it is important to create a new dictionary here so that
        # we reset the data dict after generating new data
        self.data = {'DT': self.dt,
                     'DC_gen': self.dc_gen,
                     'DC_infer': self.dc_infer,
                     'L0': self.l0,
                     'L1': self.l1,
                     'GAMMA': self.gamma,
                     'LAMBDA': self.LAMBDA,
                     'D': self.d,
                     'N_L': self.n_l,
                     'a': self.a,
                     'N_T': self.n_t,
                     'L_I': self.l_i,
                     'L_N': self.l_n,
                     # 'Rho': self.rho,
                     # 'Eps': self.eps,
                     'N_g_itr': self.n_g_itr,
                     'N_itr': self.n_itr,
                     'N_P': self.n_p,
                     'XS': self.XS, 'YS': self.YS,
                     'XE': self.XE, 'YE': self.YE,
                     'Var': self.Var,
                     'G': self.G,
                     'XR': self.XR, 'YR': self.YR,
                     'IE': self.IE,
                     'actual_motion_mode': self.pg.mode(),
                     'S_gen': self.s_gen, 'S_gen_name': self.s_gen_name,
                     'R': self.R,
                     'Ips': self.Ips, 'FP': self.FP,
                     'motion_prior': self.motion_prior}

    def save(self):
        """
        Saves information relevant to the EM run
        data.pkl - saves dictionary with all data relevant to EM run
        (Only includes dict for EM data if that was run)
        Returns the filename
        """
        if not self.save_mode:
            raise RuntimeError('Need to enable save mode to save')

        fn = os.path.join(self.output_dir,
                          'data_' + time_string() + '.pkl')
        pkl.dump(self.data, open(fn, 'wb'))
        return fn

    def calculate_inner_products(self):
        """
        Calculates the inner products used
        """
        self.Ips, self.FP = self.tc.RFS(
            self.s_gen, self.XR, self.YR, self.l0, self.l1,
            self.dt, self.G)


# def true_path_infer_image_costs(self, N_g_itr = 10):
#        """
#       Infers the image given the true path
#        Prints the costs associated with this step
#        """
#        self.reset_img_gpu()
#        print 'Original Path, infer image'
#        t = self.N_T
#        self.run_M(t)


#    def true_image_infer_path_costs(self):
#        print 'Original image, Infer Path'
#        print 'Path SNR'
#        self.t_S.set_value(self.S)
#        for _ in range(4):
#            self.run_E(self.N_T)
#
#        if self.debug:
#            self.pf.plot(self.XR[0], self.YR[0], self.DT)


