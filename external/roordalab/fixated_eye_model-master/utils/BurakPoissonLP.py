import numpy as np
import particle_filter_new as PF

# For the particle filter module, this class mediates the emission
# probabilities


class PoissonLP(PF.LikelihoodPotential):

    """
    Poisson likelihood potential for use in the burak EM
    """

    def __init__(self, D_O, D_H, L0, L1, DT, G, spike_energy):
        """
        D_O - dimension of output (number of neurons)
        D_H - dimesion of hidden state
        L0 - lower firing rate
        L1 - higher firing rate
        DT - time step
        G - gain factor
        spike_energy - function handle to spike energy function
                         see prob for desired arguments for spike_energy
        """
        if not (D_H == 2 or D_H == 4):
            raise ValueError('D_H should be 2 or 4')
        PF.LikelihoodPotential.__init__(self, D_H, D_O)
        self.L0 = L0
        self.L1 = L1
        self.DT = DT
        self.G = G
        self.spike_energy = spike_energy

    def prob(self, Yc, Xc):
        """
        Gives probability p(R_t|X_t, S),
        Yc - Observed spiking pattern - R_t in form (D_O,)
        Xc - Current hidden state in form (N_S, D_H)
        Assumes x position, y position are indices 0, 1 along the
            D_H dimension
        """
        self.N_S, _ = Xc.shape
        # Note to pass to theano function need:
        #  XR -> (N_batches, N_timesteps)
        #  R -> (N_b, N_pix, N_t)
        # here, we have one timestep, and particles correspond to batches,
        #  and given how R is, we need to broadcast the spikes over batches

        _XR = np.zeros((self.N_S, 1)).astype('float32')
        _YR = np.zeros((self.N_S, 1)).astype('float32')
        _XR[:, 0] = Xc[:, 0]
        _YR[:, 0] = Xc[:, 1]

        _R = np.zeros((self.D_O, 1)).astype('float32')
        _R[:, 0] = Yc
        # to pass to spike_prob function, N_P batches, 1 time step
        Es = - \
            self.spike_energy(_XR, _YR, _R, self.L0, self.L1, self.DT, self.G)
        Es = Es - Es.mean()
        return np.exp(Es)
