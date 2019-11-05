import numpy as np
from scipy.io import loadmat
import abc


class Center:

    def __init__(self, Lx, DC, DT):
        """
        Class that Implements a diffusing center in a box of size Lx

        Lx = linear dimension of square to diffuse in
        DC = diffusion constant
        DT = timestep size
        x = coordinates of the current location of the random walk,
            initialized as [0, 0]
        Initializes Center Object
        """
        self.Lx = Lx
        self.DC = DC
        self.m0 = np.array([0, 0], dtype='float64')  # current position
        self.DT = DT
        # The diffusion is biased towards the center by taking a product of
        # gaussians
        # A product of gaussians is also a gaussian with mean, sdev given as
        # (mn, sn)
        self.m1 = np.array([0, 0], dtype='float64')  # center of image
        # Standard deviation for diffusion
        self.s0 = np.sqrt(self.DC * self.DT)
        self.s1 = self.Lx / 4  # Standard Deviation for centering gaussian
        self.sn = 1 / np.sqrt(1 / self.s0 ** 2 + 1 / self.s1 ** 2)

    def advance(self):
        """
        Updates location according to a random walk that stays within a box
        """

        self.mn = (
            self.m0 / self.s0 ** 2 + self.m1 / self.s1 ** 2) * self.sn ** 2
        while(True):
            # Note that for 2d diffusion, each component's variance is half the
            #   variance of the overall step length
            temp = self.mn + \
                np.random.normal(size=2, scale=self.sn / np.sqrt(2))
            if (temp[0] > - self.Lx / 2
                    and temp[0] < self.Lx / 2
                    and temp[1] > - self.Lx / 2
                    and temp[1] < self.Lx / 2):
                self.m0 = temp
                break

    def get_center(self):
        return self.m0

    def reset(self):
        self.m0 = np.array([0, 0], dtype='float64')

    def __str__(self):
        return str(self.x)


class PathGenerator():
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, N_T, *args):
        """
        Initialize a 2d Path Generator
        """
        self.N_T = N_T

    @abc.abstractmethod
    def gen_path(self):
        """
        Generates a 2d path with N_T timesteps starting at (0,0)
        """
        return np.zeros((2, self.N_T))

    @abc.abstractmethod
    def mode(self):
        """
        Returns a string describing the path
        """
        return ' '


class DiffusionPathGenerator(PathGenerator):

    def __init__(self, N_T, Lx, DC, DT):
        """
        Creates a path generator that does bounded diffusion
        N_T - number of timesteps to generate path
        Lx - window that the diffusion is restricted to
        DC - diffusion constant
        DT - timestep
        """
        PathGenerator.__init__(self, N_T)
        self.c = Center(Lx, DC, DT)

    def gen_path(self):
        self.c.reset()
        path = np.zeros((2, self.N_T))
        for i in range(self.N_T):
            path[:, i] = self.c.get_center()
            self.c.advance()
        return path

    def mode(self):
        return 'Diffusion'


class ExperimentalPathGenerator(PathGenerator):

    def __init__(self, N_T, filename, DT):
        """
        Creates a path generator that uses real experimental data
        filename - filename for pkl file that contains an array of paths
              data['paths'] = (N_runs, number of timesteps, 2)
              2 - number of dimensions of path eg. x,y
        """
        PathGenerator.__init__(self, N_T)
        self.data = loadmat(filename)
        self.DT = self.data['DT'][0, 0]
        self.paths = self.data['paths']
        self.N_runs, self.N_T_data, _ = self.paths.shape

        if not self.DT == DT:
            raise ValueError('Data timestep doesnt match simulation timestep')
        if self.N_T > self.N_T_data:
            pass  # raise ValueError('Simulation has more timesteps than data')

    def gen_path(self):
        """
        Generate a path from the data
        """
        q = np.random.randint(self.N_runs)
        st = np.random.randint(self.N_T_data - self.N_T)
        pat = self.paths[q, st:(st + self.N_T), :].transpose()
        pat[0] = pat[0] - pat[0][0]
        pat[1] = pat[1] - pat[1][0]
        return pat

    def mode(self):
        return 'Experimental_data'


if __name__ == '__main__':
    fn = 'data/resampled_paths.mat'
    pg = ExperimentalPathGenerator(100, fn, 0.001)
    path = pg.gen_path()
    print path
