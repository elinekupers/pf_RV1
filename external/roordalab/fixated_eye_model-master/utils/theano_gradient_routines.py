# Helper functions for gradient descent in theano
# Implements:   -- SGD with Momentum
#               -- ADADelta

# Functions take in an energy function and a parameter
# Returns auxillary variables for descent routines
#   to minimize wrt that parameter in an updates list
#   that can be passed to a theano function


import numpy as np
import matplotlib.pyplot as plt
import theano
import theano.tensor as T
from collections import OrderedDict


def momentum(t_E, t_X, *t_params):
    """
    t_X - parameter to do gradient descent on (shared variable)
    t_E - energy to minimize (function of t_X)
    t_params - tuple of parameters (theano scalars)
                (t_Eta, t_M)
    Return
    updates - OrderedDict of updates for theano function
        X_mom - shared variable accumulating linear combination of
            past gradients
        X - appropriate update X -> X + dX
    """
    t_Eta = t_params[0]
    t_M = t_params[1]
    p_name = t_X.name
    e_name = t_E.name
    g_name = 't_' + e_name + '_' + p_name

    t_g_E_X = T.grad(t_E, t_X)
    t_g_E_X.name = g_name
    shape = t_X.get_value().shape

    t_X_mom = theano.shared(np.zeros(shape).astype('float32'),
                            p_name + '_mom')

    updates = OrderedDict()
    t_X_mom_temp = t_Eta * t_g_E_X + t_M * t_X_mom
    updates[t_X_mom] = t_X_mom_temp
    updates[t_X] = t_X - t_X_mom_temp
    return updates


def ada_delta(t_E, t_X, *t_params):
    """
    t_E - some energy function to be minimized
    t_X - some parameter of the energy function, as a shared variable
    t_params - tuple of theano scalars for adadelta parameters
            (t_Rho, t_Eps)
            t_Rho - decay rate
            t_Eps - regularization for square roots
    Returns:
    updates - OrderedDict to be passed to a theano function
        X_Eg2 - exponential moving average of gradients (variable by variable)
        X_Edx2 - EMA of gradient steps (separated by variable)
        update rule for t_X


    """
    t_Rho = t_params[0]
    t_Eps = t_params[1]

    p_name = t_X.name  # Parameter name
    e_name = t_E.name
    g_name = 't_' + e_name + '_' + p_name
    # Calculate the Gradient
    t_g_X = T.grad(t_E, t_X)
    t_g_X.name = g_name

    # Create updates dictionary (theano wants an ordered dictionary)
    updates = OrderedDict()
    shape = t_X.get_value().shape
    # Create shared variables for ADADELTA
    t_Eg2 = theano.shared(np.zeros(shape).astype('float32'),
                          p_name + '_Eg2')
    t_Edx2 = theano.shared(np.zeros(shape).astype('float32'),
                           p_name + '_Ed' + p_name + '2')

    # Start calculating updates, note all updates to shared variables
    #    happen at once, so we need to do it this way:

    t_Eg2_temp = t_Rho * t_Eg2 + (1. - t_Rho) * (t_g_X ** 2)
    updates[t_Eg2] = t_Eg2_temp

    t_dX = - T.sqrt(t_Edx2 + t_Eps) / T.sqrt(t_Eg2_temp + t_Eps) * t_g_X

    updates[t_Edx2] = t_Rho * t_Edx2 + (1. - t_Rho) * (t_dX ** 2)

    updates[t_X] = t_X + t_dX

    return updates


def main():
    """
    Gives sample code for using these routines
    """

    # Define a minimization problem
    N_D = 2
    X0 = np.ones(N_D)  # Initial value for parameters
    A = np.diagflat([1, 5])

    # Create energy function with theano variables
    t_X = theano.shared(X0, 'X')  # Parameters
    t_A = T.matrix('A')  # Fixed Data
    t_E = T.sum(t_X.dimshuffle(0, 'x') * t_A * t_X.dimshuffle('x', 0))
    t_E.name = 'E'

    # ADADelta parameters
    Rho = 0.95
    Eps = 0.001
    params_ada = (Rho, Eps)

    # ADADelta parameters as theano variables
    t_Rho = T.scalar('Rho')  # learning rate
    t_Eps = T.scalar('Eps')
    t_params_ada = (t_Rho, t_Eps)

    t_updates_ada = ada_delta(t_E, t_X, *t_params_ada)

    ada_delta_minimize = theano.function(inputs=[t_A,
                                                 t_Rho, t_Eps],
                                         outputs=[t_E] + t_updates_ada.keys(),
                                         updates=t_updates_ada)

    # Momentum parameters
    Eta = 0.01
    M = 0.5
    params_mom = (Eta, M)

    # Momentum parameters as theano variables
    t_Eta = T.scalar('eta')  # learning rate
    t_M = T.scalar('m')  # momentum term
    t_params_mom = (t_Eta, t_M)

    t_updates_mom = momentum(t_E, t_X, *t_params_mom)

    momentum_minimize = theano.function(inputs=[t_A,
                                                t_Eta, t_M],
                                        outputs=[t_E] + t_updates_mom.keys(),
                                        updates=t_updates_mom)

    # ADADelta
    t_X.set_value(X0)
    N_itr = 100
    Xs = np.zeros((N_itr, N_D))  # Store values of parameters
    for i in range(N_itr):
        E, Eg2, Edx2, X = ada_delta_minimize(A, *params_ada)
        Xs[i] = X

    plt.scatter(Xs[:, 0], Xs[:, 1])
    plt.title('ada_delta with ' + str(N_itr) + ' steps')
    plt.plot([0.], [0.], 'ro')
    plt.show()

    # Momentum
    t_X.set_value(X0)
    N_itr = 100
    Xs = np.zeros((N_itr, N_D))  # Store values of parameters
    for i in range(N_itr):
        E, X_mom, X = momentum_minimize(A, *params_mom)
        Xs[i] = X

    plt.scatter(Xs[:, 0], Xs[:, 1])
    plt.title('momentum with ' + str(N_itr) + ' steps')
    plt.plot([0.], [0.], 'ro')
    plt.show()

    return 0

if __name__ == "__main__":
    main()
