def compute_noise(n, sigma, tau=19):
    """
    This function takes as input:

    - a field bitsize n
    - a standard deviation sigma
    - a tailcut rate tau (by default set to 19)

    And computes the {RE,ARE,SD,EN}-noisiness of the random variables X and Y:

    - X is uniformly distributed in [0, ..., 2 ** n - 1]
    (hence its Hamming weight HW(X) is in [0, ..., n])
    - Y = f(X) is equal to HW(X) + a Gaussian of standard deviation sigma

    A warning is raised if parameters do not fall into the asymptotic regime.
    """

    # p_x is the probability that X = x, for x in [0, ..., 2 ** n - 1]
    p_x = (2. ** (- n))
    # p_hw(k) is the probability that HW(X) = k, for k in [0, ..., n]
    p_hw(k) = binomial(n, k) * (2. ** (- n))
    # p_y_given_k(y, k) is the probability (density) of Y at the point y,
    # given that HW(X) = k
    p_y_given_k(y, k) = exp(- (k - y) * (k - y) / (2. * sigma * sigma)) / (float(sigma * sqrt(2. * pi)))
    # p_y(y) is the probability (density) of Y at the point y
    p_y(y) = sum(p_hw(k) * p_y_given_k(y, k) for k in range(n + 1))

    # Print the parameters
    print ""
    print "Parameters"
    print "=========="
    print "Size of the field   : N     = ", 2 ** n
    print "Max. Hamming weight : n     = ", n
    print "Standard deviation  : sigma = ", sigma
    print "Tailcut rate        : tau   = ", tau
    print ""

    # Check if we are in the asymptotic regime
    if (sigma < 5 * n):
        print "Warning: not in the asymptotic regime for RE, ARE, SD and EN."
        print "Try setting (sigma > 5 * n)"
        print ""
    elif (sigma < 5 * n * tau):
        print "Warning: not in the asymptotic regime for RE"
        print "(but asymptotic regime for ARE, SD and EN is OK)."
        print "Try setting (sigma > 5 * n * tau) to be in that asymptotic regime too."
        print ""

    # Compute the noise with the relative error RE
    y0 = float(n + tau * sigma)
    k0 = n
    RE = abs(p_y_given_k(y0, k0) / p_y(y0) - 1)
    print "RE(X|Y)                        =", RE
    # In the asymptotic regime, this should be equivalent to a constant
    print "RE(X|Y) * sigma / (tau * n)    =", RE * sigma / (tau * n)
    print "1 / 2                          =", 1. / 2
    print ""

    # Compute the noise with the average relative error ARE
    k0 = 0
    k1 = n
    p_pmi_left(y) = abs(p_y_given_k(y, k0) - p_y(y))
    p_pmi_right(y) = abs(p_y_given_k(y, k1) - p_y(y))
    p_pmi(y) = max_symbolic(abs(p_y_given_k(y, k0) - p_y(y)), abs(p_y_given_k(y, k1) - p_y(y)))
    ARE = numerical_integral(p_pmi(y), - tau * sigma, n + tau * sigma, params=[y])[0]
    print "ARE(X|Y)                       =", ARE
    # In the asymptotic regime, this should be equivalent to a constant
    print "ARE(X|Y) * sigma / n           =", ARE * sigma / n
    print "1 / sqrt(2 * pi)               =", float(1 / sqrt(2 * pi))
    print ""

    # Compute the noise with the statistical distance SD
    SD_integrand = p_x * sum(binomial(n, k) * abs(p_y(y) - p_y_given_k(y, k)) for k in range(n + 1))
    SD = 0.5 * numerical_integral(SD_integrand, - tau * sigma, n + tau * sigma, params=[y])[0]
    print "SD(X|Y)                        =", SD
    # In the asymptotic regime, this should be equivalent to a constant
    print "SD(X|Y) * sigma / sqrt(n)      =", SD * sigma / sqrt(1. * n)
    print "1 / (2 * pi)                   =", float(1 / (2 * pi))
    print ""

    # Compute the noise with the Euclidean norm EN
    EN_integrand = sqrt(sum(binomial(n, k) * ((p_y(y) - p_y_given_k(y, k)) ** 2) for k in range(n + 1)))
    EN = p_x * numerical_integral(EN_integrand, -tau * sigma, n + tau * sigma, params=[y])[0]
    print "EN(X|Y)                        =", EN
    # In the asymptotic regime, this should be equivalent to a constant
    print "EN(X|Y) * sigma * sqrt(N / n)  =", EN * sigma * sqrt(2 ** n) / sqrt(1. * n)
    print "1 / sqrt(2 * pi)               =", float(1 / sqrt(2 * pi))
    print ""
