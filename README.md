# Leakage Models

## Introduction

These scripts prove experimental validation of the bounds provided in Proposition 3 of [this paper](https://doi.org/10.1007/978-3-030-26948-7_24). The bounds are proven; the scripts serve to show that the proofs are indeed correct and may be useful in other contexts as well (e.g. they can be adapted to other noise distributions than the one we considered).

## How to use

This script can only be used on [Sagemath](https://www.sagemath.org/). Here is a small example for the AES field (n = 8) and a Gaussian noise of standard deviation 50.

```python
sage: attach("leakage.sage")
sage: compute_noise(8, sigma=50)

Parameters
==========
Size of the field   : N     =  256
Max. Hamming weight : n     =  8
Standard deviation  : sigma =  50
Tailcut rate        : tau   =  19

Warning: not in the asymptotic regime for RE
(but asymptotic regime for ARE, SD and EN is OK).
Try setting (sigma > 5 * n * tau) to be in that asymptotic regime too.

RE(X|Y)                        = 2.97075320220171
RE(X|Y) * sigma / (tau * n)    = 0.977221448092669
1 / 2                          = 0.500000000000000

ARE(X|Y)                       = 0.066469770504
ARE(X|Y) * sigma / n           = 0.41543606565
1 / sqrt(2 * pi)               = 0.398942280401

SD(X|Y)                        = 0.00877843993091969
SD(X|Y) * sigma / sqrt(n)      = 0.155182360084802
1 / (2 * pi)                   = 0.159154943092

EN(X|Y)                        = 0.00141124735097700
EN(X|Y) * sigma * sqrt(N / n)  = 0.399161028722954
1 / sqrt(2 * pi)               = 0.398942280401

```

We can see that formulae for ARE(X|Y), SD(X|Y) and EN(X|Y) are close to the theoretical estimates. The formula for RE(X|Y) has a slower convergence rate than the other ones, and the code raises a warning if parameters are too low for the asymptotic regime to kick in.