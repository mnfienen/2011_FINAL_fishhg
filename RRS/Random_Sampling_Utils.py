#! /usr/local/bin/python
import random


def sample_without_replacement(n,r):
    '''
    http://code.activestate.com/recipes/272884-random-samples-without-replacement/
    
    to call this, must use LIST to get the results --> for example:
    a = list(sampling_without_replacement(150,3))
    This will result in a list of 3 integers between 0 and 150
    '''
    # select r randomly chosen, sorted integers from [0,n]
    _rand = random.random # aliasing speed hack
    popn = n
    for samp in xrange(r,0,-1):
        cumprob = 1.0
        x = _rand()
        while x < cumprob:
            cumprob -= cumprob * samp/popn
            popn -= 1
        yield n-popn-1
        