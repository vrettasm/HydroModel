# Imports:
from numpy import log, array, asarray

# This module includes:
__all__ = ['safelog']


# SAFELOG:
def safelog(x=None):
    ''' Safe logarithm
    
    Description:
        This (helper) function prevents the computation of very small or very  large
        values of logarithms that would lead to -/+ inf, by setting predefined LOWER
        and UPPER bounds. The bounds are set as follows:
        
            - LOWER = 1.0E-300
            - UPPER = 1.0E+300
            
        It is assumed that the input values lie within this range.
    
    Example:
        >> numpy.log(1.0E-350)
        >> -inf
        >>
        >> safelog(1.0E-350)
        >> -690.77552789821368
    
    Input:
        x : input array (N x M).
        
    Output:
        x : the log(x) after the values of x have been filtered (N x M).
    
    Copyright (c) Michail D. Vrettas, PhD - March 2015.
    
    Last Update: March 2015.
    
    References:
        (N/A)
    '''

    # Prevent empty input.
    if (x==None):
        print(" [safelog] in debug mode ... exiting:", end="")
        return None

    # Define LOWER and UPPER bounds.
    _LWR_Bound = 1.0E-300
    _UPR_Bound = 1.0E+300

    # Make sure input is an array.
    x = asarray(x)

    # Overflow Error can happen if the input is "int64"
    try:
        # Check for scalar.
        if (x.ndim == 0):
            if (x < _LWR_Bound):
                x = _LWR_Bound
            elif (x > _UPR_Bound):
                x = _UPR_Bound
                # _end_if
        else:
            # Check Lower/Upper bounds
            x[x < _LWR_Bound] = _LWR_Bound
            x[x > _UPR_Bound] = _UPR_Bound
    except OverflowError:
        return safelog(array(x, dtype='float'))

    # Return the log of the filtered input.
    return log(x)
    # == END OF FILE ==
