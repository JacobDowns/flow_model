import dolfin as df
from ufl import max_value, min_value

def softplus(y1,y2,alpha=1):
    """ 
    The softplus function is a differentiable approximation to the ramp 
    function.  Its derivative is the logistic function. Larger alpha 
    makes a sharper transition.
    """
    return max_value(y1,y2) + (1./alpha)*df.ln(1+df.exp(alpha*(min_value(y1,y2)-max_value(y1,y2))))

def logistic(y1, y0 = df.Constant(1.), k = df.Constant(1.0)):
    """ 
    The logistic function. y0 is the y-value of the sigmoid's midpoint.
    k controls the steepness of the curve. 
    """
    return df.Constant(1.) / (df.Constant(1.) +  df.exp(-k*(y1 - y0)))

def Abs(y):
    """ 
    Absolute value since this inexplicably dissapeared from Fenics.
    """
    return df.sqrt(y**2 + df.Constant(1e-16))
