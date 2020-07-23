
import sympy as sp


class opt_tech():
    """
    This is a class for the two main classes of optimization techniques for solving 
    1-D unconstraind maximization problems numerically

    ...

    Attributes
    ----------
    fx : sympy class obj
        the objective function you want to find its optimum point
    itr_no : int
        the number of iterations you will go through (default 1)
    prec : int
        the number of decimal places after the decimal point (default 3)
    """

    def __init__(self, fx, itr_no=1, prec=3):
        """
        The constructor for opt_tech class. 

        ...

        Parameters
        ----------
        fx : sympy class obj
            the objective function you want to find its optimum point
        itr_no : int
            the number of iterations you will go through (default 1)
        prec : int
            the number of decimal places after the decimal point (default 3)
        """

        self.fx = fx
        self.itr_no = itr_no
        self.prec = prec
        self.x = sp.Symbol("x")

    def newtons_method(self, xo=0):
        """
        A direct root finding Optimization Technique for 1-D constrained problem.

        ...

        Parameters
        ----------
        xo : int
            the initial guess needed to compute the optimum point (default 0)
        """

        print()
        print()
        print("Newton's Method :")
        print("_________________")

        # to get the 1st derivative for the objective function
        df = sp.diff(self.fx)
        # to get the 2nd derivative for the objective function
        ddf = sp.diff(df)
        # set the initial value to the initial guess
        xi_old = xo

        # find the optimal value for the given number of iterations
        for i in range(self.itr_no):
            # calculate the imporved approximation of xi_old
            xi_new = xi_old - df.subs(self.x, xi_old) / \
                ddf.subs(self.x, xi_old)
            # calculate the relative error = |xi - xi+1| / |xi+1|
            relative_error = abs(xi_new-xi_old)/abs(xi_new)
            print("x{} = {:.{prec}f}, relative error = {:.{prec}f}".format(
                i+1, xi_new, relative_error, prec=self.prec))
            # update the old approximation to be equal the new one
            xi_old = xi_new
        print()

        # calculate the 1st derivative value with our final approximation
        f_dash = df.subs(self.x, xi_old)
        # find if our final approximation is close to zero (accepted) or not
        note = "(acceptably small)" if abs(f_dash) <= 0.01 else ""
        if note:
            print("f'({:.{prec}f}) = {:.{prec}f}  {}".format(
                xi_old, f_dash, note, prec=self.prec))
            print()
            print("fmax = f({:.{prec}f}) = {:.{prec}f}".format(
                xi_old, self.fx.subs(self.x, xi_old), prec=self.prec))
        else:
            print("{:.{prec}f} is rejected".format(xi_old, prec=self.prec))
        print("_________________")

    def golden_section(self, xl, xu):
        """
        An Elimination Optimization Technique for 1-D constrained problem.

        ...

        Parameters
        ----------
        xl : int
            the lower bound for the unimodal objective function interval
        xu : int
            the upper bound for the unimodal objective function interval
        """

        print()
        print()
        print("Golden Section Method :")
        print("_______________________")
        print()
        # creating the table format for printing
        table_format = ["i", "xl", "x2", "x1", "xu",
                        "fx2", "fx1", "xu-xl", "errBound"]
        print("|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|{:>{tab}}|".format(
            *table_format, tab=2*self.prec))

        # we can calculate all the lengths directly using the formula
        # xu(i)-xl(i) = ( xu(0)-xl(0) ) * .618^(i-1)
        length = [round((xu-xl)*.618**i, self.prec)
                  for i in range(self.itr_no)]

        # we can calculate all the lengths directly using the formula
        # xu(i)-xl(i) = ( xu(0)-xl(0) ) * .618^(i+1)
        error_bound = [round((xu-xl)*.618**(i+2), self.prec)
                       for i in range(self.itr_no)]

        # find the optimal value for the given number of iterations
        for i in range(self.itr_no):
            # calculate our comparison points and their values
            x2 = xu - .618*(xu-xl)
            x1 = xl + .618*(xu-xl)

            fx2 = self.fx.subs(self.x, x2)
            fx1 = self.fx.subs(self.x, x1)

            # creating the data format for printing
            data = [i+1, xl, x2, x1, xu, fx2, fx1, length[i], error_bound[i]]
            print("|{:>{tab}}|{:>{tab}.{prec}f}|{:>{tab}.{prec}f}|{:>{tab}.{prec}f}|{:>{tab}.{prec}f}|{:>{tab}.{prec}f}|{:>{tab}.{prec}f}|{:>{tab}.{prec}f}|{:>{tab}.{prec}f}|".format(
                *data, tab=2*self.prec, prec=self.prec))

            # find which interval to eliminate
            # if f(x1) < f(x2), then [xl, x2] doesn't contain the maximum point
            if fx1 > fx2:
                xl, x2 = x2, x1
            # if f(x2) > f(x1), then [x1, xu] doesn't contain the maximum point
            # in case, f(x1) == f(x2), then both [xl, x2] and [x1, xu] doesn't
            # contain the maximum point so we can safely do either options
            # without losing the maximum point
            else:
                xu, x1 = x1, x2
        print("_________________")


# the objective function
fx = 2*sp.sin(sp.Symbol("x")) - .1*(sp.Symbol("x")**2)
# the precision value (number of decimal places)
prec = 3
# the number of iterations
itr_no = 3

# initializing our optimization techniques class
opt_tech = opt_tech(fx, itr_no, prec)

# our initial guess
xo = 2.5

# Calling our Newton's Method
opt_tech.newtons_method(xo)

# the lower bound
xl = 0
# the upper bound
xu = 4
# the precision value (number of decimal places)
opt_tech.prec = 4
# the number of iterations
opt_tech.itr_no = 8

# Calling our Golden Section Method
opt_tech.golden_section(xl, xu)
