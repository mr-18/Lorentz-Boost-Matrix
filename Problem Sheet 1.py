'''
Problem sheet 1: creating classes for: a vector of length 4;
                                       a 4 x 4 matrix;
                                       a four-vector
                                       a lorentz boost
'''
class Vector:
    '''
    Class representing a vector of length 4
    '''
    def __init__(self, v0, v1, v2, v3):
        '''
        Initialise the class with its four components
        '''
        self.v = [v0, v1, v2, v3]
    def __str__(self):
        '''
        Return a string to print for this vector
        '''
        return "{}".format(self.v)
    def __repr__(self):
        '''
        Return the representation of this vector
        '''
        return "{}({}, {}, {}, {})".format(self.__class__.__name__,
                                           self.v[0], self.v[1], self.v[2],
                                           self.v[3])
    def __getitem__(self, i):
        '''
        Return the 'i' component of the vector
        '''
        return self.v[i]
    def __setitem__(self, i, s):
        '''
        Set the 'i' component of the vector with the scalar s
        '''
        self.v[i] = s
    def __len__(self):
        '''
        Return the length of the vector
        '''
        return len(self.v)
    def __pos__(self):
        '''
        Return a copy of this vector with the '+' operator applied to each
        element
        '''
        return self.__class__(+self[0], +self[1], +self[2], +self[3])
    def __neg__(self):
        '''
        Return a copy of this vector with the '-' operator applied to each
        element
        '''
        return self.__class__(-self[0], -self[1], -self[2], -self[3])
    def __iadd__(self, v):
        '''
        Augmented assignment '+=' for adding a vector 'v' to this vector
        element-wise
        '''
        for i in range(0, len(self)):
            self[i] += v[i]
        return self
    def __isub__(self, v):
        '''
        Augmented assignment '-=' for subtracting a vector 'v' from
        this vector element-wise
        '''
        for i in range(0, len(self)):
            self[i] -= v[i]
        return self
    def __add__(self, v):
        '''
        Return the addition of this vector with the vector 'v'
        '''
        u = +self
        u += v
        return u
    def __sub__(self, v):
        '''
        Return the subtraction of vector 'v' from this vector
        '''
        u = +self
        u += -v
        return u
    def __invert__(self):
        '''
        Return the complex conjugate of this vector
        '''
        u = +self
        for i in range(0, len(u)):
            if isinstance(u[i], complex):
                u[i] = complex(u[i].real, -u[i].imag)
        return u
    def __imul__(self, v):
        '''
        Augmented assignment '*=' for multiplying this vector
        with a scalar, vector or matrix 'v'
        '''
        if isinstance(v, (int, float)):
            for i in range(0, len(self)):
                self[i] *= v
        elif isinstance(v, Vector):     
            u = +self                   #if v is a vetor, each element of this
            self = 0                    #vector is multiplied by the corresponding
            for i in range(0, len(u)):  #element of v, then summing over all elements;
                self += u[i]*v[i]       #effectively transposing vector 'v'
        elif isinstance(v, Matrix):
            u = []
            for i in range(0, 4):
                u.append(self[i])       #transposing vector
            u *= v
            self = u
        return self
    def __mul__(self, v):
        '''
        Return the multiplication of this vector with
        scalar, vector or matrix 'v'
        '''
        u = +self
        u *= v
        return u                      
    def __rmul__(self, v):
        '''
        Return the multiplication of this vector with
        scalar or vector 'v', no matter the order
        '''
        u = +self
        u *= v
        return u
    def __itruediv__(self, s):
        '''
        Augmented assignment '/=' for dividing this vector
        by a scalar 's'
        '''
        for i in range(0, len(self)):
            self[i] /= s
        return self
    def __truediv__(self, s):
        '''
        Return the division of this vector by a scalar 's'
        '''
        u = +self
        u /= s
        return u
    def __abs__(self):
        '''
        Return the Frobenius norm of this vector
        '''
        return (self*self)**(1/2)
    def __ipow__(self, p):
        '''
        Augmented assignment '**=' for implementing the
        absolute value of this vector to the power 'p'
        '''
        if p % 2 == 0:              #when p is even,
            a = abs(self) ** p      #only an integer is returned
        else:                           #when p is odd,
            a = (self ** (p-1))*self    #an integer multiplied
        return a                        #by this vector is returned
    def __pow__(self, p):
        '''
        Return the value of this vector to the power of
        a number 'p'
        '''
        u = +self
        u **= p
        return u
class Matrix:
    '''
    Class representing a matrix of size 4 x 4
    '''
    def __init__(self, m0j, m1j, m2j, m3j):
        '''
        Initialise the class with its four components
        '''
        self.m = [m0j, m1j, m2j, m3j]
    def __str__(self):
        '''
        Return a string to print for this matrix
        '''
        return "{}\n{}\n{}\n{}".format(self.m[0], self.m[1], self.m[2], self.m[3])
    def __repr__(self):
        '''
        Return the representation of this matrix
        '''
        return "{}({}, {}, {}, {})".format(self.__class__.__name__,
                                           self.m[0], self.m[1], self.m[2],
                                           self.m[3])
    def __getitem__(self, ind):
        '''
        Return the 'i, j' element of the matrix
        '''
        (i, j) = ind
        return self.m[i][j]
    def __setitem__(self, ind, s):
        '''
        Set the 'i, j' element of the matrix with the scalar s
        '''
        (i, j) = ind
        self.m[i][j] = s
        return self
    def __pos__(self):
        '''
        Return a copy of this matrix with the '+' operator applied to each
        element
        '''
        return self.__class__([+self[0, j] for j in range(0, 4)],
                              [+self[1, j] for j in range(0, 4)],
                              [+self[2, j] for j in range(0, 4)],
                              [+self[3, j] for j in range(0, 4)])
    def __neg__(self):
        '''
        Return a copy of this matrix with the '-' operator applied to each
        element
        '''
        return self.__class__([-self[0, j] for j in range(0, 4)],
                              [-self[1, j] for j in range(0, 4)],
                              [-self[2, j] for j in range(0, 4)],
                              [-self[3, j] for j in range(0, 4)])
    def __iadd__(self, m):
        '''
        Augmented assignment '+=' for adding a matrix 'm' to this matrix
        element-wise
        '''
        for i in range(0, 4):
            for j in range(0, 4):
                self[i, j] += m[i, j]
        return self
    def __isub__(self, m):
        '''
        Augmented assignment '-=' for subtracting a matrix 'm' from
        this matrix element-wise
        '''
        for i in range(0, 4):
            for j in range(0, 4):
                self[i, j] -= m[i, j]
        return self
    def __add__(self, m):
        '''
        Return the addition of this matrix with the matrix 'm'
        '''
        u = +self
        u += m
        return u
    def __sub__(self, m):
        '''
        Return the subtraction of matrix 'm' from this matrix
        '''
        u = +self
        u -= m
        return u
    def __invert__(self):
        '''
        Return the conjugate transpose for this matrix
        '''
        u = +self
        a = u - u   #creates a 4 x 4 matrix of zeros
        for i in range(0, 4):
            for j in range(0, 4):
                if isinstance(u[i, j], complex):
                    u[i, j] = complex(u[i, j].real, -u[i, j].imag)
                a[i, j] = u[j, i]
        return a
    def __imul__(self, m):
        '''
        Augmented assignment '*=' for multiplying this matrix
        with a scalar, vector (transposed or not) or matrix 'm'
        '''
        if isinstance(m, (int, float)):
            for i in range(0, 4):
                for j in range(0, 4):
                    self[i, j] *= m
        elif isinstance(m, Vector):
            u = +self
            self = m - m
            for i in range(0, 4):
                for j in range(0, 4):
                    self[i] += u[i, j]*m[j]
        elif isinstance(m, Matrix):
            u = +self
            self = self - self
            for i in range(0, 4):
                for j in range(0, 4):
                    for k in range(0, 4):
                        self[i, j] += u[i, k]*m[k, j]
        elif isinstance(m, list):       #list is transposed vector
            u = +self
            self = []
            for j in range(0, 4):
                a = 0
                for i in range(0, 4):
                    a += m[i]*u[i, j]
                self.append(a)
        return self
    def __mul__(self, m):
        '''
        Return the multiplication of this matrix with
        scalar, vector or matrix 'm'
        '''
        u = +self
        u *= m
        return u
    def __rmul__(self, m):
        '''
        Return the multiplication of this matrix with
        scalar or vector 'm', no matter the order
        '''
        u = +self
        u *= m
        return u
    def __itruediv__(self, s):
        '''
        Augmented assignment '/=' for dividing this matrix
        by a scalar 's'
        '''
        for i in range(0, 4):
            for j in range(0, 4):
                self[i, j] /= s
        return self
    def __truediv__(self, s):
        '''
        Return the division of this matrix by a scalar 's'
        '''
        u = +self
        u /= s
        return u
    def __abs__(self):
        '''
        Return the Frobenius norm of this Matrix
        '''
        conj = ~self
        u = 0
        for i in range(0, 4):
            for j in range(0, 4):
                u += conj[i, j]*self[j, i]
        return u**(1/2)
    def __ipow__(self, p):
        '''
        Augmented assignment '**=' for implementing the
        absolute value of this matrix to the power 'p'
        '''
        u = +self
        if p > 0:
            for n in range(p):
                self *= u
        elif p == 0:
            self = u - u        #creating 4 x 4 matrix of zeros
            for i in range(0, 4):
                for j in range(0, 4):
                    if i == j:
                        self[i, j] = 1
        return self
    def __pow__(self, p):
        '''
        Return the value of this matrix to the power of
        a number 'p'
        '''
        u = +self
        u **= p
        return u
class FourVector(Vector):
    '''
    Class representing a vector of length 4 which can be used to
    represent physical quantities (usually space-time or
    momentum-energy); referred to as a four-vector
    '''
    def __invert__(self):
        '''
        Lowers or raises the index of this four-vector
        '''
        u = +self
        m = Matrix([1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1])
        m *= u
        return m
    def __abs__(self):
        '''
        Return the norm of this four-vector
        '''
        return (self*~self)**(1/2)
class BoostMatrix(Matrix):
    '''
    Class representing a Lorentz Boost matrix
    '''
    def __init__(self, q_mu):
        '''
        Initialise the class with its four components
        '''
        B_x = q_mu[1]/q_mu[0]
        B_y = q_mu[2]/q_mu[0]
        B_z = q_mu[3]/q_mu[0]
        B = (B_x**2 + B_y**2 + B_z**2)**(1/2)
        gamma = 1/((1 - B**2)**(1/2))
        alpha = (gamma**2)/(1 + gamma)
        m0 = [gamma, -gamma*B_x, -gamma*B_y, -gamma*B_z]
        m1 = [-gamma*B_x, 1 + alpha*(B_x**2), alpha*B_x*B_y, alpha*B_x*B_z]
        m2 = [-gamma*B_y, alpha*B_y*B_x, 1 + alpha*(B_y**2), alpha*B_y*B_z]
        m3 = [-gamma*B_z, alpha*B_z*B_x, alpha*B_z*B_y, 1 + alpha*(B_z**2)]
        self.m = [m0, m1, m2, m3]
if __name__ == "__main__":
#    C = Vector(1, 2, 3, 4)
#    D = Vector(1+1j, 2j, 3+4j, 5)
#    E = Vector(4, -3, -2, 1)
#    print(C*E)
#    print(E*C)
#    print(abs(D))
#    print(+E)
#    A = Matrix([1, 2j, -3-4j, 4], [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4])
#    B = Matrix([4, 3, 2, 1], [4, 3, 2, 1], [4, 3, 2, 1], [4, 3, 2, 1])
#    C = Matrix([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1])
#    D = Vector(2, 1, 2, 1)
#    print(A**0)
     A = FourVector(1, 2, 3, 4)
     B = BoostMatrix(A)
     print(type(B))