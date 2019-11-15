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
                self.v[0], self.v[1], self.v[2], self.v[3])
    
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
    
    def __imul__(self, s):
        '''
        Augmented assignment '*=' for multiplying this vector
        with a scalar 's'
        '''
        for i in range(0, len(self)):
             self[i] *= s
        return self
    
    def __mul__(self, v):
        '''
        Return the multiplication of this vector with
        vector 'v'
        '''
        u = +self
        if type(v) == int or type(v) == float:
            u *= v          #if v is a number, execute __imul__ method
            return u
        else:
            a = 0                       #if v is a vetor, each element of this
            for i in range(0, len(u)):  #vector is multiplied by the corresponding
                a += u[i]*v[i]          #element of v, then summing over all elements;
            return a                    #effectively transposing vector 'v'
                
    
    def __rmul__(self, v):
        '''
        Return the multiplication of this vector with
        vector 'v', no matter the order
        '''
        u = +self
        if type(v) == int or type(v) == float:
            u *= v
            return u
        else:
            a = u*v
            return a
        
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
            a = abs(self) ** p      #only an integer
            return a                #is returned
        else:
            u = +self               #when p is odd,
            a = self ** (p-1)       #an integer multiplied
            return a*u              #by this vector is returned
        
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
                self.m[0], self.m[1], self.m[2], self.m[3])
        
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
        return self.__class__(+self[0], +self[1], +self[2], +self[3])
    
    def __neg__(self):
        '''
        Return a copy of this matrix with the '-' operator applied to each
        element
        '''
        return self.__class__(-self[0], -self[1], -self[2], -self[3])
    
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
    
  
if __name__ == "__main__":
    a = Matrix([1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4])
    b = Matrix([4, 3, 2, 1], [4, 3, 2, 1], [4, 3, 2, 1], [4, 3, 2, 1])
    print(a + b)
    
    
    
    
    
    