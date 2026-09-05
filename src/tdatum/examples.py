"""Construct explicit matrix pairs; call ``t_datum()`` to validate a pair.

The constructors are adapted from Yuma Mizuno's local SageMath library.
They construct examples; they do not test finite type or periodicity.
"""

from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.matrix.all import *
from sage.functions.generalized import kronecker_delta
from sage.misc.flatten import flatten
from sage.arith.functions import lcm
from sage.arith.misc import GCD
from sage.combinat.root_system.dynkin_diagram import DynkinDiagram
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.misc.cachefunc import cached_method
from functools import reduce

__all__=['RSG','SG','TamelyLaced','UntwistedAffine','Unknown','LengthOne','Rank2']

class TDatumConstructor(object):
    def __init__(self,variable_name='z'):
        self._variable_name =variable_name
        
    def variable_name(self):
        return self._variable_name
    
    def indices(self):
        raise NotImplementedError()
    
    def size(self):
        return len(self.indices())
        
    def _N0_entry(u,v):
        raise NotImplementedError()
    
    def _N_plus_entry(self,u,v):
        raise NotImplementedError()
    
    def _N_minus_entry(self,u,v):
        raise NotImplementedError()
    
    def variable(self):
        return LaurentPolynomialRing(QQ,self.variable_name()).gen()
    
    def pair(self):
        N0=matrix(self.size() , lambda u,v: self._N0_entry(self.indices()[u],self.indices()[v]))
        Ap = N0-matrix(self.size() , lambda u,v: self._N_plus_entry(self.indices()[u],self.indices()[v]))
        Am = N0-matrix(self.size() , lambda u,v: self._N_minus_entry(self.indices()[u],self.indices()[v]))
        return (Ap,Am)

    def t_datum(self, D='identity'):
        """Return a validated T-datum, with an explicitly supplied D if needed."""
        from .t_datum import TDatum
        return TDatum(*self.pair(), D=D)
    
        
class RSG(TDatumConstructor):
    """Reduced sine-Gordon matrix pairs for an integer continued-fraction list."""

    def __init__(self,n_list,variable_name='z'):
        super(RSG,self).__init__(variable_name)
        self._n_list = tuple(n_list)
        if (not self._n_list or any(n not in ZZ or n < 1 for n in self._n_list)
                or self._n_list[0] < 2 or sum(self._n_list) <= 2):
            raise ValueError("Use a nonempty positive integer list with n1 >= 2 and sum(n_list) > 2.")
        
    def n(self,a):
        return self._n_list[a-1]
    def F(self):
        return len(self._n_list)
    def p(self,a):
        if a==1:
            return 1
        elif a==2:
            return self.n(1)
        else:
            return self.q(a-1)
    def q(self,a):
        if a==0:
            return 1
        elif a==1:
            return self.n(1)
        else:
            return self.n(a)*self.q(a-1)+ self.p(a-1)
        
    def epsilon(self,a):
        return (-1)**(a-1)
    
    @cached_method
    def indices(self):
        am = [ range(1,self.n(1)-1) ] + [ range(1,1+self.n(a)) for a in range(2,1+self.F()) ]
        return reduce(lambda a,b: a+b,[[ (a,m) for m in am[a-1] ] for a in range(1,1+self.F())])
    
    def _N0_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(b,k):
            result += z**(2*self.p(a)) + 1
        return result
    
    def _N_plus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if self.n(1)==2 and (a,m)==(2,1) and (b,k)==(2,1):
            result+=z**(self.p(a))
        if (a,m)==(2,1):
            if (b,k) == (2,2):
                result+= z**(self.p(a))
            elif b == 1:
                result+=  z**self.p(a)*(z**(k+1) + z**(-k-1))
        elif a>2 and m==1 and self.epsilon(a)==-1:
            if (b,k) ==(a,2):
                result +=  z**(self.p(a))
            elif (b,k) == (a-2 , self.n(a-2) - 2*kronecker_delta(a,3)):
                result +=  z**(self.p(a))
            elif b==a-1:
                result +=  z**(self.p(a))*(z**(self.p(a) -(self.n(a-1) +1 -k)*self.p(a-1) ) + z**(-self.p(a) +(self.n(a-1) +1 -k)*self.p(a-1) ))
        elif a>2 and m==1 and self.epsilon(a)==1:
            if self.n(a)==1 and (b,k)==(a+1,1) :
                result+= z**self.p(a)
        elif abs(self.indices().index((a,m)) - self.indices().index((b,k)))==1:
            if b%2==0:
                result+= z**self.p(a)            
        return result
    
    def _N_minus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(2,1):
            if (b,k) == (1,1):
                result+= z**self.p(a)
            elif self.n(2)==1 and (b,k)==(3,1):
                result+= z**self.p(a)
        elif a>2 and m==1 and self.epsilon(a)==1:
            if (b,k) ==(a,2):
                result+= z**self.p(a)
            elif (b,k) == (a-2 , self.n(a-2) - 2*kronecker_delta(a,3)):
                result+= z**self.p(a)
            elif b==a-1:
                result+= z**self.p(a)*(z**(self.p(a) -(self.n(a-1) +1 -k)*self.p(a-1) ) + z**(-self.p(a) +(self.n(a-1) +1 -k)*self.p(a-1) ))
        elif a>2 and m==1 and self.epsilon(a)==-1:
            if self.n(a)==1 and (b,k)==(a+1,1) :
                result+= z**self.p(a)
        elif abs(self.indices().index((a,m)) - self.indices().index((b,k)))==1:
            if b%2==1:
                result+= z**self.p(a)               
        return result
            
        
    
    
class SG(TDatumConstructor):
    def __init__(self,n_list,variable_name='z'):
        super(SG,self).__init__(variable_name)
        self._n_list = n_list
        
    def n(self,a):
        return self._n_list[a-1]
    def F(self):
        return len(self._n_list)
    def p(self,a):
        if a==1:
            return 1
        elif a==2:
            return self.n(1)
        else:
            return self.q(a-1)
    def q(self,a):
        if a==0:
            return 1
        elif a==1:
            return self.n(1)
        else:
            return self.n(a)*self.q(a-1)+ self.p(a-1)
        
    def epsilon(self,a):
        return (-1)**(a-1)
    
    @cached_method
    def indices(self):
        am = [[-2,-1,0]+list(range(1,self.n(1)-1))] + [ list(range(1,1+self.n(a))) for a in range(2,1+self.F())]
        return reduce(lambda a,b: a+b,[[ (a,m) for m in am[a-1] ] for a in range(1,1+self.F())])
    
    
    def _N0_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(b,k):
            result += z**(2*self.p(a)) + 1
        return result
    
    def _N_plus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(2,1):
            if (b,k) == (2,2) or (b,k) == (1,-1) or (b,k) == (1,-2):
                result+= z**(self.p(a))                
            elif b == 1:
                result+=  z**self.p(a)*(z**(k+1) + z**(-k-1))
        elif a>2 and m==1 and self.epsilon(a)==-1:
            if (b,k) ==(a,2):
                result +=  z**(self.p(a))
            elif (b,k) == (a-2 , self.n(a-2) - 2*kronecker_delta(a,3)):
                result +=  z**(self.p(a))
            elif b==a-1:
                result +=  z**(self.p(a))*(z**(self.p(a) -(self.n(a-1) +1 -k)*self.p(a-1) ) + z**(-self.p(a) +(self.n(a-1) +1 -k)*self.p(a-1) ))
        elif a>2 and m==1 and self.epsilon(a)==1:
            if self.n(a)==1 and (b,k)==(a+1,1) :
                result+= z**self.p(a)
        elif abs(self.indices().index((a,m)) - self.indices().index((b,k)))==1:
            if m>=0 and k>=0:
                if b%2==0:
                    result+= z**self.p(a)   
        return result
    
    def _N_minus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(2,1):
            if self.n(2)==1 and (b,k)==(3,1):
                result+= z**self.p(a)
        elif a>2 and m==1 and self.epsilon(a)==1:
            if (b,k) ==(a,2):
                result+= z**self.p(a)
            elif (b,k) == (a-2 , self.n(a-2) - 2*kronecker_delta(a,3)):
                result+= z**self.p(a)
            elif b==a-1:
                result+= z**self.p(a)*(z**(self.p(a) -(self.n(a-1) +1 -k)*self.p(a-1) ) + z**(-self.p(a) +(self.n(a-1) +1 -k)*self.p(a-1) ))
        elif a>2 and m==1 and self.epsilon(a)==-1:
            if self.n(a)==1 and (b,k)==(a+1,1) :
                result+= z**self.p(a)
        elif (m,k)==(-1,0) or (m,k)==(0,-1) or (m,k)==(-2,0) or (m,k)==(0,-2):
            result+=z**self.p(a)
        elif abs(self.indices().index((a,m)) - self.indices().index((b,k)))==1:
            if m>=0 and k>=0:
                if b%2==1:
                    result+= z**self.p(a)       
        return result
    
    
class UntwistedAffine(TDatumConstructor):
    def __init__(self,dynkin_type,rank,level,variable_name='z'):
        super(UntwistedAffine,self).__init__(variable_name)
        self._dynkin_type = dynkin_type
        self._rank=rank
        self._level=level
        
    def dynkin_type(self):
        return self._dynkin_type
    
    def rank(self):
        return self._rank
    
    def level(self):
        return self._level
    
    def cartan_matrix(self):
        return CartanMatrix([self.dynkin_type(),self.rank()])
    
    def d(self,a):
        C=self.cartan_matrix()
        return DynkinDiagram(C).symmetrizer().list()[a-1]
    
    def t_lcm(self):
        C=self.cartan_matrix()
        return lcm(DynkinDiagram(C).symmetrizer().list())
    
    def t(self,a):
        return self.t_lcm()/self.d(a)
    
    def _adjacency(self,a,b):
        if self.cartan_matrix()[a-1][b-1]<0:
            return 1
        else:
            return 0
    
    @cached_method
    def indices(self):
        return flatten([[(a,m) for m in range(1,self.t(a)*self.level())] for a in range(1,1+self.rank())],max_level=1)
    
    def _sinh_hat(self,c):
        z=self.variable()
        return (z**c - z**(-c))/2
    
    def _N0_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(b,k):
            result += z**(2*self.d(a)) + 1
        return result
    
    def _N_minus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        LPR=z.parent()
        ta,tb=self.t(a),self.t(b)
        if tb*m==ta*k:
            result += LPR(z**self.d(a)*(self._sinh_hat(self.t_lcm()/ta) /self._sinh_hat(self.t_lcm()/max(ta,tb))))
        if tb>ta:
            for j in range(1,tb-ta+1):
                if tb/ta*(m+1)-j==k or tb/ta*(m-1)+j==k:
                    result+=LPR(z**self.d(a)*(self._sinh_hat(j*self.t_lcm()/tb) /self._sinh_hat(self.t_lcm()/tb)))
        return result*self._adjacency(a,b)
    
    def _N_plus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if a==b and abs(m-k)==1:
            result+=z**self.d(a)
        return result
    
    
class TamelyLaced(TDatumConstructor):
    def __init__(self,cartan_matrix,symmetrizer,level,variable_name='z'):
        super(TamelyLaced,self).__init__(variable_name)
        self._cartan_matrix = cartan_matrix
        self._symmetrizer = symmetrizer
        self._rank=cartan_matrix.ncols()
        self._level=level
        z=self.variable()
        if not (diagonal_matrix(self._rank,[dd for dd in self._symmetrizer]) * self._cartan_matrix).is_symmetric():
            raise ValueError("The input matrix is not a Cartan matrix with a symmetrizer {}.".format(symmetrizer))
      
    def cartan_matrix(self):
        return self._cartan_matrix
    def _is_adjacent(self,a,b):
        return self.cartan_matrix()[a][b]<0
    def rank(self):
        return self._rank
    def level(self):
        return self._level
    def d(self,a):
        return self._symmetrizer[a]
    def t_lcm(self):
        return lcm([self.d(a) for a in range(self.rank())])
    def t(self,a):
        return self.t_lcm()/self.d(a)
    
    def _N0_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        if a==b and m==k:
            return 1+z**(2*self.d(a))
        else:
            return 0
    
    def _N_minus_entry2(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        if self._is_adjacent(a,b):
            if self.d(a)>1:
                p=self.d(a)/self.d(b)
                j=k-p*m
                if -p<j<p:
                    return sum(z**(self.d(a) +p-abs(j)+1-2*kk ) for kk in range(1,p-abs(j)+1))
                else:
                    return 0
            elif self.d(a)==1:
                if k==m/self.d(b):
                    return z**(self.d(a))
                else:
                    return 0
        else:
            return 0
        
    def z_integer(self,n,degree):
        z=self.variable()
        if n<0:
            raise ValueError("n should be a non-negative integer.")
        return sum(z**((degree*(n-(2*k-1)))) for k in range(1,n+1))
    
    def _N_minus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        C=self.cartan_matrix()
        if self._is_adjacent(a,b):
            dd=lcm(self.d(a),self.d(b))
            tta=dd/self.d(a)
            ttb=dd/self.d(b)
            gcd_d=GCD(self.d(a),self.d(b))
            if m*ttb/tta in ZZ:
                p=m*ttb/tta
                j=p-k
                if -ttb<j<ttb:
                    f=-C[a][b]*z**(self.d(a))*self.z_integer(ttb-abs(j),self.d(b))/tta
                    return sum(c*z**(e) for c,e in zip(f.coefficients(),f.exponents()))
                else:
                    return 0
            else:
                return 0
        else:
            return 0
        
    def _N_plus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        if abs(m-k)==1 and a==b:
            return z**(self.d(a))
        else:
            return 0
    
    @cached_method
    def indices(self):
        return flatten([[(a,m) for m in range(1,self.t(a)*self.level())] for a in range(self.rank())],max_level=1)

class Unknown(TDatumConstructor):
    def __init__(self,quiver_type,variable_name='z'):
        super(Unknown,self).__init__(variable_name)
        self._quiver_type=quiver_type

    def indices(self):
        if self._quiver_type == 'E6':
            return [0, 1]
        elif self._quiver_type == 'E7':
            return [0, 1, 2]
        elif self._quiver_type == 'E6^11':
            return [0, 1, 2]
        
    def _triple(self,sign):
        z=self.variable()
        if self._quiver_type=='E6':
            if sign==1:
                return matrix(2,2,[0,z,z+z**5+z**9,0])
            elif sign==-1:
                return matrix(2,2,[0,0,z**3+z**7,0])
            elif sign==0:
                return diagonal_matrix([1+z**kk for kk in [2,10]])
            
        if self._quiver_type=='E7':
            if sign==0:
                return diagonal_matrix([1+z**kk for kk in [4,4,6]])
            if sign==1:
                return matrix(3,3,[0,z**2,0,z**2,0,z**2,0,z**2+z**4,0])
            elif sign==-1:
                return diagonal_matrix([z**2,z**2,0])
            
        if self._quiver_type=='E6^11':
            if sign == 0:
                return diagonal_matrix([1+z**kk for kk in [2,6,8]])
            if sign==1:
                return matrix(3,3,[0,0,z,z**3,0,0,z+z**7,z**2+z**6,0])
            elif sign==-1:
                return matrix(3,3,[0,z,0,z+z**5,0,0,0,0,0])
    def _N0_entry(self, u, v):
        return self._triple(0)[u][v]
    def _N_plus_entry(self, u, v):
        return self._triple(1)[u][v]
    def _N_minus_entry(self, u, v):
        return self._triple(-1)[u][v]
            
            
class LengthOne(TDatumConstructor):
    def __init__(self,degree_list,variable_name='z'):
        if degree_list[::-1]!=degree_list:
            raise ValueError('{} is not a palindrome.'.format(degree_list))
        super(LengthOne,self).__init__(variable_name)
        self._degree=len(degree_list)+1
        self._degree_list=degree_list
    def indices(self):
        return [0]
    def _N0_entry(self,u,v):
        z=self.variable()
        return 1+z**self._degree
    def _N_plus_entry(self,u,v):
        z=self.variable()
        return sum(c*z**(degree+1) for degree,c in enumerate(self._degree_list) if c>0)
    def _N_minus_entry(self,u,v):
        z=self.variable()
        return sum((-c)*z**(degree+1) for degree,c in enumerate(self._degree_list) if c<0)
    
class Rank2(TDatumConstructor):
    """Six stored rank-two examples, selected by labels 1 through 6."""

    def __init__(self, label, variable_name='z'):
        super(Rank2, self).__init__(variable_name)
        if label not in ZZ or not 1 <= label <= 6:
            raise ValueError("Rank2 labels are the integers 1,...,6.")
        self._label = label
        
    def indices(self):
        return [0, 1]
    
    def _N0_entry(self, u, v):
        return self._triple(0)[u][v]
    def _N_plus_entry(self, u, v):
        return self._triple(1)[u][v]
    def _N_minus_entry(self, u, v):
        return self._triple(-1)[u][v]
    
    def _triple(self, sign):
        z=self.variable()
        label = self._label
        if label == 1:
            if sign == 0:
                return diagonal_matrix(2, [1+z**2, 1+z**2])
            if sign == 1:
                return matrix(2, 2, [0, z, z, 0])
            if sign == -1:
                return zero_matrix(2)
        if label == 2:
            if sign == 0:
                return diagonal_matrix(2, [1+z**2, 1+z**6])
            if sign == 1:
                return matrix(2, 2, [0, z, z+z**5, 0])
            if sign == -1:
                return matrix(2, 2, [0, 0, z**3, 0])
        if label == 3:
            if sign == 0:
                return diagonal_matrix(2, [1+z**2, 1+z**10])
            if sign == 1:
                return matrix(2, 2, [0, z, z+z**5+z**9, 0])
            if sign == -1:
                return matrix(2, 2, [0, 0, z**3+z**7, 0])
        if label == 4:
            if sign == 0:
                return diagonal_matrix(2, [1+z**2, 1+z**2])
            if sign == 1:
                return matrix(2, 2, [0, z, z, 0])
            if sign == -1:
                return diagonal_matrix(2, [z, z])
        if label == 5:
            if sign == 0:
                return diagonal_matrix(2, [1+z**2, 1+z**3])
            if sign == 1:
                return matrix(2, 2, [0, z, z+z**2, 0])
            if sign == -1:
                return diagonal_matrix(2, [z, 0])
        if label == 6:
            if sign == 0:
                return diagonal_matrix(2, [1+z**2, 1+z**2])
            if sign == 1:
                return matrix(2, 2, [0, z, z, z])
            if sign == -1:
                return zero_matrix(2)
