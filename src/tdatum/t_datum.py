r"""
T-datum

AUTHORS:

- Yuma Mizuno (2020): initial version

"""

# ****************************************************************************
#       Copyright (C) 2020 Yuma Mizuno mizuno.y.aj@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# u should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ****************************************************************************


from sage.structure.sage_object import SageObject
from sage.structure.element import Matrix
from sage.misc.cachefunc import cached_method
from sage.rings.all import QQ,ZZ,LaurentPolynomialRing
from sage.arith.all import lcm,gcd
from sage.matrix.all import identity_matrix,diagonal_matrix,matrix,MatrixSpace
from sage.graphs.digraph import DiGraph
from sage.graphs.graph_generators import GraphGenerators
from sage.misc.misc_c import prod
from sage.misc.flatten import flatten
from sage.misc.latex import latex
from sage.functions.generalized import kronecker_delta
from sage.combinat.permutation import Permutation,Permutations
from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
from sage.rings.infinity import PlusInfinity
from sage.plot.colors import rainbow
    
from itertools import product
from functools import reduce

__all__=['TDatum','MutationLoop']

class TDatum(SageObject):
    r"""
    A T-datum.

    A T-datum of size $r$ is a triple $(A_+,A_-,D)$ where
    $A_+$ and $A_-$ are matrices in $\mathrm{Mat}_{r \times r}(\ZZ[z^{\pm}])$ and
    $D$ is a positive integer matrix.
    These matrices satisfy 
    $N_0 D= D N_0$, $D^{-1}N_{\pm}D \in $\mathrm{Mat}_{r \times r}(\ZZ[z^{\pm}])$,
    and $A_+ D A_-^{\dagger} = A_- D A_+^{\dagger}$,
    where $A_{\pm}^{\dagger} := (A_{\pm}|_{z=z^{-1}})^{\mathsf{T}}$.
    The condition $A_+ D A_-^{\dagger} = A_- D A_+^{\dagger}$ is called the symplectic relation.
    Moreover, the pair of matrices $(A_+,A_-)$ must be written as
    $A_\pm = N_0 - N_{\pm}$ using a triple $(N_0,N_+,N_-) \in \mathrm{Mat}_{r \times r}(\ZZ[z^{\pm}])^3$
    satisfying the following conditions:

    - (N1) $n_{ab;p}^{0} = \delta_{ab} \delta_{p0} + \delta_{a\sigma(b)} \delta_{pp_a}$ for some $\sigma \in \mathfrak{S}_r$ and $p_a \in \ZZ_{>0}$,

    - (N2) $n_{ab;p}^{+} \geq 0$ and $n_{ab;p}^{-} \geq 0$ for any $a,b,p$,

    - (N3) $n_{ab;p}^{+} = 0$ and $n_{ab;p}^{-}=0$ unless $0<p<p_a$,
    
    - (N4) $n_{ab;p}^{+} n_{ab;p}^{-} = 0$ for any $a,b,p$,
    
    where we write the entries of $N_{\varepsilon}$ as 
    $N_{\varepsilon}=\bigl(\sum_{p \in \ZZ_{\geq 0}} n_{ab;p}^{\varepsilon} z^p \bigr)_{1 \leq a,b \leq r}$.

    INPUT:
    
    - ``A_plus`` -- square matrix whose entries are univariate Laurent polynomial

    - ``A_minus`` -- square matrix whose entries are univariate Laurent polynomial

    - ``D`` -- positive integer diagonal matrix
    
    OUTPUT:
    
    - the T-datum for ``A_plus``, ``A_minus``, and ``D``
    
    EXAMPLES:
    
    An example of T-datum when ``D`` is the identity matrix::

        sage: z = LaurentPolynomialRing(QQ,'z').gen()
        sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
        sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
        sage: td = TDatum(A_plus,A_minus)
        sage: N0,Np,Nm = td.triple()
        sage: N0
        [1 + z^2       0]
        [      0 1 + z^2]
        sage: Np
        [0 0]
        [0 0]
        sage: Nm
        [0 z]
        [z 0]

    An example of T-datum when ``D`` is not the identity matrix::

        sage: z = LaurentPolynomialRing(QQ,'z').gen()
        sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
        sage: A_minus = matrix(2,2,[1+z^2,-z,-2*z,1+z^2])
        sage: td = TDatum(A_plus,A_minus)
        sage: N0,Np,Nm = td.triple()
        sage: N0
        [1 + z^2       0]
        [      0 1 + z^2]
        sage: Np
        [0 0]
        [0 0]
        sage: Nm
        [  0   z]
        [2*z   0]
    """
    
    def __init__(self, A_plus, A_minus, D='identity', check=True):
        if check:
            if not isinstance(A_plus,Matrix):
                raise ValueError("The input should be a pair of matrices.")
            if not isinstance(A_minus,Matrix):
                raise ValueError("The input should be a pair of matrices.")
            from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing_univariate
            if not isinstance(A_plus.base_ring(),LaurentPolynomialRing_univariate):
                raise ValueError("The base ring of A_plus and A_minus should be a univariate laurent polynomial ring.")
            if not isinstance(A_minus.base_ring(),LaurentPolynomialRing_univariate):
                raise ValueError("The base ring of A_plus and A_minus should be a univariate laurent polynomial ring.")
        
        self._size = A_plus.ncols()
        if self._size != A_plus.nrows() or self._size != A_minus.ncols() or self._size != A_minus.nrows():
             raise ValueError("The inputs should be square matrices.")
        
        lp = A_plus.base_ring()
        self._A_plus = MatrixSpace(lp,A_plus.nrows(),A_plus.ncols())(A_plus)
        self._A_minus = MatrixSpace(lp,A_minus.nrows(),A_minus.ncols())(A_minus)
        self._triple = None
        if D=='identity':
            self._D = MatrixSpace(ZZ,self._size)(identity_matrix(self._size))
        else:
            self._D = MatrixSpace(ZZ,self._size)(D)
        self._var_name = A_plus.base_ring().variable_name()
        
        N0 = self.triple()[0]
        z = self.variable()
        P=diagonal_matrix([z**(-d) for d in self.degrees()])*(N0-identity_matrix(N0.ncols()))
        sigma = (P.is_permutation_of(identity_matrix(N0.ncols()),check = True))[1]
        self._permutation = Permutations(N0.ncols())(prod(sigma))
        
        if check:
            if not self._is_t_datum():
                raise ValueError("The inputs is not a T-datum.")
    
    def symmetrizer(self):
        """
        Return the symmetrizer tuple of the T-datum.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: Ap = matrix(2, 2, [1+z^2,-z,-2*z,1+z^2])
            sage: Am = matrix(2, 2, [1+z^2,0,0,1+z^2])
            sage: D = diagonal_matrix([1,2])
            sage: td = TDatum(Ap, Am, D)
            sage: td.symmetrizer()
            (1, 2)
        """
        return tuple(self._D.diagonal())
    
    def _is_t_datum(self):
        """
        Return ``True`` if the input is a T-datum and ``False`` otherwise.
        """
        z = self.variable()
        Ap,Am = self.pair()
        N0,Np,Nm = self.triple()
        #condition (C-2)
        if not all([all(flatten([[coeff>0 for coeff in f.coefficients()] for f in X.list() if f!=0])) for X in self.triple()]):
            raise ValueError('The inputs do not satisfy (N2).')
        #condition (C-4)
        for a in range(Nm.nrows()):
            for b in range(Nm.ncols()):
                fp=Np[a][b]
                fm=Nm[a][b]
                for coeff_p,exponent_p in zip(fp.coefficients(),fp.exponents()):
                    for coeff_m,exponent_m in zip(fm.coefficients(),fm.exponents()):
                        if exponent_p==exponent_m:
                            if coeff_p*coeff_m != 0:
                                raise ValueError('The inputs do not satisfy (N4).')
        #condition (C-3)
        l=self.degrees()
        for a in range(Nm.nrows()):
            for b in range(Nm.ncols()):
                fp=Np[a][b]
                fm=Nm[a][b]
                if not all([0<exp and exp<l[a] for exp in fp.exponents()]):
                    raise ValueError('The inputs do not satisfy (N3).')
                if not all([0<exp and exp<l[a] for exp in fm.exponents()]):
                    raise ValueError('The inputs do not satisfy (N3).')
        #condition (C-1)
        P=diagonal_matrix([z**(-d) for d in self.degrees()])*(N0-identity_matrix(N0.ncols()))
        if not P.is_permutation_of(identity_matrix(N0.ncols())):
            raise ValueError('The inputs do not satisfy (N1).')

        #symplectic relation
        D=diagonal_matrix(self.symmetrizer())
        if Ap*D*(Am.transpose().subs({z:1/z})) != Am*D*(Ap.transpose().subs({z:1/z})):
            raise ValueError('The inputs do not satisfy the symplectic relation.')
        
        return True
    
    def variable(self):
        """
        Return the variables of the polynomials in the matrices in the T-datum.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.variable()
            z
        """
        return LaurentPolynomialRing(QQ,self._var_name).gen()
    
    def size(self):
        """
        Return the size of the matrices in the T-datum.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.size()
            2
        """
        return self.triple()[0].ncols()

    def permutation(self):
        r"""
        Return the permutation $\sigma$ such that
        $n_{ab;p}^{0} = \delta_{ab} \delta_{p0} + \delta_{a\sigma(b)} \delta_{pp_a}$.

        EXAMPLES::
        
            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[[1-2*z^2,-2*z+z^3],[z^3,1]])
            sage: A_minus = matrix(2,2,[[1,z^3],[-2*z+z^3,1-2*z^2]])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.permutation()
            [2, 1]
        """
        return self._permutation
    
    def pair(self):
        """
        Return the pair of the matrices $(A_+, A_-)$.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.pair()
            (
            [1 + z^2       0]  [1 + z^2      -z]
            [      0 1 + z^2], [     -z 1 + z^2]
            )
        """
        return (self._A_plus, self._A_minus)
    
    def triple(self):
        r"""
        Return the triple of the matrices $(N_0,N_+,N_-)$.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.triple()
            (
            [1 + z^2       0]  [0 0]  [0 z]
            [      0 1 + z^2], [0 0], [z 0]
            )
        """
        if self._triple:
            return self._triple

        def polynomial_positive_part(polynomial):
            if polynomial in ZZ:
                return max(0,polynomial)
            else:
                x=polynomial.variables()[0]
                f=polynomial
                return sum(max(0,c)*x**e for e,c in zip(f.exponents(),f.coefficients()))

        def matrix_positive_part(mat):
            nr=mat.nrows()
            nc=mat.ncols()
            return matrix(nr,nc,lambda a,b:polynomial_positive_part(mat[a][b]))
        
        Ap, Am = self.pair()
        self._triple = (matrix_positive_part(Ap), matrix_positive_part(-Ap), matrix_positive_part(-Am))
        return self._triple
    
    def is_indecomposable(self):
        r"""
        Return ``True`` if $(A_+,A_-,D)$ is indecomposable and ``False`` otherwise.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.is_indecomposable()
            True

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2-z,0,0,1+z^2-z])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.is_indecomposable()
            False
        """
        z = self.variable()
        N0,Np,Nm = [ N.subs({z:1}) for N in self.triple()]
        return DiGraph(N0+Np+Nm).is_connected()

    def _n(self,a,b,p,sign):
        r"""
        Return the $(a,b;p)$-th entry of the matrices $N_{\varepsilon}$ where $\varepsilon=$``sign``.
        """
        if sign==0:
            N=self.triple()[0]
        elif sign==1:
            N=self.triple()[1]
        elif sign==-1:
            N=self.triple()[2]
        else:
            raise ValueError("{} should be a sign".format(sign))
        f=N[a][b]
        ec_list=[(e,c) for e,c in zip(f.exponents(),f.coefficients())]
        for e,c in ec_list:
            if e==p:
                return c
        return 0

    def _n_check(self,a,b,p,sign):
        r"""
        Return the $(a,b;p)$-th entry of the matrices $N_{\varepsilon}^{\vee}$ where $\varepsilon=$ ``sign``.
        """
        dd=self.symmetrizer()
        return self._n(a,b,p,sign)/dd[a]*dd[b]
    
    def degrees(self):
        r"""
        Return the tuple $(p_a)$ such that
        $n_{ab;p}^{0} = \delta_{ab} \delta_{p0} + \delta_{a\sigma(b)} \delta_{pp_a}$.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,-z,-z-z^5,1+z^6])
            sage: A_minus = matrix(2,2,[1+z^2,0,-z^3,1+z^6])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.degrees()
            (2, 6)
        """
        N0=self.triple()[0]
        ncol=N0.ncols()
        return tuple((sum(N0[k][i] for i in range(ncol))-1).exponents()[0] for k in range(ncol))
    
    def langlands_dual(self):
        r"""
        Return the Langlands dual T-datum $(A_+^{\vee},A_-^{\vee},D^{\vee})$.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-2*z,1+z^2])
            sage: D=diagonal_matrix([1,2])
            sage: td = TDatum(A_plus,A_minus,D)
            sage: td.triple()
            (
            [1 + z^2       0]  [0 0]  [  0   z]
            [      0 1 + z^2], [0 0], [2*z   0]
            )
            sage: td.langlands_dual().triple()
            (
            [1 + z^2       0]  [0 0]  [  0 2*z]
            [      0 1 + z^2], [0 0], [  z   0]
            )
        """
        D=diagonal_matrix(self.symmetrizer())
        delta=lcm(self.symmetrizer())*gcd(self.symmetrizer())
        Ap,Am = self.pair()
        D_dual = D.parent()(diagonal_matrix([ delta/d for d in self.symmetrizer() ]))
        Ap_dual = Ap.parent()(D.inverse()*Ap*D)
        Am_dual = Am.parent()(D.inverse()*Am*D)
        return TDatum(Ap_dual,Am_dual,D_dual)
        
    def sign_dual(self):
        r"""
        Return the sign dual T-datum $(A_-,A_+,D)$.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.triple()
            (
            [1 + z^2       0]  [0 0]  [0 z]
            [      0 1 + z^2], [0 0], [z 0]
            )
            sage: td.sign_dual().triple()
            (
            [1 + z^2       0]  [0 z]  [0 0]
            [      0 1 + z^2], [z 0], [0 0]
            )
        """
        Ap, Am = self.pair()
        return TDatum(Am, Ap)

    def _phi(self,u,r):
        r"""
        Return the tuple $\varphi_u(r)$ defined by
        $\varphi_u(a,u+p) = (a,u+p_a)$ if $p=0$ and $\varphi_u(a,u+p) = (u,u+p)$ if $0<p<p_a$,
        where $r = (a,u+p)$.
        """
        a,up = r
        p = up - u
        sigma = self.permutation()
        if p ==0:
            return (sigma(a+1)-1,u+self.degrees()[sigma(a+1)-1])
        elif 0<p<self.degrees()[a]:
            return r
        else:
            raise ValueError('({}, {}) is not a valid input.'.format(u,r))
    
    def _phi_inv(self,u,r):
        r"""
        The inverse function of ``self._phi()``. See :func:`_phi<sage.combinat.t_datum.TDatum._phi>`.
        """
        a,up = r
        p = up - u
        sigma = self.permutation()
        if p ==self.degrees()[a]:
            return (sigma.inverse()(a+1)-1,u)
        elif 0<p<self.degrees()[a]:
            return r
        else:
            raise ValueError('({}, {}) is not a valid input.'.format(u,r))

    def _phi_inv_vec(self,u,r):
        r"""
        The inverse function of the composition $\varphi_{t-1} \circ \dots \circ \varphi_{0}$.
        See :func:`_phi<sage.combinat.t_datum.TDatum._phi>`.
        """
        a,up = r
        if u==0:
            return r
        elif u>0:
            return self._phi_inv_vec(u-1, self._phi_inv(u-1,r))
        else:
            return self._phi_inv_vec(u-1 ,self._phi(u),r)

    def _psi(sekf,r,t=1):
        r"""
        Return the tuple $(a,u+t)$ where $r = (a,u)$.
        """
        a,u = r
        return (a,u+t)

    def _R(self,R,u):
        r"""
        Return the tuple $\varphi_{t-1} \circ \dots \circ \varphi_{0}(R)$.
        """
        if u == 0:
            return R
        elif u > 0:
            return tuple(sorted([self._phi(u-1,r) for r in self._R(R,u-1)]))
        elif u < 0:
            return tuple(sorted([self._phi_inv(u+1,r) for r in self._R(R,u+1)]))
        else:
            raise ValueError('{} should be an integer.'.format(u))

    def maximal_initial_indices(self):
        r"""
        Return the tuple $((a,p) \in [1,r] \times \ZZ \mid 0 \leq p <p_a )$.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,-z,-z-z^5,1+z^6])
            sage: A_minus = matrix(2,2,[1+z^2,0,-z^3,1+z^6])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.maximal_initial_indices()
            ((0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5))
        """
        pp=self.degrees()
        return tuple(flatten([[(a,p) for p in range(pp[a])] for a in range(self.size())],max_level=1))

    def connected_component(self,a):
        r"""
        Return the indices of the connected component of the (valued) quiver associated with $(A_+,A_-,D)$
        that contains the index $(a,0)$.
        This method is defined only for indecomposable T-data.

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.connected_component(0)
            ((0, 0), (1, 1))
            sage: td.connected_component(1)
            ((0, 1), (1, 0))
        """
        if not self.is_indecomposable():
            raise ValueError('The self should be indecomposable.')
        ml_maximal = self.mutation_loop(R=None)
        B = ml_maximal.initial_b_matrix()
        R = self.maximal_initial_indices()
        a0 = R.index((a,0))
        R2 = tuple(R[i] for i in flatten([c for c in DiGraph(B).connected_components() if a0 in c]))
        return R2

    def mutation_loop(self,R=None):
        r"""
        Return the mutation loop associated with the T-datum $(A_+,A_-,D)$.

        INPUT:

            -``R`` -- tuple (default: ``None``); the tuple consisting of the indices of 
            the initial (valued) quiver in the mutation loop.
            If  ``R`` = ``None``, the tuple ``R`` is given by
            :func:`self.maximal_initial_indices()<sage.combinat.t_datum.TDatum.maximal_initial_indices()>`.

        OUTPUT: the mutation loop

        EXAMPLES::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(2,2,[1+z^2,0,0,1+z^2])
            sage: A_minus = matrix(2,2,[1+z^2,-z,-z,1+z^2])
            sage: td = TDatum(A_plus,A_minus)
            sage: gamma = td.mutation_loop()
            sage: gamma.initial_b_matrix()
            [ 0  0  0 -1]
            [ 0  0  1  0]
            [ 0 -1  0  0]
            [ 1  0  0  0]
            sage: gamma0 = td.mutation_loop(td.connected_component(1))
            sage: gamma0.initial_b_matrix()
            [ 0  1]
            [-1  0]
        """
        if not R:
            R=self.maximal_initial_indices()
        I=self._R(R,0)
        def ent(i,j):
            a,p=I[i]
            b,q=I[j]
            s1=-self._n(a,b,p-q,1)+self._n(a,b,p-q,-1)+self._n_check(b,a,q-p,1)-self._n_check(b,a,q-p,-1)
            s2=sum(sum(self._n(a,c,p-u,-1)*self._n_check(b,c,q-u,1)-self._n(a,c,p-u,1)*self._n_check(b,c,q-u,-1)\
                       for u in range(min(p,q))) for c in range(self.size()))
            return s1-s2
        B=matrix(len(I),lambda i,j: ent(i,j) )
        t = ZZ(len(self.maximal_initial_indices()) / len(R))
        i=[]
        for u in range(t):
            i.append([I.index(self._phi_inv_vec(u,(a,u1))) for a,u1 in list( self._R(R,u)) if u1==u])
        nu=Permutation([I.index(self._phi_inv_vec(t,self._psi(r,t)))+1 for r in I])
        return MutationLoop(B,i,nu)
    
    def plot_mutation_loop(self):
        r"""
        Plot the mutation loop assiciated with $(A_+,A_-,D)$.

        EXAMPLE::

            sage: z = LaurentPolynomialRing(QQ,'z').gen()
            sage: A_plus = matrix(1,1,[1-2*z**2+z**4])
            sage: A_minus = matrix(1,1,[1-z-z**3+z**4])
            sage: td = TDatum(A_plus,A_minus)
            sage: td.plot_mutation_loop()
            Digraph on 4 vertices

        .. PLOT::

            z = LaurentPolynomialRing(QQ,'z').gen()
            A_plus = matrix(1,1,[1-2*z**2+z**4])
            A_minus = matrix(1,1,[1-z-z**3+z**4])
            td = TDatum(A_plus,A_minus)
            p = td.plot_mutation_loop()
            sphinx_plot(p)
        """
        B = self.mutation_loop().initial_b_matrix()
        I = self.maximal_initial_indices()
        indices={}
        for i,v in enumerate(I):
            indices.update({i:v})
        cq = ClusterQuiver(B).relabel(indices)

        def plot_quiver(cluster_quiver, greens):
            n, m = cluster_quiver._n, cluster_quiver._m
            nlist = list(cluster_quiver._vertex_dictionary.values())[:n]
            mlist = list(cluster_quiver._vertex_dictionary.values())[n:]
            colors = rainbow(11)
            color_dict = { colors[0]:[], colors[1]:[], colors[6]:[], colors[5]:[] }

            dg = DiGraph( cluster_quiver._digraph )

            # For each edge in our graph we assign a color
            for edge in dg.edges():
                v1,v2,(a,b) = edge

                if v1 in nlist and v2 in nlist:
                    if (a,b) == (1,-1):
                        color_dict[ colors[0] ].append((v1,v2))
                    else:
                        color_dict[ colors[6] ].append((v1,v2))
                else:
                    if (a,b) == (1,-1):
                        color_dict[ colors[1] ].append((v1,v2))
                    else:
                        color_dict[ colors[5] ].append((v1,v2))
                if a == -b:
                    if a == 1:
                        dg.set_edge_label(v1, v2, '')
                    else:
                        dg.set_edge_label(v1, v2, a)

            # Partition out the green vertices
            for i in greens:
                if i in nlist:
                    nlist.remove(i)
                else:
                    mlist.remove(i)
            partition = (nlist, mlist, greens)

            vertex_color_dict = {'tomato': partition[0],
                                'lightblue': partition[1],
                                'lightgreen': partition[2]}

            options = {
                'graph_border' : True,
                'edge_colors': color_dict,
                'vertex_colors': vertex_color_dict,
                'edge_labels' : True,
                'vertex_labels': True,
                'vertex_size': 560
            }

            return dg.plot( **options )
        return plot_quiver(cq,[(a,p) for (a,p) in I if p == 0])
    
class MutationLoop(SageObject):
    r"""
    A mutation loop.

    A mutation loop is a tuple $\gamma = (B,d,\mathbf{i},\nu)$ where

    - $B$ is an $I \times I$ integer skew-symmetrizable matrix where $I=[1,\dots,n]$,

    - $d$ is a right symmetrizer of $B$,

    - $\mathbf{i}$ is a sequence $\mathbf{i} = (\mathbf{i}(0),\dots, \mathbf{i}(t-1))$ where each $\mathbf{i}(u)$ is a sequence of $I$,

    - $\nu$ is a permutation of $I$ such that $\nu(B(t)) = B$ and $\nu(d)=d$.

    From these data, we have the following sequence of matrix mutations:

    .. MATH::

        B =: B(0) \mapsto
        B(1)  \mapsto
        \cdots
        \mapsto
        B(t), \quad
        \mu_{\mathbf{i}(u)}: B(u) \mapsto B(u+1).

    We impose $B_{ij}(u)=0$ for any $i,j \in \mathbf{i}(u)$ so that the sequence of mutations
    $\mu_{\mathbf{i}(u)}$ is independent of the order of $\mathbf{i}(u)$.
    
    INPUT:
    
    - ``b_matrix`` -- a skew-symmetrizable matrix $B$,

    - ``sequence_of_indices`` -- a sequence of indices $\mathbf{i}$

    - ``permutation`` -- a permutation $\nu$

    - ``symmetrizer`` -- a right symmetrizer $d$ (default: the tuple $(1, \dots, 1)$)
    
    OUTPUT:
    
    - a mutation loop
    
    EXAMPLES::
    
        sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
        sage: i = [[0,2],[1]]
        sage: nu = Permutation([1,2,3])
        sage: gamma = MutationLoop(b,i,nu)

        sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-2,0])
        sage: d = [1,1,2]
        sage: i = [[0,2],[1]]
        sage: nu = Permutation([1,2,3])
        sage: gamma = MutationLoop(b,i,nu,d)
    
    """
    def __init__(self, b_matrix, sequence_of_indices, permutation,symmetrizer=None):
        self._b_matrix=b_matrix
        self._vertices=sequence_of_indices[:]
        self._permutation=permutation
        size=self._b_matrix.ncols()
        if not symmetrizer:
            self._symmetrizer = tuple([1]*size)
        else:
            self._symmetrizer = tuple(symmetrizer)
        cq=ClusterQuiver(self._b_matrix)
        cq.mutate(flatten(self._vertices))
        if cq.b_matrix() != matrix(size,lambda i,j:self._b_matrix[permutation.inverse()(i+1)-1][permutation.inverse()(j+1)-1]):
            raise ValueError('The first quiver and the last quiver are distinct.')

    def __eq__(self, other):
        b1, b2 = self._b_matrix, other._b_matrix
        i1, i2 = tuple(flatten(self._vertices)), tuple(flatten(self._vertices))
        nu1, nu2 = self._permutation, other._permutation
        d1, d2 = self._symmetrizer, other._symmetrizer
        return b1 == b2 and i1 == i2 and nu1 == nu2 and d1 == d2

    def b_matrix(self,u):
        r"""
        Return the exchange matrix at the time ``u``.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.b_matrix(1)
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]
            sage: gamma.b_matrix(2)
            [ 0 -1  0]
            [ 1  0  1]
            [ 0 -1  0]

        """
        if u==0:
            return self._b_matrix
        cq = ClusterQuiver(self._b_matrix)
        t=len(self._vertices)
        p=self.permutation()
        if u>0:
            m=u//t
            k=u%t
            for tt in range(k):
                cq.mutate(self.indices(tt))
            pm=Permutation(range(1,1+self.size())) if m==0 else reduce(lambda a,b:a*b,[p]*m)
            return matrix(self.size(), lambda i,j: cq.b_matrix()[pm.inverse()(i+1)-1][pm.inverse()(j+1)-1])
        elif u<0:
            k=1
            while True:
                if u+k*self._order_of_permutation()*t>=0:
                    return self.b_matrix(u+k*self._order_of_permutation()*t)
                k+=1

    def initial_b_matrix(self):
        r"""
        Return the initial exchange matrix $B$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.initial_b_matrix()
            [ 0 -1  0]
            [ 1  0  1]
            [ 0 -1  0]

        """
        return self.b_matrix(0)

    def symmetrizer(self):
        r"""
        Return the symmetriser $d$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.symmetrizer()
            (1, 1, 1)

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-2,0])
            sage: d = [1,1,2]
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu,d)
            sage: gamma.symmetrizer()
            (1, 1, 2)
        
        """
        return self._symmetrizer

    
    def vertices(self):
        r"""
        Return the sequence $\mathbf{i} = (\mathbf{i}(0), \dots, \mathbf{i}(t-1))$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.vertices()
            [[0, 2], [1]]

        """
        return self._vertices[:]

    def size(self):
        r"""
        Return the size of exchange matrices in the mutation loop.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.size()
            3
        
        """
        return self._b_matrix.ncols()

    def length(self):
        r"""
        Return the integer $t$ where $\mathbf{i} = (\mathbf{i}(0),\dots,\mathbf{i}(t-1))$.

        EXAMPLES::
            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.length()
            2

        """
        return len(self._vertices)

    def whole_length(self):
        r"""
        Return the length of $\mathbf{i}$, which is the sum of the lenghs $\mathbf{i}(0),\dots,\mathbf{i}(t-1)$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.whole_length()
            3

        """
        return len(flatten(self._vertices))
    
    def inverse(self):
        r"""
        Return the mutation loop $(B,d,(\mathbf{t-1},\dots,\mathbf{i}(0)),\nu^{-1})$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,2,-2,-2,0,2,2,-2,0])
            sage: MutationLoop(b,[[1],[0]],Permutation([2,3,1])).inverse() == MutationLoop(b,[[1],[2]],Permutation([3,1,2]))
            True

        """
        i_inv = [[self.permutation()(v+1)-1 for v in vs] for vs in self._vertices[::-1]]
        return MutationLoop(self.b_matrix(0),i_inv,self.permutation().inverse())

    def permutation(self):
        r"""
        Return the permutation $\nu$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.permutation()
            [1, 2, 3]
            sage: nu2 = Permutation([3,2,1])
            sage: gamma2 = MutationLoop(b,i,nu2)
            sage: gamma2.permutation()
            [3, 2, 1]

        """
        return self._permutation
    
    def indices(self,u):
        r"""
        Return the sequence $\mathbf{i}(u)$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.indices(0)
            [0, 2]

        """
        return self._vertices[u]
    
    @cached_method
    def _order_of_permutation(self):
        r"""
        Return the order of the permutation $\nu$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu1 = Permutation([1,2,3])
            sage: nu2 = Permutation([3,2,1])
            sage: gamma1 = MutationLoop(b,i,nu1)
            sage: gamma1._order_of_permutation()
            1
            sage: gamma2 = MutationLoop(b,i,nu2)
            sage: gamma2._order_of_permutation()
            2

        """
        p=self.permutation()
        k=1
        while True:
            if p==Permutation(range(1,1+self.size())):
                return k
            p = p * self.permutation()
            k+=1

    def _time_of_mutation(self,a):
        r"""
        Return the integer $u$ such that $i_a \in \mathbf{i}(u)$
        where $i_a$ is the $a$-th element of $\mathbf{i}$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: [gamma._time_of_mutation(a) for a in range(3)]
            [0, 0, 1]

        """
        length_of_vertices = [len(vs) for vs in self._vertices]
        t=0
        whole_l=0
        for l in length_of_vertices:
            whole_l+=l
            if a<whole_l:
                return t
            t+=1

    @staticmethod
    @cached_method
    def _permutation_power(permutation,exponent):
        r"""
        Return the ``exponent``-th power of ``permutation``.

        EXAMPLES::
            sage: MutationLoop._permutation_power(Permutation([2,3,1]),3)
            [1, 2, 3]

        """
        if exponent>=0:
                return reduce(lambda a,b:a*b,[permutation]*exponent,Permutation(range(1,permutation.size()+1)))
        else:
            return MutationLoop._permutation_power(permutation.inverse(),-exponent)
    
    def latency(self,i,u):
        r"""
        Return the latency of the pair (``i``, ``u``).

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.latency(0,0)
            0
            sage: gamma.latency(1,0)
            1

        """
        for v in range(self._order_of_permutation()*self.length()+3):
            k=(v+u)%self.length()
            m=(v+u)//self.length()
            if i in [MutationLoop._permutation_power(self.permutation(),m)(w+1)-1 for w in self._vertices[k]]:
                return v
        return PlusInfinity()
    
    def _colatency(self,i,u):
        r"""
        Return the colatency of the pair (``i``, ``u``).

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma._colatency(0,0)
            0
            sage: gamma._colatency(1,0)
            1
            
        """
        for v in range(0,-self._order_of_permutation()*self.length()-3,-1):
            k=(v+u)%self.length()
            m=(v+u)//self.length()
            if i in [MutationLoop._permutation_power(self.permutation(),m)(w+1)-1 for w in self._vertices[k]]:
                return -v
        return PlusInfinity()
    
    def _strict_latency(self,i,u):
        r"""
        Return the strict latency of the pair (``i``, ``u``).

        EXAMPLES::
        
            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma._strict_latency(0,0)
            2
            sage: gamma._strict_latency(1,0)
            1

        """
        if self.latency(i,u)!=0:
            return self.latency(i,u)
        elif self.latency(i,u)==0:
            return self.latency(i,u+1)+1
        
    def _strict_colatency(self,i,u):
        r"""
        Return the strict colatency of the pair (``i``, ``u``).

        EXAMPLES::
        
            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma._strict_colatency(0,0)
            2
            sage: gamma._strict_colatency(1,0)
            1

        """
        if self._colatency(i,u)!=0:
            return self._colatency(i,u)
        elif self._colatency(i,u)==0:
            return self._colatency(i,u-1)+1
    
    def is_complete(self):
        r"""
        Return ``True`` if all latencies are finite.

        EXAMPLES::
            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.is_complete()
            True
            sage: i2 = [[1],[1]]
            sage: gamma2 = MutationLoop(b,i2,nu)
            sage: gamma2.is_complete()
            False

        """
        return all([self.latency(i,0) < PlusInfinity() for i in range(self.size())])
            
    def _mutation_points_equivalence(self,i,u,j,v):
        r"""
        Return ``True``` if there is a integer $g$ such that $j=\nu^g(i)$ and $v = u+gt$,
        where $nu$ = ``self.permutation()`` and $t$ = ``self.length()``.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([3,2,1])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma._mutation_points_equivalence(0,0,2,2)
            True
            sage: gamma._mutation_points_equivalence(0,0,0,2)
            False

        """
        nu=self.permutation()
        if (v-u)%self.length()!=0:
            return False
        time = (v-u)//self.length()
        return  j==MutationLoop._permutation_power(nu,time)(i+1)-1
        
    def _pi(self,i,u):
        r"""
        Return the integer $\pi(i,u) \in \{ 1, \dors, r \}$.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma._pi(0,0)
            0
            sage: gamma._pi(2,0)
            1
            sage: gamma._pi(1,1)
            2
            sage: gamma._pi(0,2)
            0

        """
        if self.latency(i,u) != 0:
            raise ValueError('({},{}) is not a mutation point.'.format(i,u))
        elif self.latency(i,u) == 0:
            for a,ia in enumerate(flatten(self._vertices)):
                ua=self._time_of_mutation(a)
                if self._mutation_points_equivalence(i,u,ia,ua):
                    return a
    
    def next_mutation_point(self,i,u):
        r"""
        Return the next mutation point of (``i``, ``u``).

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.next_mutation_point(0,0)
            (0, 2)
            sage: gamma.next_mutation_point(1,1)
            (1, 3)

        """
        if self.latency(i,u)==PlusInfinity():
            raise ValueError('The latency of ({},{}) is inifinity.'.format(i,u))
        else:
            return (i,u+self._strict_latency(i,u))
        
    def _previous_mutation_point(self,i,u):
        r"""
        Return the previous mutation point of (``i``, ``u``).

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma._previous_mutation_point(0,0)
            (0, -2)
            sage: gamma._previous_mutation_point(1,1)
            (1, -1)

        """
        if self._colatency(i,u)==PlusInfinity():
            raise ValueError('The latency of ({},{}) is inifinity.'.format(i,u))
        else:
            return (i,u-self._strict_colatency(i,u))

    def _symmetrizer_at_mutation_points(self,i,u):
        r"""
        Return the symmetrizing positive integer at the mutation point (``i``, ``u``).  

        """
        return self.symmetrizer()[i]
    
    def t_datum(self,variable_name='z'):
        r"""
        Return the T-datum of the mutation loop.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-1,0])
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu)
            sage: gamma.t_datum().pair()
            (
            [1 + z^2       0       0]  [1 + z^2       0      -z]
            [      0 1 + z^2       0]  [      0 1 + z^2      -z]
            [      0       0 1 + z^2], [     -z      -z 1 + z^2]
            )

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-2,0])
            sage: d = [1,1,2]
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu,d)
            sage: gamma.t_datum().pair()
            (
            [1 + z^2       0       0]  [1 + z^2       0      -z]
            [      0 1 + z^2       0]  [      0 1 + z^2    -2*z]
            [      0       0 1 + z^2], [     -z      -z 1 + z^2]
            )

        """
        NT0,NTp,NTm = self._T_system_triple()
        Ap=NT0-NTp
        Am=NT0-NTm
        D=diagonal_matrix(self.symmetrizer()[i] for i in  flatten(self._vertices))
        return TDatum(Ap,Am,D)
    
    def _T_system_triple(self,variable_name='z'):
        r"""
        Return the triple of matrices that describe the T-system of the mutation loop.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-2,0])
            sage: d = [1,1,2]
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu,d)
            sage: gamma._T_system_triple()
            (
            [1 + z^2       0       0]  [0 0 0]  [  0   0   z]
            [      0 1 + z^2       0]  [0 0 0]  [  0   0 2*z]
            [      0       0 1 + z^2], [0 0 0], [  z   z   0]
            )

        """
        var = LaurentPolynomialRing(QQ,variable_name).gen()
        flt_vs=flatten(self._vertices)
        def entry(a,b,sign):
            k = flt_vs[a]
            u = self._time_of_mutation(a)
            result = 0
            if sign==0:
                if a==b:
                    result += 1
                if self._pi(*self.next_mutation_point(k,u))==b:
                    result += var**(self._strict_latency(k,u))
                return result
            elif sign==1 or -1:
                return sum([max(0,-sign*self.b_matrix(u)[j][k]) * var**(self.latency(j,u)) for j in range(self.size()) if self.latency(j,u)<PlusInfinity() and self._pi(*self.next_mutation_point(j,u))==b])
            else:
                raise ValueError('{} is not a sign.'.format(sign))
        N0 = matrix(LaurentPolynomialRing(QQ,variable_name),self.whole_length(),lambda a,b:entry(b,a,0))        
        Np = matrix(LaurentPolynomialRing(QQ,variable_name),self.whole_length(),lambda a,b:entry(b,a,1))
        Nm = matrix(LaurentPolynomialRing(QQ,variable_name),self.whole_length(),lambda a,b:entry(b,a,-1))
        return (N0,Np,Nm)
        
    def _Y_system_triple(self,variable_name='z'):
        r"""
        Return the triple of matrices that describe the Y-system of the mutation loop.

        EXAMPLES::

            sage: b = matrix(3,3,[0,-1,0,1,0,1,0,-2,0])
            sage: d = [1,1,2]
            sage: i = [[0,2],[1]]
            sage: nu = Permutation([1,2,3])
            sage: gamma = MutationLoop(b,i,nu,d)
            sage: gamma._Y_system_triple()
            (
            [1 + z^2       0       0]  [0 0 0]  [  0   0   z]
            [      0 1 + z^2       0]  [0 0 0]  [  0   0   z]
            [      0       0 1 + z^2], [0 0 0], [  z 2*z   0]
            )
            
        """
        var = LaurentPolynomialRing(QQ,variable_name).gen()
        flt_vs=flatten(self._vertices)
        def entry(a,b,sign):
            k = flt_vs[a]
            u = self._time_of_mutation(a)
            result = 0
            if sign==0:
                if a==b:
                    result += 1
                if self._pi(*self._previous_mutation_point(k,u))==b:
                    result += var**(self._strict_colatency(k,u))
                return result
            elif sign==1 or -1:
                from itertools import product
                return sum([max(0,sign*self.b_matrix(v)[j][k]) * var**(self.latency(k,v)) for j,v in product(range(self.size()),range(u+1-self._strict_colatency(k,u),u)) if self.latency(j,v)==0 and self._pi(j,v)==b])
            else:
                raise ValueError('{} is not a sign.'.format(sign))
        N0 = matrix(LaurentPolynomialRing(QQ,variable_name),self.whole_length(),lambda a,b:entry(a,b,0))
        Np = matrix(LaurentPolynomialRing(QQ,variable_name),self.whole_length(),lambda a,b:entry(a,b,1))
        Nm = matrix(LaurentPolynomialRing(QQ,variable_name),self.whole_length(),lambda a,b:entry(a,b,-1))
        return (N0,Np,Nm)
