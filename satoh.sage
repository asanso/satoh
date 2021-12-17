def _miller_(P, Q, n):
        if Q.is_zero():
            raise ValueError("Q must be nonzero.")
        if n.is_zero():
            raise ValueError("n must be nonzero.")
        n_is_negative = False
        if n < 0:
            n = n.abs()
            n_is_negative = True

        one = E.base_field().one()
        t = one
        V = P
        nbin = n.bits()
        i = n.nbits() - 2
        while i > -1:
            S = 2*V
            ell = V._line_(V, Q)
            vee = S._line_(-S, Q)
            t = (t**2)*(ell/vee)
            V = S
            if nbin[i] == 1:
                S = V+P
                ell = V._line_(P, Q)
                vee = S._line_(-S, Q)
                t = t*(ell/vee)
                V = S
            i = i-1
        if n_is_negative:
            vee = V._line_(-V, Q)
            t = 1/(t*vee)
        return t

def recover(v, P):
    u = v^70
    x1 = P[0] + u
    x2 = P[0] -u
    L1 =  E.lift_x(x1)
    L2 =  E.lift_x(x2)
    if _miller_(P,L1,140) == v:
        return L1
    if _miller_(P,-L1,140) == v:
        return -L1
    if _miller_(P,L2,140) == v:
        return L2
    if _miller_(P,-L2,140) == v:
        return -L2

p =139
_.<I> = GF(p)[]
K.<i> = GF(p^2, modulus=I^2+4)
E = EllipticCurve(K, [-13, -7])

P = E.lift_x(67) 
P = -P

a =  25*i + 109 
Pa = recover(a,P)

b = 112*i +22
Pb = recover(b,P)

recover(a,Pb) == recover(b,Pa)
