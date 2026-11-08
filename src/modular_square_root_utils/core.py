"""Modular square roots modulo an odd prime.

Implements the Tonelli-Shanks algorithm.  The public entry point is
:func:`modular_sqrt`, which returns a square root of ``a`` modulo the prime
``p`` or ``None`` when ``a`` is a quadratic non-residue.

Why Tonelli-Shanks?  It is the standard general-purpose algorithm for odd
primes.  The special cases ``p == 2`` and ``p % 4 == 3`` are handled directly
because they are common and much faster than the general loop.
"""

from __future__ import annotations

__all__ = ["modular_sqrt", "legendre_symbol"]


def legendre_symbol(a: int, p: int) -> int:
    """Return the Legendre symbol ``(a | p)`` for an odd prime ``p``.

    The result is ``1`` if ``a`` is a quadratic residue modulo ``p``,
    ``-1`` if it is a non-residue, and ``0`` if ``a`` is divisible by ``p``.

    This uses Euler's criterion: ``a ** ((p - 1) // 2) mod p``.  For large
    primes that exponentiation is fast enough, and it avoids needing a
    separate factoring step.
    """
    if p <= 2:
        raise ValueError("p must be an odd prime")
    a %= p
    if a == 0:
        return 0
    value = pow(a, (p - 1) // 2, p)
    # pow returns 1 or p-1 (which is -1 mod p).
    return 1 if value == 1 else -1


def modular_sqrt(a: int, p: int) -> int | None:
    """Return one square root of ``a`` modulo the odd prime ``p``.

    If ``a`` is a quadratic non-residue modulo ``p`` the function returns
    ``None``.  If ``a`` is ``0`` modulo ``p`` the result is ``0``.  When two
    roots exist (``r`` and ``p - r``) the smaller one is returned so that the
    result is deterministic and independent of the internal search order.

    The prime ``p`` must be an odd integer greater than ``2``.  Passing a
    composite modulus is not supported and is not detected; the algorithm
    assumes primality.
    """
    if p <= 2:
        raise ValueError("p must be an odd prime")
    a %= p
    if a == 0:
        return 0
    if legendre_symbol(a, p) != 1:
        return None

    # Fast path: p % 4 == 3.  Then a root is a ** ((p + 1) // 4) mod p.
    if p % 4 == 3:
        root = pow(a, (p + 1) // 4, p)
        return min(root, p - root)

    # General Tonelli-Shanks.
    # Write p - 1 = q * 2**s with q odd.
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1

    # Find a quadratic non-residue z.  The smallest one is almost always
    # tiny, so a linear scan is fine in practice.
    z = 2
    while legendre_symbol(z, p) != -1:
        z += 1

    m = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)

    while t != 1:
        # Find the least i, 0 < i < m, with t ** (2**i) == 1.
        i = 0
        t2i = t
        while t2i != 1:
            t2i = (t2i * t2i) % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p

    return min(r, p - r)
