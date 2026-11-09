# modular_square_root_utils

Computes a square root of an integer modulo an odd prime, returning `None` when the integer is a quadratic non-residue.

```python
from modular_square_root_utils import modular_sqrt, legendre_symbol

# 2 is a quadratic residue modulo 7, with roots 3 and 4.
assert modular_sqrt(2, 7) == 3

# 3 is a non-residue modulo 7.
assert modular_sqrt(3, 7) is None

# The Legendre symbol tells you which case you are in.
assert legendre_symbol(2, 7) == 1
assert legendre_symbol(3, 7) == -1
```

## Why this exists

Solving `x*x = a (mod p)` for an odd prime `p` comes up in elliptic-curve point decompression, in quadratic-residue-based hashing, and in number-theory exercises. Python's standard library has no modular square root, and the usual workarounds either pull in a big dependency or restrict you to the easy case `p % 4 == 3`. This package implements Tonelli-Shanks directly, so it works for every odd prime, and keeps the fast path for `p % 4 == 3` because that case is common in practice.

The trade-off is deliberate: the code assumes `p` is prime and does not verify it. Primality testing is a separate, more expensive problem, and callers who need it can test beforehand. Passing a composite modulus produces results that are meaningless rather than an error.

## One awkward edge

When a residue has two roots, `r` and `p - r`, the function returns the smaller one. That makes the result deterministic, but it means you cannot rely on the function returning a particular root from the pair; if you need the other one, compute `p - modular_sqrt(a, p)`.

## API

- `modular_sqrt(a: int, p: int) -> int | None` — returns a root of `a` modulo the odd prime `p`, or `None` if `a` is a non-residue. Returns `0` when `a` is `0` modulo `p`. Raises `ValueError` if `p <= 2`.
- `legendre_symbol(a: int, p: int) -> int` — returns `1`, `-1`, or `0` following Euler's criterion. Raises `ValueError` if `p <= 2`.
