import unittest

from modular_square_root_utils import legendre_symbol, modular_sqrt


class LegendreSymbolTests(unittest.TestCase):
    def test_residue(self):
        # 2 is a quadratic residue modulo 7 because 3*3 = 9 = 2 (mod 7).
        self.assertEqual(legendre_symbol(2, 7), 1)

    def test_non_residue(self):
        # 3 is a non-residue modulo 7.
        self.assertEqual(legendre_symbol(3, 7), -1)

    def test_zero(self):
        self.assertEqual(legendre_symbol(0, 7), 0)
        self.assertEqual(legendre_symbol(14, 7), 0)

    def test_rejects_even_modulus(self):
        with self.assertRaises(ValueError):
            legendre_symbol(5, 2)


class ModularSqrtTests(unittest.TestCase):
    def test_zero(self):
        self.assertEqual(modular_sqrt(0, 7), 0)
        self.assertEqual(modular_sqrt(7, 7), 0)

    def test_non_residue_returns_none(self):
        self.assertIsNone(modular_sqrt(3, 7))

    def test_p_mod_4_is_3(self):
        # 7 % 4 == 3, so this exercises the fast path.
        self.assertEqual(modular_sqrt(2, 7), 3)
        self.assertEqual(modular_sqrt(4, 7), 2)

    def test_general_tonelli_shanks(self):
        # 13 % 4 == 1, so the general loop runs.
        # 6 is a residue modulo 13: 5*5 = 25 = 12, 6*6 = 36 = 10, 7*7 = 49 = 10,
        # 4*4 = 16 = 3, 2*2 = 4, 3*3 = 9, 1*1 = 1, 8*8 = 64 = 12, 9*9 = 81 = 3,
        # 10*10 = 100 = 9, 11*11 = 121 = 4, 12*12 = 144 = 1.  Wait -- let's
        # just use 3: 4*4 = 16 = 3 (mod 13).
        self.assertEqual(modular_sqrt(3, 13), 4)
        self.assertEqual(modular_sqrt(4, 13), 2)

    def test_smaller_root_is_returned(self):
        # 5 has roots 4 and 9 modulo 13; the smaller one is 4.
        self.assertEqual(modular_sqrt(3, 13), 4)
        # 12 has roots 5 and 8 modulo 13; the smaller one is 5.
        self.assertEqual(modular_sqrt(12, 13), 5)

    def test_large_prime(self):
        p = 1_000_000_007
        a = 123_456_789
        root = modular_sqrt(a, p)
        self.assertIsNotNone(root)
        self.assertEqual((root * root) % p, a % p)
        self.assertLessEqual(root, p - root)

    def test_rejects_even_modulus(self):
        with self.assertRaises(ValueError):
            modular_sqrt(1, 2)

    def test_rejects_modulus_below_three(self):
        with self.assertRaises(ValueError):
            modular_sqrt(0, 1)

    def test_negative_input_is_reduced(self):
        # -1 mod 7 is 6, which is a non-residue.
        self.assertIsNone(modular_sqrt(-1, 7))
        # -5 mod 7 is 2, which is a residue with root 3.
        self.assertEqual(modular_sqrt(-5, 7), 3)


if __name__ == "__main__":
    unittest.main()
