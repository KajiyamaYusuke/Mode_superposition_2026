import importlib.util
from decimal import Decimal
from pathlib import Path
import unittest


path = Path(__file__).resolve().parents[1] / 'body_ratio_jitter_261103' / 'inventory_modes.py'
spec = importlib.util.spec_from_file_location('body_ratio_inventory', path)
inventory = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inventory)


class TokenTests(unittest.TestCase):
    def test_decimal_names(self):
        for token, expected in [('b1p2c3', '1.2'), ('b11p4c3', '11.4'),
                                ('b10p8c3', '10.8'), ('b12c3', '12'),
                                ('b1p20c3', '1.20')]:
            self.assertEqual(inventory.parse_token(token), (Decimal(expected), 3))

    def test_invalid_tokens(self):
        for token in ('b1pc3', 'b1p2c', 'prefixb1p2c3', 'b1p2c3extra'):
            with self.assertRaises(ValueError):
                inventory.parse_token(token)

    def test_equivalent_names_and_cover_remain_distinct(self):
        self.assertEqual(inventory.parse_token('b3c3'), inventory.parse_token('b3p0c3'))
        self.assertNotEqual(inventory.parse_token('b3c2'), inventory.parse_token('b3c3'))


if __name__ == '__main__':
    unittest.main()
