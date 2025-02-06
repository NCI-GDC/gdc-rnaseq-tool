#!/usr/bin/env python3

import unittest

from click.testing import CliRunner
from gdc_rnaseq_tool import __main__ as MOD


class ThisTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.runner = CliRunner()

    def test_pass(self) -> None:
        result = self.runner.invoke(MOD.main)
        self.assertEqual(result.exit_code, 0)


# __END__
