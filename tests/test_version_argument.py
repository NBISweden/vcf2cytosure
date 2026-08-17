from unittest.mock import patch

import pytest

from vcf2cytosure.vcf2cytosure import main


def test_version_argument():
    with patch("sys.argv", ["vcf2cytosure", "--version"]):
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 0
