import pytest

from unittest.mock import patch

import vcf2cytosure

def test_version_argument():
    with patch('sys.argv', ['vcf2cytosure.py','--version']):
        with pytest.raises(SystemExit) as excinfo:
            vcf2cytosure.main()
        assert excinfo.value.code == 0
