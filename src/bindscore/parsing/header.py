import sys
import pathlib

# test_enthalpy.py'nin bulunduğu yerden root klasörüne (BindScore) çıkıyoruz
ROOT = pathlib.Path(__file__).parent.parent.resolve()

# Gitmesini istediğin tüm alt klasörleri buraya tek seferde tanımlıyoruz
SRC_PATHS = [
    ROOT / 'src' / 'bindscore' / 'pdb_file_treatment',
    ROOT / 'src' / 'bindscore' / 'parsing',
    ROOT / 'src' / 'bindscore' / 'scoring'
]

# Hepsini Python'un arama listesinin en başına çakıyoruz
for path in SRC_PATHS:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))