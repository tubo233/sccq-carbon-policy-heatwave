"""Join the numbered parts (.001, .002, ...) of this record into the three zip files and check them.

Usage: put all downloaded files in one folder and run  python join_parts.py
"""
import glob
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ['Research_package_main.zip', 'Research_package_tmax.zip', 'Research_package_tmin.zip']


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    expected = {}
    with open(os.path.join(HERE, 'MD5SUMS.txt'), encoding='utf-8') as f:
        for line in f:
            digest, name = line.split()
            expected[name] = digest
    ok = True
    for name in NAMES:
        parts = sorted(glob.glob(os.path.join(HERE, name + '.[0-9][0-9][0-9]')))
        if not parts:
            print(f'{name}: no parts found'); ok = False; continue
        for p in parts:
            if md5(p) != expected.get(os.path.basename(p)):
                print(f'{os.path.basename(p)}: checksum mismatch, download it again'); ok = False
        out = os.path.join(HERE, name)
        with open(out, 'wb') as w:
            for p in parts:
                with open(p, 'rb') as r:
                    w.write(r.read())
        good = md5(out) == expected[name]
        ok &= good
        print(f'{name}: {len(parts)} parts joined, checksum {"OK" if good else "MISMATCH"}')
    print('All files restored; extract the three zip files into one folder.' if ok else 'Some files are missing or damaged.')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
