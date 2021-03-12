#!/usr/bin/env python
# small script to add new telescope
# makes an entry to the data/config.txt

from argparse import ArgumentParser
from ROOT import TFile
from numpy import frombuffer
from os.path import dirname, realpath, join


def get_n_rocs(filename):
    f = TFile(filename)
    t = f.Get('tree')
    t.SetEstimate(10 * t.GetEntries())
    n = t.Draw('plane', '', 'goff')
    return int(max(frombuffer(t.GetVal(0), count=n)) + 1)


def get_last_tel(filename):
    with open(filename) as f:
        return max([int(line.split()[0]) for line in f.readlines() if not line.startswith('#') and len(line) > 5])


def get_z_pos(year):
    if year < 2016:
        return 0
    elif year == 2016:
        return 1
    elif year > 2016:
        return 2


def save_config(filename, n_rocs, tel_id):
    last_tel = get_last_tel(filename)
    new_tel = input('Enter the telescope number, press enter for default ({}): '.format(last_tel + 1)) if tel_id is None else tel_id
    tel = last_tel + 1 if not new_tel else int(new_tel)
    mask = input('Enter the outer mask file number, press enter for default: ')
    mask = 0 if not mask else int(mask)
    cal = input('Enter the calibration number (press enter for default ({})): '.format(tel))
    cal = tel if not cal else int(cal)
    year = int(input('Enter the year (YYYY): '))
    z_pos = get_z_pos(year)
    cmt = input('Enter a comment (Month, DUT, ...): ')
    type_ = ['PAD', 'PIX', 'BCM'][int(input('Enter the type of the DUT (0 for Pad, 1 for Pixel and 2 for BCM\': '))]
    with open(filename, 'a') as f:
        f.write('\n{: 3d}{: 4d}{: 8d}{: 5d}{: 6d}{: 6d} {}  # {}'.format(tel, n_rocs, mask, cal, z_pos, year, type_, cmt))


if __name__ == '__main__':

    d = dirname(dirname(realpath(__file__)))

    p = ArgumentParser()
    p.add_argument('run', help='path to the ROOT file that is used for alignment')
    p.add_argument('id', nargs='?', help='telescope id', default=None)
    args = p.parse_args()

    config = join(d, 'config', 'telescopes.txt')
    save_config(config, get_n_rocs(args.run), args.id)
