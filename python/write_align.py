#!/usr/bin/env python
# small script to convert the numerous alignment files into a single one
# the old alignment data can be found in data/old_alignments.tar.gz

from glob import glob


def remove_letters(word):
    return ''.join([c for c in word if c.isdigit()])


file_names = glob('telescope*.dat')
file_names = sorted(file_names, key=lambda x: int(remove_letters(x)))

lines = ['#TEL CH ROC  RZ           RY           X            Y            Z           (dX)         (dY)\n']  # header

for name in file_names:
    with open(name) as f:
        tel_id = remove_letters(name).rjust(3)
        lines.append('# TELESCOPE {} {}\n'.format(tel_id.strip(), '-' * 88))
        for line in f.readlines():
            if len(line) > 10 and not line.startswith('#'):
                data = line.split()
                if data[1] == '-1' and float(data[2]) == 0:
                    continue  # skip zero offset lines
                align = ['{:+1.4e}'.format(float(word)) for word in data[2:]]
                info = '{}   {}'.format(tel_id, '   '.join(data[:2]))
                lines.append('{}  {}  {}{}\n'.format(info, align[0], ' ' * 13 if len(align) in [4, 6] else '', '  '.join(align[1:])))

with open('alignments.txt', 'w') as f:
    f.writelines(lines)

