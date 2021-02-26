import sys

with open(sys.argv[1]+'_ref.fasta', 'r') as f:
    lines = []
    for i in f:
        lines.append(i.rstrip('\n'))
    combined = ''.join(lines[1:])

with open(sys.argv[1]+'_ref_clean.fasta', 'w') as out:
    out.write('>'+sys.argv[1]+'\n')
    out.write(combined+'\n')
    out.close()
