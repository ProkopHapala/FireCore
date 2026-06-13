#!/usr/bin/env python3
import re, sys

def convert(text):
    # Replace unicode subscripts/superscripts
    sub_map = {'₀':'_0','₁':'_1','₂':'_2','₃':'_3','₄':'_4','₅':'_5','₆':'_6','₇':'_7','₈':'_8','₉':'_9',
               '⁰':'^0','¹':'^1','²':'^2','³':'^3','⁴':'^4','⁵':'^5','⁶':'^6','⁷':'^7','⁸':'^8','⁹':'^9'}
    for k,v in sub_map.items(): text = text.replace(k, v)
    # Inline math: backtick-enclosed greek/math patterns
    pattern = re.compile(r'`([^`]*(?:[Ψψ]|\\[a-zA-Z]+|\^|_|\\nabla|\\partial|\\times|\\dagger|=|\(|\)|\+|\\|\*|/)[^`]*)`')
    text = pattern.sub(lambda m: '$%s$' % m.group(1), text)
    return text

if __name__ == '__main__':
    if len(sys.argv)<2:
        print('Usage: convert_md_math.py <file.md>')
        sys.exit(1)
    fname = sys.argv[1]
    txt = open(fname, 'r', encoding='utf-8').read()
    new = convert(txt)
    open(fname, 'w', encoding='utf-8').write(new)
