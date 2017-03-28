import re
import json
import argparse

from bs4 import BeautifulSoup


class Strings:
    def __init__(self, value=""):
        if isinstance(value, Strings):
            self.value = value.value
        else:
            self.value = value

    def replace(self, pattern, repl):
        self.value = re.sub(pattern, repl, self.value)
        return self

class PhylotreeParser:
    class DuplicateDetectionError(Exception): pass

    class Pattern:
        ATGC = re.compile("[atgcATGC]")

        REGULAR = re.compile('''
        \(?   # unstable
            {}  # 祖先
            \d+ # 位置
            {}  # 子孫
            !*  # reversion
        \)?'''.format(ATGC.pattern, ATGC.pattern), re.X)

        IRREGULAR = re.compile('''\A
        \(?                     # unstable
          (?:
            # ex: C459d
            {}?
            \d+
            d                   # deletion

          |
            # ex: 960.XC
            \d+
            \.[\dX]?{}+         # insertion
            d?                  # deletion

          |
            # ex: 59-60d
            \d+\-\d+d           # range

          |
            # ネームスペースだけ予約？
            reserved
          )
        !*                      # reversion
        \)?
        \Z'''.format(ATGC.pattern, ATGC.pattern), re.X)

        EXCEPTIONS = re.compile("|".join([re.escape(i) for i in re.split("\s+", '''T65d G71d A249d C299d C309d A337d
          T455d C456d C459d C498d C960d
          A1409d A1656d A2074d A2395d A4317d
          A5752d A5894d C5899d C7471d T15944d
          A16166d C16187d C16193d C16257d T16325d

          59-60d 105-110d 106-111d 290-291d 291-294d 8281-8289d

          573.XC 960.XC 965.XC 5899.XC 8278.XC  5899.XCd!

          44.1C 5899.1C 459.1C 15944.1T
          374.1A 455.1T 2232.1A 16259.1A
          595.1C 1719.1G 3172.1C 191.1A
          2156.1A 55.1T 65.1T 60.1T 42.1G
          12310.1A 291.1A 5752.1A 310.1T
          597.1T 3229.1A 2405.1C 8276.1C
          2484.1C 745.1T 498.1C 16169.1C
          5899.1C 8279.1T 5740.1A 960.1C
          3307.1A 93.1T 456.1T 356.1C
          3158.1T

          44.1C 55.1T 60.1T 65.1T 93.1T
          191.1A 291.1A 310.1T 356.1C 374.1A
          455.1T 456.1T 459.1C 498.1C 595.1C
          597.1T 745.1T 960.1C 1719.1G 2156.1A
          2232.1A 2405.1C 2484.1C 3158.1T 3172.1C
          3229.1A 3307.1A 5740.1A 5752.1A 5899.1C
          5899.1C 8276.1C 8279.1T 12310.1A 15944.1T
          16169.1C 16259.1A
          C5899.1d!  459.1Cd!
          60.1TT 292.1AT 368.1AGAA
          8289.1CCCCCTCTA 8289.1CCCCCTCTACCCCCTCTA

          455.2T 2232.2A

          (573.XC) (745.1T) (960.1C)
          (C965d) (C16193d)
          reserved''')]), re.X)

        BRANCH_CONDITION = re.compile('''(?:
            {}
            |
            {}
        )'''.format(REGULAR.pattern, EXCEPTIONS.pattern), re.X)

        BRANCH_CONDITIONS = re.compile('''\A
            {}
            (?:\s{})*
        \Z'''.format(BRANCH_CONDITION.pattern, BRANCH_CONDITION.pattern), re.X)

        TABLE_TITLE = re.compile("mt\-MRCA")

        @staticmethod
        def is_irregular(text):
            return PhylotreeParser.Pattern.IRREGULAR.search(text)

        @staticmethod
        def is_table_title(text):
            return PhylotreeParser.Pattern.TABLE_TITLE.search(text)

        @staticmethod
        def is_branch_conditions(text):
            return PhylotreeParser.Pattern.BRANCH_CONDITIONS.search(text)

    BLANK_BRANCH_NAME = '__BRANCH__'
    SELF_BRANCH_NAME = '__SELF__'

    def __init__(self, filename):
        self.file = filename
        self.in_useful_range = False
        self.queue = {} # []
        self.raw_tree = {}

    def parse(self):
        for table in self.parse_html().find_all('table'):
            self.process_table(table)
        return self

    def parse_html(self):
        return BeautifulSoup(self.read_file(), 'html.parser')

    def read_file(self):
        with open(self.file, 'r', encoding='windows-1252') as f:
            return f.read()

    def process_table(self, table):
        for tr in table.find_all('tr'):
            if not self.is_in_useful_range(tr.get_text()):
                continue 
            self.process_tr(tr)

    def process_tr(self, tr):
        cols = []
        conditions = []
        depth = 0
        haplogroup_candidates = []

        for (i, td) in enumerate(tr.find_all('td')):
            raw = td.get_text()
            text = self.trim_white_space(raw)
            cols.append(text)
            if len(text) == 0:
                continue

            if PhylotreeParser.Pattern.is_branch_conditions(text):
                conditions = text.split(' ')
                depth = i - 1
            else:
                haplogroup_candidates.append(text)

        if len(conditions) == 0:
            return

        example_accessions = self.extract_example_accessions(cols)
        haplogroup = self.detect_haplogroup(haplogroup_candidates, example_accessions)
        self.grow_tree(haplogroup, conditions, example_accessions, depth)

    def extract_example_accessions(self, cols):
        return [cols[i] for i in [-2, -1] if cols[i]]

    def detect_haplogroup(self, candidates, example_accessions):
        haplogroup = ''
        for candidate in candidates:
            if candidate in example_accessions:
                continue
            if haplogroup:
                raise DuplicateDetectionError
            haplogroup = candidate
        if not haplogroup:
            haplogroup = PhylotreeParser.BLANK_BRANCH_NAME

        return haplogroup

    def grow_tree(self, haplogroup, conditions, example_accessions, depth):
        self.queue[depth] = haplogroup
        node = self.get_deep_hash(self.raw_tree, self.queue, list(range(depth + 1)))
        node[PhylotreeParser.SELF_BRANCH_NAME] = {
                    'conditions': conditions,
                    'example_accessions': example_accessions
                }

    def get_deep_hash(self, dic, arr, itr):
        idx, itr = itr[0], itr[1:]
        key = arr[idx]
        if not key in dic:
            dic[key] = {}
        next_dic = dic[key]
        if len(itr) == 0:
            return next_dic
        return self.get_deep_hash(next_dic, arr, itr) 

    def prettify(self):
        return self._prettify(self.raw_tree)

    def _prettify(self, raw_tree):
        self_branch = None
        if PhylotreeParser.SELF_BRANCH_NAME in raw_tree:
            self_branch = raw_tree[PhylotreeParser.SELF_BRANCH_NAME]
        descentants = {}
        for (name, branches) in raw_tree.items():
            if name != PhylotreeParser.SELF_BRANCH_NAME:
                descentants[name] = self._prettify(branches)

        pretty_tree = {}
        if self_branch:
            pretty_tree['conditions'] = self_branch['conditions']
            pretty_tree['example_accessions'] = self_branch['example_accessions']

        if len(descentants) != 0:
            pretty_tree['descentants'] = descentants
            
        return pretty_tree

    def prettify_to_array(self):
        return self._prettify_to_array(self.prettify())

    def _prettify_to_array(self, raw_tree, pretty_tree=[], parent=[]):
        for name, branches in raw_tree['descentants'].items():
            current = [*parent, { 'name': name, 'conditions': branches['conditions'] }]
            pretty_tree.append(current)
            if branches['descentants']:
                self._prettify_to_array(branches, pretty_tree, parent)

        return pretty_tree
        
    def trim_white_space(self, text):
        return Strings(text).replace('^\s+', '').replace('\s+$', '').replace('\s+', ' ').value

    def is_in_useful_range(self, text):
        if not self.in_useful_range and self.Pattern.is_table_title(text):
            self.in_useful_range = True
            return False
        return self.in_useful_range

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--input', '-i', required=True,
                           help='Path to input file')
    args = argparser.parse_args()
    parser = PhylotreeParser(args.input)
    parser.parse()
    pretty = parser.prettify()
    print(json.dumps(pretty, sort_keys=True, indent=4))
